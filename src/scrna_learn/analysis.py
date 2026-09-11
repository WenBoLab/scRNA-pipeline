from __future__ import annotations

import copy
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse

from .common import outpaths, save_adata, write_json


def qc(cfg):
    root = outpaths(cfg)
    for old in (root / "figures").glob("doublets_*.png"):
        old.unlink()
    a = sc.read_h5ad(root / "objects/01_counts.h5ad")
    a.var["mt"] = a.var_names.str.startswith(cfg["qc"]["mitochondrial_prefix"])
    if not a.var["mt"].any():
        raise ValueError("No mitochondrial genes found. Check gene symbols and mitochondrial_prefix.")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    audited, kept, summaries = [], [], []
    for sample in a.obs["sample_id"].unique():
        s = a[a.obs["sample_id"] == sample].copy()
        q = {**cfg["qc"], **cfg["qc"].get("sample_overrides", {}).get(str(sample), {})}
        if not 0 <= q["min_genes"] < q["max_genes"] or not 0 < q["max_pct_mt"] <= 100:
            raise ValueError(f"Invalid QC override for {sample}")
        obs = s.obs.copy()
        obs["fail_low_genes"] = obs.n_genes_by_counts < q["min_genes"]
        obs["fail_high_genes"] = obs.n_genes_by_counts >= q["max_genes"]
        obs["fail_mt"] = obs.pct_counts_mt >= q["max_pct_mt"]
        obs["pass_basic_qc"] = ~obs[["fail_low_genes", "fail_high_genes", "fail_mt"]].any(axis=1)
        obs["doublet_score"] = np.nan
        obs["predicted_doublet"] = False
        obs["doublet_tested"] = False
        s = s[obs.pass_basic_qc].copy()
        if s.n_obs < 3:
            raise ValueError(f"{sample}: fewer than 3 cells after QC; inspect thresholds")
        threshold = None
        if cfg["doublets"]["enabled"]:
            if s.n_obs < 50:
                raise ValueError(f"{sample}: too few cells for teaching Scrublet workflow (<50)")
            # Scrublet must see unnormalized counts, independently for each capture.
            sc.pp.scrublet(s, expected_doublet_rate=cfg["doublets"]["expected_rate"],
                          threshold=cfg["doublets"]["threshold"], random_state=cfg["seed"],
                          n_prin_comps=min(30, s.n_obs - 2, s.n_vars - 2), verbose=False)
            if "predicted_doublet" not in s.obs or "threshold" not in s.uns["scrublet"]:
                raise ValueError("Scrublet did not determine a threshold; inspect scores and set doublets.threshold")
            threshold = float(s.uns["scrublet"]["threshold"])
            obs.loc[s.obs_names, "doublet_score"] = s.obs.doublet_score.to_numpy()
            obs.loc[s.obs_names, "predicted_doublet"] = s.obs.predicted_doublet.to_numpy()
            obs.loc[s.obs_names, "doublet_tested"] = True
            fig, ax = plt.subplots(figsize=(6, 3.4))
            ax.hist(s.obs.doublet_score, bins=40, alpha=.75, label="Observed cells")
            ax.hist(s.uns["scrublet"]["doublet_scores_sim"], bins=40, alpha=.45, label="Simulated doublets")
            ax.axvline(threshold, color="#af2538", linestyle="--", label=f"Threshold {threshold:.3f}")
            ax.set(xlabel="Scrublet score", ylabel="Count", title=f"Doublet detection: {sample}")
            ax.legend(fontsize=8)
            # Numeric index avoids using user sample IDs as filesystem paths.
            fig.tight_layout(); fig.savefig(root / f"figures/doublets_{len(summaries) + 1}.png", dpi=160); plt.close(fig)
        obs["retained"] = obs.pass_basic_qc & ~(obs.predicted_doublet & cfg["doublets"]["remove"])
        reasons = []
        for row in obs.itertuples():
            flags = [name for name, value in [("low_genes", row.fail_low_genes),
                     ("high_genes", row.fail_high_genes), ("high_mito", row.fail_mt),
                     ("predicted_doublet", row.predicted_doublet and cfg["doublets"]["remove"])] if value]
            reasons.append(";".join(flags) if flags else "retained")
        obs["decision"] = reasons
        retained = a[obs.index[obs.retained]].copy()
        retained.obs = obs.loc[retained.obs_names].copy()
        kept.append(retained)
        audited.append(obs)
        summaries.append(dict(sample_id=str(sample), before=len(obs), after_basic_qc=int(obs.pass_basic_qc.sum()),
                              predicted_doublets=int(obs.predicted_doublet.sum()), retained=int(obs.retained.sum()),
                              doublet_threshold=threshold, doublet_enabled=cfg["doublets"]["enabled"],
                              min_genes=q["min_genes"], max_genes=q["max_genes"], max_pct_mt=q["max_pct_mt"]))
    audit = pd.concat(audited)
    audit.to_csv(root / "tables/cell_qc.csv", index_label="cell_id")
    b = ad.concat(kept, merge="first", uns_merge="first")
    if b.n_obs < 10:
        raise ValueError("Fewer than 10 retained cells; inspect QC before proceeding")
    sc.pp.filter_genes(b, min_cells=cfg["qc"]["min_cells_per_gene"])
    if b.n_vars < 20:
        raise ValueError("Fewer than 20 retained genes; inspect input and QC")
    write_json(root / "tables/qc_summary.json", {"input_cells": a.n_obs, "retained_cells": b.n_obs,
              "input_genes": a.n_vars, "retained_genes": b.n_vars, "samples": summaries})
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    for ax, metric, title in zip(axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                                 ["Detected genes", "UMI counts", "Mitochondrial UMI (%)"]):
        values = [audit[metric].values, audit.loc[audit.retained, metric].values]
        ax.boxplot(values, tick_labels=["Before", "Retained"], showfliers=False)
        ax.set_title(title)
    fig.tight_layout(); fig.savefig(root / "figures/qc.png", dpi=170); plt.close(fig)
    save_adata(b, root / "objects/02_qc.h5ad")


def normalize(cfg):
    root = outpaths(cfg)
    a = sc.read_h5ad(root / "objects/02_qc.h5ad")
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=cfg["normalization"]["target_sum"])
    sc.pp.log1p(a)
    # 'seurat' dispersion HVGs expect log1p data; seurat_v3 would require raw counts.
    batch_key = "sample_id" if a.obs.sample_id.nunique() > 1 else None
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=min(cfg["normalization"]["n_hvg"], a.n_vars),
                                batch_key=batch_key, subset=False)
    a.var.to_csv(root / "tables/gene_metrics.csv", index_label="gene")
    sc.pl.highly_variable_genes(a, show=False)
    plt.gcf().savefig(root / "figures/highly_variable_genes.png", bbox_inches="tight", dpi=160); plt.close("all")
    a.uns["matrix_semantics"] = "X: log1p(CP10k or configured target_sum); layers/counts: raw integer UMI counts"
    save_adata(a, root / "objects/03_normalized.h5ad")


def embed(cfg):
    root = outpaths(cfg)
    a = sc.read_h5ad(root / "objects/03_normalized.h5ad")
    h = a[:, a.var.highly_variable].copy()
    if h.n_vars < 3:
        raise ValueError("Too few HVGs for PCA")
    sc.pp.scale(h, max_value=10)
    n_pcs = min(cfg["embedding"]["n_pcs"], h.n_obs - 1, h.n_vars - 1)
    sc.tl.pca(h, n_comps=n_pcs, svd_solver="arpack", random_state=cfg["seed"])
    a.obsm["X_pca"] = h.obsm["X_pca"].copy()
    a.uns["pca"] = copy.deepcopy(h.uns["pca"])
    a.varm["PCs"] = np.zeros((a.n_vars, n_pcs), dtype=np.float32)
    a.varm["PCs"][a.var.highly_variable.values] = h.varm["PCs"]
    del h
    representation = "X_pca"
    if cfg["embedding"]["integration"] == "harmony":
        key = cfg["embedding"]["batch_key"]
        if key not in a.obs or a.obs[key].nunique() < 2:
            raise ValueError("Harmony requires at least two real batches")
        cross = pd.crosstab(a.obs[key], a.obs.condition)
        if len(cross.columns) > 1 and ((cross > 0).sum(axis=1) == 1).all():
            raise ValueError("Batch is confounded with condition; Harmony cannot resolve this design")
        sc.external.pp.harmony_integrate(a, key=key, random_state=cfg["seed"],
                                          nclust=min(50, max(2, a.n_obs // 30)))
        representation = "X_pca_harmony"
    sc.pp.neighbors(a, n_neighbors=min(cfg["embedding"]["n_neighbors"], a.n_obs - 1),
                    n_pcs=n_pcs, use_rep=representation, random_state=cfg["seed"])
    sc.tl.leiden(a, resolution=cfg["embedding"]["resolution"], random_state=cfg["seed"],
                 flavor="igraph", n_iterations=2, directed=False, key_added="cluster")
    sc.tl.umap(a, random_state=cfg["seed"])
    if cfg["embedding"]["tsne"]:
        sc.tl.tsne(a, use_rep=representation, n_pcs=n_pcs, random_state=cfg["seed"],
                    perplexity=min(30, max(2, (a.n_obs - 1) / 3)), n_jobs=1)
    sc.pl.pca_variance_ratio(a, n_pcs=n_pcs, show=False)
    plt.gcf().savefig(root / "figures/pca_variance.png", bbox_inches="tight", dpi=160); plt.close("all")
    sc.pl.umap(a, color=["cluster", "pct_counts_mt"], wspace=.35, show=False)
    plt.gcf().savefig(root / "figures/umap_clusters.png", bbox_inches="tight", dpi=170); plt.close("all")
    if cfg["embedding"]["tsne"]:
        sc.pl.tsne(a, color="cluster", show=False)
        plt.gcf().savefig(root / "figures/tsne_clusters.png", bbox_inches="tight", dpi=170); plt.close("all")
    else:
        (root / "figures/tsne_clusters.png").unlink(missing_ok=True)
    if a.obs.sample_id.nunique() > 1:
        sc.pl.umap(a, color=["sample_id", "condition", "batch"], show=False)
        plt.gcf().savefig(root / "figures/umap_samples.png", bbox_inches="tight", dpi=160); plt.close("all")
    else:
        (root / "figures/umap_samples.png").unlink(missing_ok=True)
    a.obs.groupby(["sample_id", "cluster"], observed=True).size().rename("n_cells").to_csv(root / "tables/cluster_counts.csv")
    save_adata(a, root / "objects/04_clustered.h5ad")


def markers(cfg):
    root = outpaths(cfg)
    a = sc.read_h5ad(root / "objects/04_clustered.h5ad")
    sizes = a.obs.cluster.value_counts()
    if len(sizes) < 2 or sizes.min() < 2:
        raise ValueError("Marker comparison needs >=2 clusters, with >=2 cells in each")
    # Test all QC-passing genes on unscaled log-normalized data, not only HVGs.
    sc.tl.rank_genes_groups(a, groupby="cluster", method="wilcoxon", use_raw=False,
                           pts=True, tie_correct=True, corr_method="benjamini-hochberg")
    table = sc.get.rank_genes_groups_df(a, group=None)
    table.to_csv(root / "tables/markers_all.csv", index=False)
    positive = table[(table.logfoldchanges >= cfg["markers"]["min_log2fc"]) &
                     (table.pvals_adj < cfg["markers"]["max_adjusted_p"])]
    positive.to_csv(root / "tables/markers_positive.csv", index=False)
    positive.groupby("group", observed=True).head(10).to_csv(root / "tables/markers_top10.csv", index=False)
    save_adata(a, root / "objects/05_markers.h5ad")


def annotate(cfg):
    root = outpaths(cfg)
    a = sc.read_h5ad(root / "objects/05_markers.h5ad")
    panels = yaml.safe_load(Path(cfg["annotation"]["markers"]).read_text())
    valid, columns = {}, {}
    for idx, (label, genes) in enumerate(panels.items()):
        present = list(dict.fromkeys(g for g in genes if g in a.var_names))
        if len(present) >= 2:
            key = f"marker_score_{idx}"
            sc.tl.score_genes(a, present, score_name=key, random_state=cfg["seed"], use_raw=False)
            valid[label], columns[label] = present, key
    if len(valid) < 2:
        raise ValueError("Fewer than two usable marker panels; provide organism/tissue-appropriate markers")
    means = a.obs.groupby("cluster", observed=True)[list(columns.values())].mean()
    means.columns = list(columns)
    means.to_csv(root / "tables/annotation_scores.csv", index_label="cluster")
    rows, label_map = [], {}
    for cluster, row in means.iterrows():
        ranked = row.sort_values(ascending=False)
        best, score = ranked.index[0], float(ranked.iloc[0])
        margin = score - float(ranked.iloc[1])
        confident = score >= cfg["annotation"]["minimum_score"] and margin >= cfg["annotation"]["minimum_margin"]
        candidate = best if confident else "Uncertain"
        # Cytotoxic T cells also express NKG7/PRF1; an NK module alone is insufficient.
        t_genes = [g for g in ["CD3D", "CD3E"] if g in a.var_names]
        t_fraction = float((a[a.obs.cluster == cluster, t_genes].X > 0).mean()) if t_genes else 0.0
        if candidate == "NK cells" and t_fraction > 0.5:
            candidate = "T/NK unresolved"
        label_map[str(cluster)] = candidate
        rows.append(dict(cluster=str(cluster), n_cells=int((a.obs.cluster == cluster).sum()),
                         candidate_cell_type=candidate, top_panel=best, runner_up=ranked.index[1],
                         score=score, margin=margin, t_lineage_fraction=t_fraction, present_markers=";".join(valid[best]),
                         missing_markers=";".join(g for g in panels[best] if g not in a.var_names),
                         annotation_status="candidate_requires_review"))
    evidence = pd.DataFrame(rows)
    a.obs["candidate_cell_type"] = a.obs.cluster.astype(str).map(label_map).astype("category")
    a.obs["cell_type"] = a.obs.candidate_cell_type.copy()
    a.obs["annotation_status"] = "candidate_requires_review"
    manual_path = cfg["annotation"]["manual_labels"]
    if manual_path:
        manual = pd.read_csv(manual_path, sep="\t", dtype=str).fillna("")
        if not {"cluster", "cell_type", "evidence"}.issubset(manual.columns):
            raise ValueError("Manual annotation needs cluster, cell_type, evidence columns")
        if not manual.cluster.is_unique or set(manual.cluster) != set(a.obs.cluster.astype(str)):
            raise ValueError("Manual annotation must cover each current cluster exactly once")
        if (manual[["cell_type", "evidence"]].apply(lambda s: s.str.strip()) == "").any().any():
            raise ValueError("Manual labels and supporting evidence cannot be empty")
        a.obs["cell_type"] = a.obs.cluster.astype(str).map(manual.set_index("cluster").cell_type).astype("category")
        a.obs["annotation_status"] = "user_reviewed"
        evidence = evidence.merge(manual, on="cluster", validate="one_to_one")
        evidence["annotation_status"] = "user_reviewed"
    evidence.to_csv(root / "tables/annotation_evidence.csv", index=False)
    # Expression evidence is provided separately from marker-module scores.
    expression_rows = []
    for cluster in a.obs.cluster.cat.categories:
        sub = a[a.obs.cluster == cluster]
        for label, genes in valid.items():
            for gene in genes:
                x = sub[:, gene].X
                expression_rows.append(dict(cluster=str(cluster), panel=label, gene=gene,
                                             mean_log_expression=float(x.mean()),
                                             fraction_expressing=float((x > 0).mean())))
    pd.DataFrame(expression_rows).to_csv(root / "tables/annotation_expression.csv", index=False)
    sc.pl.dotplot(a, valid, groupby="cluster", use_raw=False, standard_scale="var", show=False)
    plt.gcf().savefig(root / "figures/marker_dotplot.png", bbox_inches="tight", dpi=170); plt.close("all")
    sc.pl.umap(a, color=["cluster", "cell_type"], wspace=.45, show=False,
               title=["Leiden clusters", "Cell types (candidate unless reviewed)"])
    plt.gcf().savefig(root / "figures/umap_celltypes.png", bbox_inches="tight", dpi=170); plt.close("all")
    feature_genes = [g for g in ["CD3D", "MS4A1", "NKG7", "LYZ", "FCER1A", "PPBP"] if g in a.var_names]
    if feature_genes:
        sc.pl.umap(a, color=feature_genes, ncols=3, use_raw=False, show=False)
        plt.gcf().savefig(root / "figures/feature_markers.png", bbox_inches="tight", dpi=160); plt.close("all")
    a.obs.to_csv(root / "tables/cell_metadata.csv", index_label="cell_id")
    counts = a.obs.groupby(["sample_id", "cell_type"], observed=True).size().rename("n_cells").reset_index()
    counts["fraction_within_sample"] = counts.n_cells / counts.groupby("sample_id", observed=True).n_cells.transform("sum")
    counts.to_csv(root / "tables/celltype_counts.csv", index=False)
    save_adata(a, root / "objects/06_annotated.h5ad")
