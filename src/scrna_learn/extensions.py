from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import hypergeom

from .common import outpaths, save_adata, sha256, validate_counts, write_json


def adjust_bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    result = np.empty_like(p)
    result[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return result


def read_gmt(path):
    sets = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, source, *genes = line.split("\t")
        if not genes or name in sets:
            raise ValueError("GMT needs unique set names and >=1 gene per set")
        sets[name] = (source, {g.strip() for g in genes if g.strip()})
    if not sets:
        raise ValueError("Empty GMT file")
    return sets


def ora(query, universe, gene_sets, min_genes=5, max_genes=500):
    """Hypergeometric ORA with all eligible terms included in BH correction."""
    universe = set(universe)
    query = set(query) & universe
    rows = []
    for name, (source, genes) in gene_sets.items():
        measured = set(genes) & universe
        if not min_genes <= len(measured) <= max_genes:
            continue
        hits = query & measured
        p = float(hypergeom.sf(len(hits) - 1, len(universe), len(measured), len(query))) if query else 1.0
        rows.append(dict(term=name, source=source, overlap=len(hits), set_size=len(measured),
                         query_size=len(query), universe_size=len(universe), pvalue=p,
                         overlap_genes=";".join(sorted(hits))))
    result = pd.DataFrame(rows, columns=["term", "source", "overlap", "set_size", "query_size", "universe_size", "pvalue", "overlap_genes"])
    result["padj"] = adjust_bh(result.pvalue)
    return result.sort_values("padj")


def enrichment(cfg):
    root = outpaths(cfg)
    path = cfg["enrichment"]["gmt"]
    if not path:
        pd.DataFrame(columns=["cluster", "term", "padj"]).to_csv(root / "tables/enrichment.csv", index=False)
        pd.DataFrame(columns=["cluster", "term", "padj"]).to_csv(root / "tables/enrichment_top.csv", index=False)
        (root / "figures/enrichment.png").unlink(missing_ok=True)
        write_json(root / "tables/enrichment_status.json", {
            "status": "skipped", "reason": "未提供有来源和版本的 GMT 基因集。默认演示不生成通路显著性结论。"})
        return
    a = sc.read_h5ad(root / "objects/05_markers.h5ad")
    markers = pd.read_csv(root / "tables/markers_positive.csv", dtype={"group": str})
    sets = read_gmt(path)
    tables = []
    for cluster in a.obs.cluster.cat.categories:
        query = markers[markers.group == str(cluster)].head(cfg["enrichment"]["top_n"])["names"]
        result = ora(query, a.var_names, sets, cfg["enrichment"]["min_genes"], cfg["enrichment"]["max_genes"])
        result.insert(0, "cluster", str(cluster))
        tables.append(result)
    combined = pd.concat(tables, ignore_index=True)
    combined.to_csv(root / "tables/enrichment.csv", index=False)
    top = combined[combined.padj < .05].groupby("cluster", observed=True).head(5)
    top.to_csv(root / "tables/enrichment_top.csv", index=False)
    import matplotlib.pyplot as plt
    if len(top):
        fig, ax = plt.subplots(figsize=(10, max(4, len(top) * .28)))
        labels = [f"C{r.cluster} | {r.term}" for r in top.itertuples()]
        y = np.arange(len(top))
        ax.barh(y, -np.log10(np.maximum(top.padj.to_numpy(), 1e-300)), color="#278d84")
        ax.set_yticks(y, [s if len(s) <= 68 else s[:65] + "..." for s in labels], fontsize=7)
        ax.invert_yaxis()
        ax.set(xlabel="-log10(BH adjusted p)", title="Exploratory cluster-marker ORA (top 5 / cluster)")
        fig.tight_layout(); fig.savefig(root / "figures/enrichment.png", bbox_inches="tight", dpi=160); plt.close(fig)
    else:
        (root / "figures/enrichment.png").unlink(missing_ok=True)
    write_json(root / "tables/enrichment_status.json", {
        "status": "completed", "gmt": str(path), "sha256": sha256(path), "tested_rows": len(combined),
        "significant_rows": int((combined.padj < .05).sum()),
        "universe": "all tested QC-passing genes", "multiple_testing": "BH across all eligible terms within each cluster",
        "interpretation": "Exploratory cluster-marker enrichment; not a condition-level disease test"})


def aggregate_pseudobulk(a, min_cells=20):
    """Sum raw UMI counts per donor x condition x cell type, pooling technical captures."""
    required = ["donor", "condition", "cell_type"]
    if not set(required).issubset(a.obs) or "counts" not in a.layers:
        raise ValueError("Need donor/condition/cell_type metadata and layers['counts']")
    if a.obs[required].isna().any().any():
        raise ValueError("Pseudobulk metadata cannot be missing")
    counts = sparse.csr_matrix(a.layers["counts"])
    checked = ad.AnnData(counts)
    validate_counts(checked)
    data, rows, audit = [], [], []
    for (donor, condition, cell_type), frame in a.obs.groupby(required, observed=True, sort=True):
        positions = a.obs_names.get_indexer(frame.index)
        retained = len(frame) >= min_cells
        audit.append(dict(donor=donor, condition=condition, cell_type=cell_type, n_cells=len(frame), retained=retained))
        if not retained:
            continue
        data.append(sparse.csr_matrix(counts[positions].sum(axis=0), dtype=np.int64))
        rows.append(dict(donor=str(donor), condition=str(condition), cell_type=str(cell_type),
                         n_cells=len(frame), sample_id=f"pb_{len(rows):04d}"))
    if not rows:
        raise ValueError("No donor-condition-celltype group has enough cells")
    obs = pd.DataFrame(rows).set_index("sample_id")
    result = ad.AnnData(sparse.vstack(data, format="csr"), obs=obs, var=a.var.copy())
    result.uns["aggregation"] = "sum raw counts by donor x condition x cell_type; technical captures pooled"
    return result, pd.DataFrame(audit)


def check_de_design(metadata, case, control, paired):
    if case == control or set(metadata.condition) != {case, control}:
        raise ValueError("Contrast must contain exactly the requested case and control")
    counts = metadata.groupby("condition", observed=True).donor.nunique()
    if (counts < 3).any():
        raise ValueError("At least 3 independent donors per condition are required by this teaching workflow")
    overlap = set(metadata.loc[metadata.condition == case, "donor"]) & set(metadata.loc[metadata.condition == control, "donor"])
    if paired:
        if len(overlap) != metadata.donor.nunique():
            raise ValueError("Paired mode requires every donor in both conditions")
        if metadata.duplicated(["donor", "condition"]).any():
            raise ValueError("Aggregate technical replicates before paired analysis")
    elif overlap:
        raise ValueError("Donors occur in both conditions; use --paired for paired design")


def pseudobulk_de(input_path, output_dir, cell_type, case, control, min_cells=20, paired=False):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    a = sc.read_h5ad(input_path)
    if "annotation_status" not in a.obs or not a.obs.annotation_status.eq("user_reviewed").all():
        raise ValueError("Review cell annotations first; pseudobulk DE requires annotation_status=user_reviewed")
    pb, audit = aggregate_pseudobulk(a, min_cells=min_cells)
    selected = pb[(pb.obs.cell_type == cell_type) & pb.obs.condition.isin([case, control])].copy()
    if selected.n_obs == 0:
        raise ValueError("No matching cell type and conditions")
    check_de_design(selected.obs, case, control, paired)
    counts = pd.DataFrame(selected.X.toarray(), index=selected.obs_names, columns=selected.var_names)
    expressed = (counts >= 10).sum(axis=0) >= 3
    counts = counts.loc[:, expressed]
    if counts.shape[1] < 20:
        raise ValueError("Too few expressed genes for DESeq2")
    design = "~ donor + condition" if paired else "~ condition"
    dds = DeseqDataSet(counts=counts, metadata=selected.obs.copy(), design=design,
                       refit_cooks=True, n_cpus=1)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", case, control], n_cpus=1)
    stats.summary()
    root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
    save_adata(pb, root / "pseudobulk_counts.h5ad")
    audit.to_csv(root / "aggregation_audit.csv", index=False)
    selected.obs.to_csv(root / "design_metadata.csv")
    stats.results_df.sort_values("padj").to_csv(root / "differential_expression.csv", index_label="gene")
    write_json(root / "design.json", {"design": design, "cell_type": cell_type, "case": case, "control": control,
               "positive_log2FC": f"higher in {case}", "min_cells": min_cells, "input_sha256": sha256(input_path),
               "unit": "independent donor, not individual cell", "limitations": "Simple balanced design only; model additional covariates explicitly for real studies."})
