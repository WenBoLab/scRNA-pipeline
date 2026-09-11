from __future__ import annotations

import json
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy import sparse

from .common import save_adata, sha256, validate_counts, write_json

DEMO_URL = "https://falexwolf.de/data/pbmc3k_raw.h5ad"
DATASET_PAGE = "https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0"


def download_demo(path):
    """Stage verified, unnormalized PBMC3k counts, preferably from the bundled copy."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    checksum_file = Path(__file__).with_name("pbmc3k.sha256")
    expected = checksum_file.read_text().strip().split()[0]
    resource_dir = Path(__file__).resolve().parents[2] / "resources"
    bundle = resource_dir / "pbmc3k_counts.h5ad"
    metadata_path = resource_dir / "pbmc3k_counts.source.json"
    metadata = json.loads(metadata_path.read_text()) if metadata_path.exists() else {}
    allowed = {expected}
    if metadata.get("sha256"):
        allowed.add(metadata["sha256"])
    if path.exists():
        if sha256(path) not in allowed:
            raise ValueError(f"Checksum mismatch in {path}. Remove this file and retry.")
    else:
        tmp = path.with_suffix(".download")
        try:
            if bundle.exists():
                if sha256(bundle) != metadata.get("sha256"):
                    raise ValueError("Bundled PBMC3k checksum mismatch; restore resources from the repository")
                shutil.copyfile(bundle, tmp)
                selected_checksum = metadata["sha256"]
            else:
                with urllib.request.urlopen(DEMO_URL, timeout=120) as response, tmp.open("wb") as out:
                    while chunk := response.read(1024 * 1024):
                        out.write(chunk)
                selected_checksum = expected
            if sha256(tmp) != selected_checksum:
                raise ValueError("PBMC3k download checksum mismatch; do not use this file")
            tmp.replace(path)
        finally:
            tmp.unlink(missing_ok=True)
    actual_checksum = sha256(path)
    write_json(path.with_suffix(".source.json"), {
        "dataset": "PBMC3k", "provider": "10x Genomics", "dataset_page": DATASET_PAGE,
        "download_url": DEMO_URL, "mirror": "Scanpy pbmc3k raw-count mirror",
        "sha256": actual_checksum, "original_download_sha256": expected,
        "packaging": metadata.get("packaging") if actual_checksum == metadata.get("sha256") else "Original Scanpy mirror file",
        "license": "CC BY 4.0",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    })
    return path


def read_counts(path, fmt):
    if fmt == "10x_mtx":
        a = sc.read_10x_mtx(path, var_names="gene_symbols", make_unique=True, gex_only=True)
    elif fmt == "10x_h5":
        a = sc.read_10x_h5(path, gex_only=True)
    elif fmt == "h5ad":
        a = sc.read_h5ad(path)
    else:
        raise ValueError(f"Unknown input format: {fmt}")
    validate_counts(a)
    a.var["gene_symbol"] = a.var_names.astype(str)
    a.var_names_make_unique()
    # Imported annotations/embeddings may originate from another analysis: recompute them.
    clean = ad.AnnData(sparse.csr_matrix(a.X, dtype="float32"), obs=a.obs.iloc[:, :0].copy(), var=a.var.copy())
    return clean


def ingest(cfg):
    if cfg["input"]["mode"] == "demo":
        path = download_demo(cfg["input"]["path"])
        samples = pd.DataFrame([dict(sample_id="pbmc3k", path=str(path), format="h5ad",
                                     donor="donor1", condition="healthy", batch="batch1")])
    else:
        samples = pd.read_csv(cfg["input"]["samples"], sep="\t", comment="#", dtype=str)
    required = {"sample_id", "path", "format", "donor", "condition", "batch"}
    if not required.issubset(samples.columns) or samples.empty:
        raise ValueError("samples.tsv requires at least one sample and: " + ", ".join(sorted(required)))
    if samples[list(required)].isna().any().any() or (samples[list(required)].eq("").any().any()):
        raise ValueError("Sample metadata cannot be empty")
    if not samples.sample_id.is_unique:
        raise ValueError("sample_id must be unique (one capture/library per row)")
    items, audit = {}, []
    for row in samples.itertuples(index=False):
        a = read_counts(row.path, row.format)
        a.obs["barcode"] = a.obs_names.astype(str)
        for column in ["sample_id", "donor", "condition", "batch"]:
            a.obs[column] = str(getattr(row, column))
        # Hash matrix inputs for later auditing, including MTX feature and barcode files.
        p = Path(row.path)
        files = [p] if p.is_file() else sorted(f for f in p.iterdir() if f.is_file())
        audit.append(dict(sample_id=row.sample_id, n_cells=a.n_obs, n_genes=a.n_vars,
                          files={str(f): sha256(f) for f in files}))
        items[row.sample_id] = a
    common_genes = set.intersection(*(set(a.var_names) for a in items.values()))
    if len(common_genes) < 3:
        raise ValueError("Samples do not share sufficient gene identifiers; harmonize annotation first")
    a = ad.concat(items, join="inner", index_unique="::", merge="first")
    if not a.obs_names.is_unique:
        raise ValueError("Sample/barcode combination is not unique")
    a.uns["input_provenance"] = json.dumps(audit)
    a.uns["feature_join"] = "inner; only genes measured in all samples are retained"
    root = Path(cfg["output_dir"])
    provenance = {"samples": audit, "shared_genes": a.n_vars}
    if cfg["input"]["mode"] == "demo":
        provenance["dataset_source"] = json.loads(path.with_suffix(".source.json").read_text())
    write_json(root / "provenance/inputs.json", provenance)
    save_adata(a, root / "objects/01_counts.h5ad")
