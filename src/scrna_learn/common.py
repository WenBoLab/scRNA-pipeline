from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str) + "\n")


def load_config(path):
    cfg = yaml.safe_load(Path(path).read_text())
    return validate_config(cfg)


def validate_config(cfg):
    if cfg["input"]["mode"] not in {"demo", "samples"}:
        raise ValueError("input.mode must be demo or samples")
    q = cfg["qc"]
    if not 0 <= q["min_genes"] < q["max_genes"] or not 0 < q["max_pct_mt"] <= 100:
        raise ValueError("Invalid QC bounds")
    if cfg["embedding"]["integration"] not in {"none", "harmony"}:
        raise ValueError("embedding.integration must be none or harmony")
    if not 0 < cfg["doublets"]["expected_rate"] < 1:
        raise ValueError("doublets.expected_rate must be between zero and one")
    if cfg["normalization"]["target_sum"] <= 0 or cfg["normalization"]["n_hvg"] < 3:
        raise ValueError("Invalid normalization settings")
    return cfg


def validate_counts(adata):
    """Reject transformed, negative or invalid inputs before any analysis."""
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError("Need at least 3 cells and 3 genes")
    vals = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(vals).all() or (vals < 0).any():
        raise ValueError("Counts must be finite and nonnegative")
    if not np.allclose(vals, np.rint(vals), atol=1e-6, rtol=0):
        raise ValueError("Input is not integer UMI counts; do not use normalized h5ad")
    if not adata.obs_names.is_unique:
        raise ValueError("Duplicated cell barcodes within a sample")
    return adata


def save_adata(adata, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.h5ad")
    adata.write_h5ad(temporary, compression="gzip")
    temporary.replace(path)


def outpaths(cfg):
    root = Path(cfg["output_dir"])
    for sub in ["objects", "tables", "figures", "provenance"]:
        (root / sub).mkdir(parents=True, exist_ok=True)
    return root
