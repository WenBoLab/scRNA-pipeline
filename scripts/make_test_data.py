"""Create explicitly synthetic counts for offline CI, never biological example results."""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


def make_fixture(destination):
    root = Path(destination).resolve(); root.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path("config/config.yaml").read_text())
    panels = yaml.safe_load(Path(cfg["annotation"]["markers"]).read_text())
    panel_genes = list(dict.fromkeys(g for genes in panels.values() for g in genes))
    genes = ["MT-TEST"] + panel_genes + [f"GENE{i}" for i in range(160)]
    rng = np.random.default_rng(42)
    x = rng.poisson(.15, (210, len(genes))).astype("float32")
    for group, symbols in enumerate(panels.values()):
        inds = [genes.index(g) for g in symbols]
        x[np.ix_(np.arange(group * 30, (group + 1) * 30), inds)] += rng.poisson(10, (30, len(inds)))
    x[:, 0] = 0
    a = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame(index=[f"test{i}" for i in range(len(x))]),
                   var=pd.DataFrame(index=genes))
    a.uns["synthetic"] = "Software validation fixture; not real biological data"
    a.write_h5ad(root / "synthetic.h5ad")
    pd.DataFrame([dict(sample_id="synthetic", path=str(root / "synthetic.h5ad"), format="h5ad",
                       donor="synthetic_donor", condition="synthetic", batch="test")]).to_csv(root / "samples.tsv", sep="\t", index=False)
    cfg.update(project="SYNTHETIC_SOFTWARE_TEST_ONLY", output_dir=str(root / "output"))
    cfg["input"] = {"mode": "samples", "samples": str(root / "samples.tsv")}
    cfg["qc"].update(min_genes=5, max_genes=500, max_pct_mt=50)
    cfg["doublets"]["enabled"] = False
    cfg["enrichment"]["gmt"] = None
    cfg["normalization"]["n_hvg"] = 100
    cfg["embedding"].update(n_pcs=15, n_neighbors=10, resolution=.6, tsne=False)
    cfg["annotation"]["markers"] = str(Path("config/markers.yaml").resolve())
    (root / "config.yaml").write_text(yaml.safe_dump(cfg))
    return cfg


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--out", default="results/ci")
    make_fixture(p.parse_args().out)
