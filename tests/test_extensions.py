import argparse
import gzip
import importlib.util
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.io import mmwrite

from scrna_learn.data import read_counts
from scrna_learn.extensions import pseudobulk_de


@pytest.mark.parametrize("fmt", ["10x_mtx", "10x_h5"])
def test_tenx_import_preserves_orientation_and_gene_expression(fmt, tmp_path):
    # 4 cells x 3 genes; on disk 10x is features x barcodes.
    counts = np.array([[1, 2, 3], [4, 1, 2], [3, 4, 1], [2, 3, 4]])
    feature_matrix = sparse.csc_matrix(counts.T)
    genes = ["MT-CO1", "CD3D", "MS4A1"]
    cells = [f"bc{i}" for i in range(4)]
    if fmt == "10x_mtx":
        with gzip.open(tmp_path / "matrix.mtx.gz", "wb") as f:
            mmwrite(f, feature_matrix)
        with gzip.open(tmp_path / "features.tsv.gz", "wt") as f:
            f.write("".join(f"id{i}\t{g}\tGene Expression\n" for i, g in enumerate(genes)))
        with gzip.open(tmp_path / "barcodes.tsv.gz", "wt") as f:
            f.write("\n".join(cells) + "\n")
        path = tmp_path
    else:
        path = tmp_path / "matrix.h5"
        with h5py.File(path, "w") as h:
            g = h.create_group("matrix")
            for name in ["data", "indices", "indptr"]:
                g[name] = getattr(feature_matrix, name)
            g["shape"] = feature_matrix.shape
            g["barcodes"] = np.array(cells, dtype="S")
            feat = g.create_group("features")
            feat["id"] = np.array([f"id{i}" for i in range(3)], dtype="S")
            feat["name"] = np.array(genes, dtype="S")
            feat["feature_type"] = np.array(["Gene Expression"] * 3, dtype="S")
            feat["genome"] = np.array(["GRCh38"] * 3, dtype="S")
    a = read_counts(path, fmt)
    np.testing.assert_array_equal(a.X.toarray(), counts)
    assert list(a.var_names) == genes
    assert list(a.obs_names) == cells


def test_pydeseq2_recovers_large_known_synthetic_effect(tmp_path):
    rng = np.random.default_rng(71)
    pieces, metadata = [], []
    for donor in range(8):
        case = donor >= 4
        means = rng.uniform(1, 2, 120)
        means[:10] *= 6 if case else 1
        x = rng.poisson(means, (25, 120))
        pieces.append(x)
        metadata.extend([dict(donor=f"synthetic_{donor}", condition="case" if case else "control",
                              cell_type="synthetic_T", annotation_status="user_reviewed") for _ in range(25)])
    a = ad.AnnData(sparse.csr_matrix(np.vstack(pieces)), obs=pd.DataFrame(metadata, index=[f"c{i}" for i in range(200)]),
                   var=pd.DataFrame(index=[f"g{i}" for i in range(120)]))
    a.layers["counts"] = a.X.copy()
    path = tmp_path / "synthetic.h5ad"; a.write_h5ad(path)
    pseudobulk_de(path, tmp_path / "de", "synthetic_T", "case", "control")
    result = pd.read_csv(tmp_path / "de/differential_expression.csv", index_col="gene")
    assert result.loc[[f"g{i}" for i in range(10)], "log2FoldChange"].median() > 1.5
    assert (result.loc[[f"g{i}" for i in range(10)], "padj"] < .05).sum() >= 8
    assert pd.read_csv(tmp_path / "de/design_metadata.csv").shape[0] == 8


def test_cellranger_command_keeps_paths_as_literal_arguments():
    spec = importlib.util.spec_from_file_location("cr", "scripts/cellranger_count.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    args = argparse.Namespace(sample_id="test_A", fastqs="test path/fastqs", fastq_sample="sample A",
                              transcriptome="test path/ref", create_bam=False, cores=8, memory_gb=64)
    cmd = module.build_command(args)
    assert cmd[0:2] == ["cellranger", "count"]
    assert "--sample=sample A" in cmd
    assert "--create-bam=false" in cmd
    args.sample_id = "../unsafe"
    with pytest.raises(ValueError):
        module.build_command(args)


def test_harmony_branch_on_explicitly_synthetic_batches(tmp_path):
    spec = importlib.util.spec_from_file_location("fixture", "scripts/make_test_data.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    cfg = module.make_fixture(tmp_path)
    original = ad.read_h5ad(tmp_path / "synthetic.h5ad")
    rows = []
    for i in range(2):
        path = tmp_path / f"synthetic_batch{i}.h5ad"
        original[i::2].copy().write_h5ad(path)
        rows.append(dict(sample_id=f"test{i}", path=str(path), format="h5ad",
                         donor=f"synthetic{i}", condition="synthetic", batch=f"batch{i}"))
    pd.DataFrame(rows).to_csv(tmp_path / "samples.tsv", sep="\t", index=False)
    cfg["embedding"]["integration"] = "harmony"
    from scrna_learn.cli import run_stage
    for stage in ["ingest", "qc", "normalize", "embed"]:
        run_stage(stage, cfg)
    root = Path(cfg["output_dir"])
    result = ad.read_h5ad(root / "objects/04_clustered.h5ad")
    assert result.obs.sample_id.nunique() == 2
    assert result.obsm["X_pca_harmony"].shape == result.obsm["X_pca"].shape
    assert np.isfinite(result.obsm["X_pca_harmony"]).all()
    assert result.uns["neighbors"]["params"]["use_rep"] == "X_pca_harmony"
    before = ad.read_h5ad(root / "objects/03_normalized.h5ad")
    assert (before.layers["counts"] != result.layers["counts"]).nnz == 0
