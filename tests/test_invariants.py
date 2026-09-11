import importlib.util
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from scrna_learn.common import validate_counts
from scrna_learn.extensions import aggregate_pseudobulk, check_de_design, ora


@pytest.mark.parametrize("bad", [np.nan, -1., .25, np.inf])
def test_rejects_invalid_counts(bad):
    x = np.ones((3, 3)); x[0, 0] = bad
    with pytest.raises(ValueError):
        validate_counts(ad.AnnData(sparse.csr_matrix(x)))


def test_pseudobulk_preserves_counts_and_pools_technical_replicates():
    x = np.arange(18).reshape(6, 3)
    obs = pd.DataFrame({"donor": ["d1"] * 3 + ["d2"] * 3, "condition": ["a"] * 6,
                       "cell_type": ["T"] * 6, "sample_id": ["lib1", "lib1", "lib2", "lib3", "lib3", "lib3"]},
                       index=[f"c{i}" for i in range(6)])
    a = ad.AnnData(sparse.csr_matrix(x), obs=obs)
    a.layers["counts"] = a.X.copy()
    a.X = sparse.csr_matrix(np.log1p(x))
    pb, audit = aggregate_pseudobulk(a, min_cells=2)
    assert pb.n_obs == 2
    np.testing.assert_array_equal(pb.X.toarray(), np.stack([x[:3].sum(0), x[3:].sum(0)]))
    assert int(pb.X.sum()) == int(x.sum())
    assert audit.retained.all()


def test_de_refuses_pseudoreplicates_and_requires_pairing():
    one = pd.DataFrame({"condition": ["case"] * 10 + ["control"] * 10,
                        "donor": ["d1"] * 10 + ["d2"] * 10})
    with pytest.raises(ValueError, match="3 independent"):
        check_de_design(one, "case", "control", False)
    paired = pd.DataFrame({"donor": ["d1", "d2", "d3"] * 2,
                           "condition": ["case"] * 3 + ["control"] * 3})
    with pytest.raises(ValueError, match="paired"):
        check_de_design(paired, "case", "control", False)
    check_de_design(paired, "case", "control", True)


def test_enrichment_uses_measured_background_and_corrects_zero_hit_terms():
    sets = {"hit": ("test", {"a", "b", "unmeasured"}), "zero": ("test", {"c", "d"})}
    result = ora({"a", "b", "unmeasured"}, {"a", "b", "c", "d"}, sets, 1, 10).set_index("term")
    assert result.loc["hit", "universe_size"] == 4
    assert result.loc["hit", "query_size"] == 2
    assert result.loc["hit", "set_size"] == 2
    assert result.loc["hit", "pvalue"] == pytest.approx(1 / 6)
    assert result.loc["hit", "padj"] == pytest.approx(1 / 3)
    assert result.loc["zero", "pvalue"] == 1


def test_offline_end_to_end_preserves_umi_and_exports_results(tmp_path):
    spec = importlib.util.spec_from_file_location("fixture", "scripts/make_test_data.py")
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    cfg = module.make_fixture(tmp_path)
    from scrna_learn.cli import run_stage, STAGES
    for stage in STAGES:
        run_stage(stage, cfg)
    root = Path(cfg["output_dir"])
    original = ad.read_h5ad(root / "objects/01_counts.h5ad")
    final = ad.read_h5ad(root / "objects/06_annotated.h5ad")
    assert final.n_obs == 210
    assert final.obs.cluster.nunique() >= 4
    assert final.var.highly_variable.sum() < final.n_vars
    expected = original[final.obs_names, final.var_names].X
    np.testing.assert_array_equal(final.layers["counts"].toarray(), expected.toarray())
    assert np.isfinite(final.obsm["X_umap"]).all()
    assert not np.allclose(final.X.toarray(), expected.toarray())
    markers = pd.read_csv(root / "tables/markers_all.csv")
    assert markers["names"].nunique() == final.n_vars  # markers use all genes, not just HVGs
    assert (root / "report.html").stat().st_size > 10000
    assert (root / "provenance/run_manifest.json").is_file()

