"""The teaching dataset must work when the upstream mirror is unavailable."""
import json
from pathlib import Path

import anndata as ad

from scrna_learn.common import sha256
from scrna_learn.data import download_demo


def test_bundled_counts_without_network(tmp_path, monkeypatch):
    def unavailable(*args, **kwargs):
        raise AssertionError("Default teaching data must not require a network request")

    monkeypatch.setattr("urllib.request.urlopen", unavailable)
    result = download_demo(tmp_path / "counts.h5ad")
    source = json.loads(result.with_suffix(".source.json").read_text())
    expected = json.loads((Path(__file__).parents[1] / "resources/pbmc3k_counts.source.json").read_text())
    assert sha256(result) == expected["sha256"] == source["sha256"]
    a = ad.read_h5ad(result)
    assert a.shape == (2700, 32738)
    assert int(a.X.sum()) == 6390631
    assert (a.X.data == a.X.data.round()).all()
    assert download_demo(result) == result
