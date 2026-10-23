"""scripts/de_stats.py: the DEG panel's numbers, kept in the table."""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import de_stats  # noqa: E402


def table(tmp_path, adata):
    root = zarr.open_group(tmp_path / "t.zarr", mode="w", zarr_format=3)
    ad.io.write_elem(root, "table", adata)
    return root["table"]


def make(n=60):
    rng = np.random.default_rng(3)
    x = rng.poisson(1.2, (n, 5)).astype(np.float32)
    x[7] = 0  # a cell with no counts
    obs = pd.DataFrame(
        {
            "celltype": pd.Categorical(["a", "b", "c"] * (n // 3)),
            "one": pd.Categorical(["only"] * n),
            "id": pd.Categorical([f"c{i}" for i in range(n)]),
            "area": np.linspace(1, 2, n),
        },
        index=[str(i) for i in range(n)],
    )
    return ad.AnnData(X=sp.csr_matrix(x), obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(5)]))


def test_the_numbers(tmp_path):
    adata = make()
    group = table(tmp_path, adata)
    assert de_stats.write_de_stats(group, lambda *_: None) == ["celltype", "id"]

    written = ad.read_zarr(tmp_path / "t.zarr" / "table")
    mean, pct = written.varm["de_mean__celltype"], written.varm["de_pct__celltype"]
    record = written.uns["merfisheyes_de"]["celltype"]
    assert list(record["categories"]) == ["a", "b", "c"]
    assert list(record["cell_counts"]) == [20, 20, 20]

    # By hand, category "b"
    x = adata.X.toarray()
    totals = x.sum(1)
    target = np.median(totals[totals > 0])
    lognorm = np.log1p(np.divide(x * target, totals[:, None], out=np.zeros_like(x), where=totals[:, None] > 0))
    b = (adata.obs["celltype"] == "b").to_numpy()
    np.testing.assert_allclose(mean[:, 1], lognorm[b].mean(0), rtol=1e-5)
    np.testing.assert_allclose(pct[:, 1], (x[b] > 0).mean(0), rtol=1e-6)


def test_matches_scanpy_normalisation():
    sc = pytest.importorskip("scanpy")
    x = make().X
    reference = ad.AnnData(x.copy())
    sc.pp.normalize_total(reference)
    sc.pp.log1p(reference)
    np.testing.assert_allclose(de_stats.log_normalize(x).toarray(), reference.X.toarray(), atol=1e-6)


def test_which_columns(tmp_path):
    adata = make(n=1200)
    # Too many categories to be a grouping
    adata.obs["id"] = pd.Categorical([f"c{i}" for i in range(1200)])
    group = table(tmp_path, adata)
    group.attrs["instance_key"] = "celltype"  # a join key is not summarised either
    assert de_stats.categorical_columns(group) == []


def test_unassigned_cells_are_left_out(tmp_path):
    adata = make()
    adata.obs["celltype"] = adata.obs["celltype"].cat.add_categories([])
    adata.obs.loc[["0", "1"], "celltype"] = np.nan
    group = table(tmp_path, adata)
    de_stats.write_de_stats(group, lambda *_: None, columns=["celltype"])
    counts = ad.read_zarr(tmp_path / "t.zarr" / "table").uns["merfisheyes_de"]["celltype"]["cell_counts"]
    assert int(np.sum(counts)) == 58
