"""The stores tests/unit/preview-store-parity.test.ts compares against are
what scripts/prepare.py writes today.

They are committed (the unit tests run without Python), so a change to
prepare.py could leave them describing an older pipeline. Each is written
afresh here and compared with the committed one by content — the table's
cells, genes, values, columns, positions and colours, and mapping.json —
not by bytes, which differ between zarr and blosc versions. When this
fails, regenerate them:

  python scripts/testdata/make-tiny-parity-stores.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

pytest.importorskip("spatialdata")
ad = pytest.importorskip("anndata")
zarr = pytest.importorskip("zarr")

spec = importlib.util.spec_from_file_location(
    "make_tiny_parity_stores", ROOT / "scripts" / "testdata" / "make-tiny-parity-stores.py"
)
stores = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stores)


def table(store: Path):
    return ad.io.read_elem(zarr.open_group(store / "tables" / "table", mode="r"))


def dense(x) -> np.ndarray:
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x)


@pytest.mark.parametrize("case", sorted(stores.CASES))
def test_the_committed_store_is_what_prepare_writes_today(case, tmp_path):
    fresh = tmp_path / case
    stores.build(case, stores.CASES[case], fresh)
    kept = stores.OUT / case
    regenerate = "regenerate: python scripts/testdata/make-tiny-parity-stores.py"

    assert json.loads((fresh / "mapping.json").read_text()) == json.loads(
        (kept / "mapping.json").read_text()
    ), regenerate

    a, b = table(fresh / f"{case}.zarr"), table(kept / f"{case}.zarr")
    assert list(a.obs_names) == list(b.obs_names), regenerate
    assert list(a.var_names) == list(b.var_names), regenerate
    np.testing.assert_allclose(dense(a.X), dense(b.X), err_msg=regenerate)
    assert list(a.obs.columns) == list(b.obs.columns), regenerate
    for column in a.obs.columns:
        assert a.obs[column].astype(str).tolist() == b.obs[column].astype(str).tolist(), (
            f"obs/{column} differs; {regenerate}"
        )
    assert sorted(a.obsm) == sorted(b.obsm), regenerate
    for key in a.obsm:
        np.testing.assert_allclose(
            np.asarray(a.obsm[key], dtype=float),
            np.asarray(b.obsm[key], dtype=float),
            err_msg=f"obsm/{key} differs; {regenerate}",
        )
    colours = lambda t: {k: list(v) for k, v in t.uns.items() if k.endswith("_colors")}  # noqa: E731
    assert colours(a) == colours(b), regenerate
