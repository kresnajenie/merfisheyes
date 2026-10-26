"""scripts/transcripts_to_points.py, judged from outside: the element it writes
passes the validator, holds every molecule of every gene (blanks left out),
with the coordinates the table had, from a parquet file, a CSV, and a
SpatialData points element's parquet directory.

  python -m pytest tests/python/test_transcripts_to_points.py -q
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import points_zarr as pz  # noqa: E402
import transcripts_to_points as t2p  # noqa: E402
import validate_store  # noqa: E402

N = 6000


def table(seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    genes = [f"Gene{i}" for i in range(12)] + ["BLANK-3", "NegControlProbe_1"]
    return pd.DataFrame(
        {
            "feature_name": rng.choice(genes, N),
            "x_location": rng.uniform(0, 500, N).round(2),
            "y_location": rng.uniform(0, 400, N).round(2),
            "z_location": rng.uniform(0, 10, N).round(2),
        }
    )


def store_at(path: Path) -> Path:
    import zarr

    root = zarr.open_group(path, mode="w", zarr_format=3)
    root.attrs["spatialdata_attrs"] = {"version": "0.2"}
    return path


def check_element(store: Path, df: pd.DataFrame, gene_col: str, x_col: str, y_col: str, dims: int):
    r = validate_store.validate(store, sample=5000)
    assert not r.failed, [m for lv, _, m in r.lines if lv == "FAIL"]
    out = store / "points" / "molecules"
    _, attrs = pz.open_element(out)
    assert attrs["axes"] == ["x", "y", "z"][:dims]
    real = df[~df[gene_col].str.contains("BLANK|NegControl")]
    assert attrs["points:n_points"] == len(real)
    for gene in sorted(real[gene_col].unique())[:3]:
        got = pz.feature_points(out, gene, level=0)
        want = real[real[gene_col] == gene]
        assert len(got["x"]) == len(want)
        # stored as float32: equal to the table's values within a hundredth
        assert np.allclose(np.sort(got["x"]), np.sort(want[x_col].to_numpy()), atol=0.01)
        assert np.allclose(np.sort(got["y"]), np.sort(want[y_col].to_numpy()), atol=0.01)
    return attrs


def test_parquet_xenium_columns_recognised(tmp_path):
    df = table()
    df.to_parquet(tmp_path / "transcripts.parquet")
    store = store_at(tmp_path / "s.zarr")
    t2p.convert(tmp_path / "transcripts.parquet", store, log=lambda *_: None)
    attrs = check_element(store, df, "feature_name", "x_location", "y_location", 2)
    assert attrs["spatialdata_attrs"]["feature_key"] == "gene"


def test_csv_merscope_columns_and_z(tmp_path):
    df = table(1).rename(columns={"feature_name": "gene", "x_location": "global_x", "y_location": "global_y", "z_location": "global_z"})
    df.to_csv(tmp_path / "detected_transcripts.csv", index=False)
    store = store_at(tmp_path / "s.zarr")
    t2p.convert(tmp_path / "detected_transcripts.csv", store, keep_z=True, batch_size=1000, log=lambda *_: None)
    check_element(store, df, "gene", "global_x", "global_y", 3)


def test_custom_columns_and_an_unknown_header(tmp_path):
    df = table(2).rename(columns={"feature_name": "g", "x_location": "px", "y_location": "py"}).drop(columns=["z_location"])
    df.to_parquet(tmp_path / "t.parquet")
    store = store_at(tmp_path / "s.zarr")
    with pytest.raises(SystemExit, match="cannot tell the gene and coordinate columns"):
        t2p.convert(tmp_path / "t.parquet", store, log=lambda *_: None)
    t2p.convert(tmp_path / "t.parquet", store, columns="g,px,py", log=lambda *_: None)
    check_element(store, df, "g", "px", "py", 2)


def test_a_points_element_directory_with_its_feature_key(tmp_path):
    """spatialdata-io's layout: points/<name>/points.parquet/part.*.parquet, the
    feature key named in the element's zarr.json."""
    df = table(3).rename(columns={"feature_name": "feature_name", "x_location": "x", "y_location": "y", "z_location": "z"})
    store = store_at(tmp_path / "s.zarr")
    el = store / "points" / "transcripts"
    (el / "points.parquet").mkdir(parents=True)
    (store / "points" / "zarr.json").write_text(json.dumps({"zarr_format": 3, "node_type": "group"}))
    # Microns scaled into a pixel `global`, as a Xenium store declares it; the
    # element written must carry the same transform (the validator checks it)
    axes = [{"name": a, "type": "space"} for a in "xyz"]
    (el / "zarr.json").write_text(json.dumps({"zarr_format": 3, "node_type": "group", "attributes": {
        "encoding-type": "ngff:points", "axes": ["x", "y", "z"],
        "coordinateTransformations": [{"type": "scale", "scale": [4.7, 4.7, 1.0], "input": {"name": "xyz", "axes": axes},
                                       "output": {"name": "global", "axes": axes}}],
        "spatialdata_attrs": {"feature_key": "feature_name", "version": "0.2"}}}))
    half = len(df) // 2
    df.iloc[:half].to_parquet(el / "points.parquet" / "part.0.parquet")
    df.iloc[half:].to_parquet(el / "points.parquet" / "part.1.parquet")
    t2p.convert(el, store, transformations=[{"type": "scale", "scale": [4.7, 4.7], "output": {"name": "global", "axes": [
        {"name": "x", "type": "space"}, {"name": "y", "type": "space"}]}}], log=lambda *_: None)
    attrs = check_element(store, df, "feature_name", "x", "y", 2)
    assert attrs["coordinateTransformations"][0]["scale"] == [4.7, 4.7]


# ── assigned and unassigned (points-zarr-spec §8.1) ─────────────────────────


def unassigned_of(store: Path, gene: str):
    """A gene's `unassigned` values, read from feature_major as the viewer
    does, with that order's column names and copied columns."""
    g, _ = pz.open_element(store / "points" / "molecules")
    fm = g["0"]["feature_major"]
    codes = list(pz._col(fm, "gene")["enum_values"][:])
    f = codes.index(gene)
    lo, hi = (int(v) for v in g["0"]["index"]["feature_rows"][f : f + 2])
    return pz._col(fm, "unassigned")[lo:hi], list(fm["columns"][:]), fm.attrs["points:copied_columns"]


def test_merscope_cell_id_splits_the_molecules(tmp_path):
    df = table(4).rename(columns={"feature_name": "gene", "x_location": "global_x", "y_location": "global_y", "z_location": "global_z"})
    # MERSCOPE marks a molecule in no cell by -1
    df["cell_id"] = np.where(np.arange(N) % 3 == 0, -1, np.arange(N) + 100)
    df.to_parquet(tmp_path / "detected_transcripts.parquet")
    store = store_at(tmp_path / "s.zarr")
    t2p.convert(tmp_path / "detected_transcripts.parquet", store, batch_size=1000, log=lambda *_: None)
    check_element(store, df, "gene", "global_x", "global_y", 2)
    for gene in ("Gene0", "Gene5"):
        got, names, copied = unassigned_of(store, gene)
        assert "unassigned" in names and copied == ["unassigned"]
        assert len(got) == int((df.gene == gene).sum())
        assert int(got.sum()) == int(((df.gene == gene) & (df.cell_id == -1)).sum())


def test_a_named_cell_column_splits_by_the_shared_rule(tmp_path):
    df = table(5).rename(columns={"feature_name": "g", "x_location": "px", "y_location": "py"}).drop(columns=["z_location"])
    # Xenium's own word, CosMx's 0 and blanks are no cell; ids are cells
    words = np.array(["UNASSIGNED", "aaabbbcc-1", "0", "", "cell_7", "none"])
    df["owner"] = words[np.arange(N) % len(words)]
    df.to_parquet(tmp_path / "t.parquet")
    store = store_at(tmp_path / "s.zarr")
    with pytest.raises(SystemExit, match="not in the table"):
        t2p.convert(tmp_path / "t.parquet", store, columns="gene=g,x=px,y=py,cell=nope", log=lambda *_: None)
    t2p.convert(tmp_path / "t.parquet", store, columns="gene=g,x=px,y=py,cell=owner", log=lambda *_: None)
    check_element(store, df, "g", "px", "py", 2)
    got, _, _ = unassigned_of(store, "Gene1")
    no_cell = df.owner.isin(["UNASSIGNED", "0", "", "none"])
    assert int(got.sum()) == int(((df.g == "Gene1") & no_cell).sum())


def test_xenium_and_unnamed_tables_are_not_split(tmp_path):
    """The browser preview doesn't split Xenium's molecules, and a cell column
    is never guessed: no `unassigned` column unless one is named."""
    df = table(6)
    df["cell_id"] = "UNASSIGNED"
    df.to_parquet(tmp_path / "transcripts.parquet")
    store = store_at(tmp_path / "s.zarr")
    t2p.convert(tmp_path / "transcripts.parquet", store, log=lambda *_: None)
    g, _ = pz.open_element(store / "points" / "molecules")
    assert "unassigned" not in list(g["0"]["tile_major"]["columns"][:])
    assert g["0"]["feature_major"].attrs["points:copied_columns"] == []


def test_the_validator_wants_unassigned_in_feature_major(tmp_path, monkeypatch):
    df = table(7).rename(columns={"feature_name": "gene", "x_location": "global_x", "y_location": "global_y"})
    df["cell_id"] = np.where(np.arange(N) % 2 == 0, -1, 5)
    df.to_parquet(tmp_path / "detected_transcripts.parquet")
    # A writer that forgets to copy it
    real = pz.PointsWriter.__init__

    def forgetful(self, *a, **kw):
        kw["copied"] = []
        real(self, *a, **kw)

    monkeypatch.setattr(pz.PointsWriter, "__init__", forgetful)
    store = store_at(tmp_path / "s.zarr")
    t2p.convert(tmp_path / "detected_transcripts.parquet", store, log=lambda *_: None)
    with pytest.raises(AssertionError, match="unassigned is not copied"):
        pz.validate(store / "points" / "molecules")
