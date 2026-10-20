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
