"""scripts/prepare_merscope.py on a MERSCOPE-shaped export made here (the
instrument's CSVs and a vpt boundaries parquet), judged from outside: the
validator passes, the store is one the spatialdata library reads back, the
blanks are gone, the polygons match the cells, the molecules sit on them.

  python -m pytest tests/python/test_prepare_merscope.py -q
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

spatialdata = pytest.importorskip("spatialdata")
gpd = pytest.importorskip("geopandas")

import prepare_merscope  # noqa: E402

CELLS, GENES, PLANES = 60, 20, 3


def make_export(folder: Path, *, boundaries=True, categories=False, shuffled=True) -> Path:
    """A region folder as MERSCOPE writes it, for a 400 × 300 µm piece."""
    import shapely

    rng = np.random.default_rng(3)
    folder.mkdir(parents=True)
    ids = [str(1_000_000_000 + i) for i in range(CELLS)]
    centres = rng.uniform([20, 20], [380, 280], size=(CELLS, 2))
    meta = pd.DataFrame(
        {
            "EntityID": ids,
            "fov": rng.integers(0, 4, CELLS),
            "volume": rng.uniform(500, 2000, CELLS),
            "center_x": centres[:, 0],
            "center_y": centres[:, 1],
        }
    )
    meta.to_csv(folder / "cell_metadata.csv", index=False)

    genes = [f"Gene{i}" for i in range(GENES)] + ["Blank-1", "Blank-2"]
    counts = pd.DataFrame(rng.poisson(0.9, size=(CELLS, len(genes))), columns=genes)
    counts.insert(0, "cell", ids)
    if shuffled:  # a different row order than the metadata
        counts = counts.sample(frac=1, random_state=1)
    counts.to_csv(folder / "cell_by_gene.csv", index=False)

    n = 3000
    tx = pd.DataFrame(
        {
            "barcode_id": rng.integers(0, GENES, n),
            "global_x": rng.uniform(0, 400, n),
            "global_y": rng.uniform(0, 300, n),
            "global_z": rng.uniform(0, 10, n),
            "gene": rng.choice(genes[:GENES], n),
            "cell_id": -1,
        }
    )
    tx.to_csv(folder / "detected_transcripts.csv")

    if boundaries:
        rows = []
        for i, (cx, cy) in enumerate(centres):
            for zi in range(PLANES):
                r = 6 + zi
                square = shapely.MultiPolygon([shapely.box(cx - r, cy - r, cx + r, cy + r)])
                rows.append({"ID": i * PLANES + zi, "EntityID": int(ids[i]), "ZIndex": zi, "Type": "cell",
                             "ZLevel": 1.5 * zi, "Geometry": square})
        gdf = gpd.GeoDataFrame(pd.DataFrame(rows), geometry="Geometry")
        gdf.to_parquet(folder / "cell_boundaries.parquet")

    if categories:
        pd.DataFrame({"EntityID": ids, "celltype": [["a", "b", "c"][i % 3] for i in range(CELLS)]}).to_csv(
            folder / "cell_categories.csv", index=False
        )
    return folder


def fails(report):
    return [m for lv, _, m in report.lines if lv == "FAIL"]


def test_an_export_becomes_a_store_with_polygons_and_molecules(tmp_path):
    export = make_export(tmp_path / "region_R1")
    folder = tmp_path / "sample"
    r = prepare_merscope.prepare_merscope(export, folder, z=1, log=lambda *_: None)
    assert not r.failed, fails(r)

    store = folder / "region_R1.zarr"
    nodes = json.loads((store / "zarr.json").read_text())["consolidated_metadata"]["metadata"]
    table = nodes["tables/table"]["attributes"]
    assert table["region"] == "cell_boundaries" and table["instance_key"] == "cell_id"
    assert nodes["tables/table/X"]["attributes"]["encoding-type"] == "csc_matrix"
    assert nodes["tables/table/X"]["attributes"]["shape"] == [CELLS, GENES]  # the blanks are gone
    assert (store / "shapes" / "cell_boundaries" / "shapes.parquet").is_file()
    # The molecules went inside the store as a tiled points element; the folder to sources/
    assert nodes["points/molecules"]["attributes"]["points:version"]
    assert (folder / "sources" / "molecules" / "manifest.json.gz").exists()
    assert json.loads((folder / "mapping.json").read_text()) == {"cells": "region_R1.zarr"}

    # Read back by the library that wrote it: the polygons are the cells', one each, on plane 1
    sdata = spatialdata.read_zarr(store, selection=("tables", "shapes"))
    adata = sdata.tables["table"]
    polys = sdata.shapes["cell_boundaries"]
    assert adata.n_obs == len(polys) == CELLS
    assert set(adata.obs["cell_id"]) == set(polys.index)
    assert all(abs(p.bounds[2] - p.bounds[0] - 14) < 1e-6 for p in polys.geometry)  # r = 6 + 1
    # The counts follow the metadata's order, not the shuffled file's
    meta = pd.read_csv(export / "cell_metadata.csv", dtype={"EntityID": str})
    assert list(adata.obs_names) == list(meta["EntityID"])
    counts = pd.read_csv(export / "cell_by_gene.csv", index_col=0, dtype={0: str}).loc[meta["EntityID"]]
    assert (adata.X.toarray() == counts[[f"Gene{i}" for i in range(GENES)]].to_numpy()).all()


def test_without_boundaries_and_with_categories(tmp_path):
    export = make_export(tmp_path / "region_R2", boundaries=False, categories=True, shuffled=False)
    folder = tmp_path / "sample"
    r = prepare_merscope.prepare_merscope(export, folder, cluster=True, log=lambda *_: None)
    assert not r.failed, fails(r)
    store = folder / "region_R2.zarr"
    nodes = json.loads((store / "zarr.json").read_text())["consolidated_metadata"]["metadata"]
    assert "shapes/cell_boundaries" not in nodes
    assert nodes["tables/table/obs/celltype"]["attributes"]["encoding-type"] == "categorical"
    assert not (folder / "annotations").exists()  # cell types exist: no Leiden


def test_a_missing_plane_is_refused(tmp_path):
    export = make_export(tmp_path / "region_R3")
    with pytest.raises(SystemExit, match="no polygons on ZIndex 7"):
        prepare_merscope.prepare_merscope(export, tmp_path / "s", z=7, log=lambda *_: None)
