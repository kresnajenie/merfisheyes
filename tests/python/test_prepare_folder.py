"""scripts/prepare_folder.py, judged from outside: the folder it leaves
behind passes scripts/validate_store.py, the table it writes is a Zarr v3
SpatialData store with a CSC copy of X, and the parts it names in
mapping.json are the ones in the folder.

  python -m pytest tests/python/test_prepare_folder.py -q
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import zarr

ROOT = Path(__file__).resolve().parents[2]
TINY = ROOT / "tests" / "data" / "tiny"
sys.path.insert(0, str(ROOT / "scripts"))

import prepare_folder  # noqa: E402

CELLS, GENES = 40, 6


def make_table(*, categorical=True) -> ad.AnnData:
    rng = np.random.default_rng(1)
    dense = rng.poisson(0.7, size=(CELLS, GENES)).astype(np.float32)
    obs = pd.DataFrame(
        {"cell_id": [f"cell-{i}" for i in range(CELLS)]},
        index=[str(i) for i in range(CELLS)],
    )
    if categorical:
        obs["celltype"] = pd.Categorical([["a", "b"][i % 2] for i in range(CELLS)])
    return ad.AnnData(
        X=sp.csr_matrix(dense),
        obs=obs,
        var=pd.DataFrame(index=[f"gene{i}" for i in range(GENES)]),
        obsm={"spatial": rng.uniform(0, 100, size=(CELLS, 2))},
    )


def write_table(adata: ad.AnnData, path: Path, zarr_format: int) -> Path:
    with ad.settings.override(zarr_write_format=zarr_format):
        adata.write_zarr(path)
    return path


def molecules_into(folder: Path) -> Path:
    import process_single_molecule as psm

    out = folder / "molecules"
    psm.process_single_molecule_data(
        input_file=str(TINY / "single-molecule" / "tiny.parquet"),
        output_folder=str(out),
        dataset_type="xenium",
    )
    return out


def image_into(folder: Path) -> Path:
    """The tiny fixture's image group, taken out of its store: an OME-Zarr
    image on its own (its zarr.json carries the multiscales)."""
    out = folder / "dapi.zarr"
    shutil.copytree(TINY / "zarr" / "spatialdata_v3.zarr" / "images" / "dapi", out)
    return out


def fails(report):
    return [m for lv, _, m in report.lines if lv == "FAIL"]


def consolidated(store: Path) -> dict:
    return json.loads((store / "zarr.json").read_text())["consolidated_metadata"]["metadata"]


# ── a lone table, in each form it is written ─────────────────────────────

@pytest.mark.parametrize("form", ["zarr2", "zarr3", "h5ad"])
def test_lone_table_becomes_a_store_with_a_csc_copy(tmp_path, form):
    folder = tmp_path / "sample"
    folder.mkdir()
    adata = make_table()
    if form == "h5ad":
        adata.write_h5ad(folder / "cells.h5ad")
    else:
        write_table(adata, folder / "cells.zarr", 2 if form == "zarr2" else 3)
    molecules_into(folder)
    image_into(folder)

    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)

    store = folder / "sample.zarr"
    root = json.loads((store / "zarr.json").read_text())
    assert root["zarr_format"] == 3
    assert root["attributes"]["spatialdata_attrs"] == {"version": "0.2"}
    nodes = consolidated(store)
    assert nodes["tables/table/X"]["attributes"]["encoding-type"] == "csr_matrix"
    assert nodes["tables/table/layers/X_csc"]["attributes"]["encoding-type"] == "csc_matrix"
    assert nodes["tables/table/layers/X_csc"]["attributes"]["shape"] == [CELLS, GENES]
    # The copy holds the same counts
    table = zarr.open_group(store, mode="r")["tables/table"]
    assert (ad.io.read_elem(table["layers/X_csc"]).toarray() == ad.io.read_elem(table["X"]).toarray()).all()

    mapping = json.loads((folder / "mapping.json").read_text())
    assert mapping == {
        "cells": "sample.zarr",
        "linkColumn": "__all__",
        "links": {"__all__": "molecules"},
        "image": "dapi.zarr",
    }
    # A Zarr v3 table is moved, not copied; the others are left in place
    original = folder / ("cells.h5ad" if form == "h5ad" else "cells.zarr")
    assert original.exists() == (form != "zarr3")


def test_the_guides_show_the_script_as_it_is():
    """The tutorial is generated from the script, not written beside it:
    the command the pages show is the script's own usage line."""
    usage = prepare_folder.__doc__.split("\n")[2].strip()
    assert usage.startswith("python scripts/prepare_folder.py <folder>")
    for page in (ROOT / "README.md", ROOT / "app" / "docs" / "prepare" / "page.tsx"):
        assert "python scripts/prepare_folder.py " in page.read_text(), page


def test_url_makes_the_links_absolute(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    write_table(make_table(), folder / "t.zarr", 3)
    molecules_into(folder)
    r = prepare_folder.prepare(folder, url="https://bucket/data/s/", log=lambda *_: None)
    assert not r.failed, fails(r)
    mapping = json.loads((folder / "mapping.json").read_text())
    assert mapping["cells"] == "https://bucket/data/s/s.zarr"
    assert mapping["links"]["__all__"] == "https://bucket/data/s/molecules"


# ── a store ───────────────────────────────────────────────────────────────

def test_a_store_with_csr_gets_the_copy_in_place_and_a_rerun_changes_nothing(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    write_table(make_table(), folder / "cells.zarr", 3)
    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)
    store = folder / "s.zarr"
    before = (store / "tables" / "table" / "layers" / "X_csc" / "zarr.json").stat().st_mtime_ns

    # Second run: the folder now holds the store only; nothing is rewritten
    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)
    assert (store / "tables" / "table" / "layers" / "X_csc" / "zarr.json").stat().st_mtime_ns == before
    assert sorted(p.name for p in folder.iterdir()) == ["mapping.json", "s.zarr"]
    assert json.loads((folder / "mapping.json").read_text()) == {"cells": "s.zarr"}


def test_a_csc_table_is_left_alone(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    adata = make_table()
    adata.X = sp.csc_matrix(adata.X)
    write_table(adata, folder / "cells.zarr", 3)
    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)
    nodes = consolidated(folder / "s.zarr")
    assert nodes["tables/table/X"]["attributes"]["encoding-type"] == "csc_matrix"
    assert "tables/table/layers/X_csc" not in nodes


def test_x_spatial_is_accepted(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    adata = make_table()
    adata.obsm["X_spatial"] = adata.obsm.pop("spatial")
    write_table(adata, folder / "cells.zarr", 3)
    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)
    assert "tables/table/obsm/spatial" in consolidated(folder / "s.zarr")


# ── what it refuses ───────────────────────────────────────────────────────

def test_two_cell_candidates_need_a_choice(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    write_table(make_table(), folder / "a.zarr", 3)
    write_table(make_table(), folder / "b.zarr", 3)
    with pytest.raises(SystemExit, match="several cells candidates"):
        prepare_folder.prepare(folder, log=lambda *_: None)
    r = prepare_folder.prepare(folder, cells=folder / "a.zarr", log=lambda *_: None)
    assert not r.failed, fails(r)


def test_no_coordinates_is_an_error(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    adata = make_table()
    del adata.obsm["spatial"]
    write_table(adata, folder / "cells.zarr", 3)
    with pytest.raises(SystemExit, match="obsm/spatial"):
        prepare_folder.prepare(folder, log=lambda *_: None)


def test_molecules_and_image_without_cells(tmp_path):
    folder = tmp_path / "s"
    folder.mkdir()
    molecules_into(folder)
    image_into(folder)
    r = prepare_folder.prepare(folder, log=lambda *_: None)
    assert not r.failed, fails(r)
    mapping = json.loads((folder / "mapping.json").read_text())
    assert mapping == {"linkColumn": "__all__", "links": {"__all__": "molecules"}, "image": "dapi.zarr"}


# ── clustering ────────────────────────────────────────────────────────────

def test_cluster_runs_only_when_no_categorical_column(tmp_path):
    pytest.importorskip("scanpy")
    folder = tmp_path / "s"
    folder.mkdir()
    write_table(make_table(categorical=True), folder / "cells.zarr", 3)
    r = prepare_folder.prepare(folder, cluster=True, log=lambda *_: None)
    assert not r.failed, fails(r)
    assert not (folder / "annotations").exists()

    other = tmp_path / "u"
    other.mkdir()
    write_table(make_table(categorical=False), other / "cells.zarr", 3)
    r = prepare_folder.prepare(other, cluster=True, log=lambda *_: None)
    assert not r.failed, fails(r)
    notes = json.loads((other / "annotations" / "annotations.json").read_text())
    assert notes["num_cells"] == CELLS and notes["de_stats"] == ["leiden"]
    codes = np.fromfile(other / "annotations" / "leiden.bin", dtype="<u2")
    assert len(codes) == CELLS and len(notes["columns"]["leiden"]["categories"]) == codes.max() + 1
