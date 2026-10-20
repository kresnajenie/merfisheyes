"""scripts/prepare_xenium.py on a real Xenium export, judged from outside:
the validator passes, the store holds what the profiles need, the
spatialdata library reads it back, and the molecules sit on the cells.

Runs only when MERFISHEYES_XENIUM names a Xenium output folder (an export
is gigabytes; none is committed):

  MERFISHEYES_XENIUM=/path/to/Xenium_outs python -m pytest tests/python/test_prepare_xenium.py -q
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

EXPORT = os.environ.get("MERFISHEYES_XENIUM")
pytestmark = pytest.mark.skipif(not EXPORT, reason="MERFISHEYES_XENIUM not set")


@pytest.fixture(scope="module")
def prepared(tmp_path_factory):
    pytest.importorskip("spatialdata_io")
    import prepare_xenium

    folder = tmp_path_factory.mktemp("xenium") / "sample"
    report = prepare_xenium.prepare_xenium(Path(EXPORT), folder, name="sample", workers=4, log=lambda *_: None)
    return folder, report


def test_the_folder_passes_the_validator(prepared):
    folder, report = prepared
    assert not report.failed, [m for lv, _, m in report.lines if lv == "FAIL"]
    assert any(w.startswith("molecules/") and "covering" in m for lv, w, m in report.lines if lv == "OK")


def test_the_store_holds_what_the_profiles_need(prepared):
    folder, _ = prepared
    store = folder / "sample.zarr"
    nodes = json.loads((store / "zarr.json").read_text())["consolidated_metadata"]["metadata"]
    table = nodes["tables/table"]["attributes"]
    assert table["region"] == "cell_boundaries" and table["instance_key"] == "cell_id"
    assert nodes["tables/table/layers/X_csc"]["attributes"]["encoding-type"] == "csc_matrix"
    assert (store / "shapes" / "cell_boundaries" / "shapes.parquet").is_file()
    images = [k for k in nodes if k.startswith("images/") and k.count("/") == 1]
    assert images, "no image written"
    for level in (k for k in nodes if k.startswith("images/") and k.count("/") == 2):
        chunks = nodes[level]["chunk_grid"]["configuration"]["chunk_shape"]
        assert max(chunks[-2:]) <= 1024, (level, chunks)
    assert not any(k.startswith("labels/") for k in nodes)
    mapping = json.loads((folder / "mapping.json").read_text())
    assert mapping == {"cells": "sample.zarr", "linkColumn": "__all__", "links": {"__all__": "molecules"}}


def test_spatialdata_reads_the_store_back(prepared):
    spatialdata = pytest.importorskip("spatialdata")
    folder, _ = prepared
    sdata = spatialdata.read_zarr(folder / "sample.zarr", selection=("tables", "shapes"))
    assert "table" in sdata.tables and "cell_boundaries" in sdata.shapes
    assert sdata.tables["table"].n_obs == len(sdata.shapes["cell_boundaries"])
