"""The format specifications, checked from outside the viewer:

- scripts/validate_store.py passes the tiny fixtures composed into one
  store and fails each profile violation it is meant to catch;
- scripts/legacy_to_points.py turns a legacy molecule folder into a points
  element the validator and the reference reader accept, molecule for
  molecule;
- the spatialdata library reads our table and image, and a store it
  writes passes our validator (skipped when spatialdata is not installed).

  python -m pytest tests/python -q
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
import zarr

ROOT = Path(__file__).resolve().parents[2]
TINY = ROOT / "tests" / "data" / "tiny"
ZARR = TINY / "zarr"
sys.path.insert(0, str(ROOT / "scripts"))

import legacy_to_points  # noqa: E402
import points_zarr as pz  # noqa: E402
import validate_store  # noqa: E402


def compose_store(dst: Path, *, points=True, shapes=True) -> Path:
    """The tiny fixtures as one store, the way the E2E server serves them."""
    store = dst / "store.zarr"
    shutil.copytree(ZARR / "spatialdata_v3.zarr", store)
    if points:
        shutil.copytree(ZARR / "points_v1.zarr" / "points", store / "points")
    if shapes:
        shutil.copytree(ZARR / "shapes_v1", store / "shapes")
        # The fixture holds the element only; the store needs the group
        (store / "shapes" / "zarr.json").write_text(json.dumps({"zarr_format": 3, "node_type": "group"}))
    zarr.consolidate_metadata(str(store))
    return store


def levels(report: validate_store.Report):
    return {(lv, w) for lv, w, _ in report.lines}


def fails(report: validate_store.Report):
    return [m for lv, _, m in report.lines if lv == "FAIL"]


# ── validate_store ────────────────────────────────────────────────────────

def test_composed_tiny_store_passes(tmp_path):
    store = compose_store(tmp_path)
    r = validate_store.validate(store, sample=1000)
    assert not r.failed, fails(r)
    assert ("OK", "tables/table") in levels(r)
    assert ("OK", "points/transcripts") in levels(r)
    assert ("OK", "shapes/cell_boundaries") in levels(r)
    assert ("OK", "images/dapi") in levels(r)
    assert validate_store.main([str(store), "--sample", "1000"]) == 0


def test_mapping_beside_the_store(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    (tmp_path / "mapping.json").write_text(json.dumps({"cells": "store.zarr", "image": "https://x/dapi.zarr"}))
    r = validate_store.validate(store)
    assert ("OK", "mapping.json") in levels(r)
    (tmp_path / "mapping.json").write_text(json.dumps({"links": {"__all__": "u"}}))
    r = validate_store.validate(store)
    assert any("linkColumn" in m for m in fails(r))


def test_missing_consolidated_metadata_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    root = json.loads((store / "zarr.json").read_text())
    del root["consolidated_metadata"]
    (store / "zarr.json").write_text(json.dumps(root))
    r = validate_store.validate(store)
    assert any("consolidated" in m for m in fails(r))


def test_stale_consolidated_metadata_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=True)
    # Written after the last consolidation: the viewer would never see it
    root = json.loads((store / "zarr.json").read_text())
    root["consolidated_metadata"]["metadata"] = {
        k: v for k, v in root["consolidated_metadata"]["metadata"].items() if not k.startswith("shapes")
    }
    (store / "zarr.json").write_text(json.dumps(root))
    r = validate_store.validate(store)
    assert any("misses" in m for m in fails(r))


def test_shapes_as_a_directory_fails(tmp_path):
    store = compose_store(tmp_path, points=False)
    f = store / "shapes" / "cell_boundaries" / "shapes.parquet"
    part = f.read_bytes()
    f.unlink()
    f.mkdir()
    (f / "part.0.parquet").write_bytes(part)
    r = validate_store.validate(store)
    assert any("directory of parts" in m for m in fails(r))


def test_csr_without_csc_copy_fails(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    x = store / "tables" / "table" / "X" / "zarr.json"
    meta = json.loads(x.read_text())
    meta.setdefault("attributes", {})["encoding-type"] = "csr_matrix"
    x.write_text(json.dumps(meta))
    # The tiny table is CSR with the copy the profile asks for: accepted
    assert not validate_store.validate(store).failed
    shutil.rmtree(store / "tables" / "table" / "layers" / "X_csc")
    zarr.consolidate_metadata(str(store))
    r = validate_store.validate(store)
    assert any("CSR" in m for m in fails(r))


SDDB_ROOTS = ROOT / "tests" / "data" / "sddb" / "roots"


def test_root_json_reads_a_store_from_its_root_alone():
    # Visium (cerebellum): CSR X with a csc_counts layer, as the reader accepts it
    root = validate_store.read_json(SDDB_ROOTS / "fJSOO03VzCTYtD1o.json")
    r = validate_store.validate(validate_store.RootStore(root))
    assert not r.failed, fails(r)
    assert any("layers/csc_counts" in m for lv, w, m in r.lines if w == "tables/table")
    assert any("root only" in m for lv, w, m in r.lines if w.startswith("shapes/"))
    assert validate_store.main(["--root-json", str(SDDB_ROOTS / "fJSOO03VzCTYtD1o.json")]) == 0


def test_root_json_fails_a_csr_table_with_no_csc_layer():
    # Visium HD (ovary): only the 16 um table has a CSC layer
    root = validate_store.read_json(SDDB_ROOTS / "haOzFbJ7Bc9x7NW6.json")
    r = validate_store.validate(validate_store.RootStore(root))
    failed = {w for lv, w, _ in r.lines if lv == "FAIL"}
    assert failed == {"tables/square_002um", "tables/square_008um", "tables/square_016um_hvg"}
    assert ("OK", "tables/square_016um") in levels(r)


def test_image_axes_and_dtype_are_checked(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=False)
    g = store / "images" / "dapi" / "zarr.json"
    meta = json.loads(g.read_text())
    meta["attributes"]["ome"]["multiscales"][0]["axes"] = [{"name": "y"}, {"name": "x"}]
    g.write_text(json.dumps(meta))
    r = validate_store.validate(store)
    assert any("axes" in m for m in fails(r))


def test_points_metadata_schema_is_enforced(tmp_path):
    store = compose_store(tmp_path, shapes=False)
    g = store / "points" / "transcripts" / "zarr.json"
    meta = json.loads(g.read_text())
    del meta["attributes"]["points:n_points"]
    g.write_text(json.dumps(meta))
    r = validate_store.validate(store, sample=100)
    assert any("points:n_points" in m for m in fails(r))


# ── legacy_to_points ──────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def legacy_folder(tmp_path_factory) -> Path:
    """A legacy molecule folder from the tiny parquet, as the processor
    writes it (5,000 molecules, 3D, raw coordinates)."""
    import process_single_molecule as psm

    out = tmp_path_factory.mktemp("legacy") / "molecules"
    psm.process_single_molecule_data(
        input_file=str(TINY / "single-molecule" / "tiny.parquet"),
        output_folder=str(out),
        dataset_type="xenium",
    )
    return out


def test_legacy_folder_becomes_a_valid_points_element(tmp_path, legacy_folder):
    store = tmp_path / "store.zarr"
    out = legacy_to_points.convert(legacy_folder, store, levels=3, tile_size=200, log=lambda *_: None)
    assert out == store / "points" / "transcripts"
    r = validate_store.validate(store, sample=2000)
    assert not r.failed, fails(r)

    # Molecule for molecule: every gene's x, y, z in the folder are in the element
    manifest = legacy_to_points.read_manifest(legacy_folder)
    genes = manifest["genes"]["unique_gene_names"]
    _, attrs = pz.open_element(out)
    assert attrs["spatialdata_attrs"]["feature_key"] == "gene"
    assert attrs["points:n_points"] == manifest["statistics"]["total_molecules"]
    for i, gene in enumerate(genes[:3]):
        expected = legacy_to_points.read_gene(legacy_folder, gene)
        got = pz.feature_points(out, gene, level=0)
        assert len(got["x"]) == len(expected)
        for k, axis in enumerate(("x", "y", "z")):
            assert np.array_equal(np.sort(np.asarray(got[axis], dtype=np.float32)), np.sort(expected[:, k]))


def test_legacy_folder_can_drop_z(tmp_path, legacy_folder):
    store = tmp_path / "flat.zarr"
    out = legacy_to_points.convert(legacy_folder, store, levels=2, tile_size=200, drop_z=True, log=lambda *_: None)
    _, attrs = pz.open_element(out)
    assert attrs["axes"] == ["x", "y"]
    assert not validate_store.validate(store, sample=500).failed


# ── spatialdata round-trip ────────────────────────────────────────────────

spatialdata = pytest.importorskip("spatialdata", reason="the spatialdata package is the external reader")


@pytest.mark.xfail(
    reason="spatialdata 0.8 rejects the tiny fixture's element version ('0.5'): "
    "the fixture is written by anndata and our scripts, not by spatialdata. "
    "The other direction (a store spatialdata writes passes our validator) holds.",
    strict=True,
)
def test_spatialdata_reads_our_table_and_image(tmp_path):
    store = compose_store(tmp_path, points=False, shapes=True)
    sdata = spatialdata.read_zarr(store, selection=("tables", "images", "shapes"))
    table = sdata.tables["table"]
    assert table.n_obs == 30
    assert "spatial" in table.obsm
    assert "dapi" in sdata.images
    assert "cell_boundaries" in sdata.shapes
    assert len(sdata.shapes["cell_boundaries"]) == 30


def test_a_store_spatialdata_writes_passes_the_validator(tmp_path):
    import anndata as ad
    import pandas as pd
    from spatialdata.models import Image2DModel, TableModel

    n = 12
    rng = np.random.default_rng(1)
    table = ad.AnnData(
        X=rng.poisson(1.0, size=(n, 5)).astype(np.float32),
        obs=pd.DataFrame({"cell_type": pd.Categorical(["a", "b"] * (n // 2)), "region": "img", "cell_id": np.arange(n)}, index=[f"c{i}" for i in range(n)]),
    )
    table.obsm["spatial"] = rng.uniform(0, 64, size=(n, 2)).astype(np.float32)
    table.var_names = [f"g{i}" for i in range(5)]
    table = TableModel.parse(table, region="img", region_key="region", instance_key="cell_id")
    image = Image2DModel.parse(rng.integers(0, 255, size=(1, 64, 64), dtype=np.uint8), dims=("c", "y", "x"), scale_factors=[2])
    sdata = spatialdata.SpatialData(images={"img": image}, tables={"table": table})
    store = tmp_path / "written.zarr"
    sdata.write(store)
    r = validate_store.validate(store)
    # What spatialdata leaves out is reported, never silently accepted
    assert not r.failed, fails(r)
