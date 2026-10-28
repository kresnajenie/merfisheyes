"""scripts/slices_to_points.py, judged from outside: many slices' molecules,
one points element, in the same frame as the combined cells.

This is the step the old pipeline did not have. It used to write one
single-molecule folder per slice, each in its own coordinates, while the
cells of every slice were combined into one frame — so the two halves could
not be drawn together. Here the molecules are moved by the same per-slice
shift the cells were moved by, and land in one store beside them.

  python -m pytest tests/bil/test_slices_to_points.py -q
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from conftest import BLANKS, ROOT, SCRIPTS

sys.path.insert(0, str(ROOT / "scripts"))

import points_zarr as pz  # noqa: E402
import validate_store  # noqa: E402


def run(script: str, *args) -> subprocess.CompletedProcess:
    r = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / script), *map(str, args)],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert r.returncode == 0, f"{script} exited {r.returncode}\n{r.stdout[-3000:]}{r.stderr[-3000:]}"
    return r


@pytest.fixture(scope="module")
def store(combined, tmp_path_factory) -> Path:
    """The whole chain: combined cells into a store, then every slice's
    molecules into the points element of that same store."""
    out = tmp_path_factory.mktemp("store")
    run("prepare.py", combined, out, "--name", "test_set")
    store = out / "test_set.zarr"
    run("slices_to_points.py", combined, store)
    return store


@pytest.fixture(scope="module")
def element(store) -> Path:
    return store / "points" / "molecules"


def offsets(combined: Path) -> dict:
    return {s["sample_id"]: s for s in json.loads((combined / "sample_offsets.json").read_text())["samples"]}


def test_the_store_is_valid_with_cells_and_molecules_in_it(store):
    r = validate_store.validate(store, sample=5000)
    assert not r.failed, [m for lv, _, m in r.lines if lv == "FAIL"]
    assert (store / "tables" / "table").exists()
    assert (store / "points" / "molecules").exists()


def test_stats_json_counts_the_molecules(store, dataset):
    """prepare.py writes stats.json before this stage adds the molecules, so
    the points stage has to bring it back in line. It is what a catalog
    reads: a store of hundreds of millions that reports 0 is worse than no
    sidecar at all."""
    stats = json.loads((store.parent / "stats.json").read_text())
    _, attrs = pz.open_element(store / "points" / "molecules")
    assert stats["molecules"] == attrs["points:n_points"] > 0
    assert stats["cells"] == dataset["cells"]


def test_several_label_files_merge(combined, tmp_path, dataset):
    """The cluster passes two: MapMyCells' labels and the artifact mask the
    combine computed. Both are one row per cell keyed by cell_id."""
    import anndata as ad

    mask = tmp_path / "artifact_mask_p25.csv"
    meta = pd.read_csv(combined / "cell_metadata.csv", dtype=str)
    idc = "id" if "id" in meta.columns else meta.columns[0]
    pd.DataFrame({"cell_id": meta[idc], "is_artifact": ["True"] * 10 + ["False"] * (len(meta) - 10)}).to_csv(mask, index=False)
    extra = tmp_path / "labels.csv"
    pd.DataFrame({"cell_id": meta[idc], "celltype": ["Astro"] * len(meta)}).to_csv(extra, index=False)

    out = tmp_path / "two_labels"
    run("prepare.py", combined, out, "--name", "two", "--labels", extra, "--labels", mask)
    obs = ad.read_zarr(out / "two.zarr" / "tables" / "table").obs
    assert "celltype" in obs.columns and "is_artifact" in obs.columns
    assert (obs["is_artifact"].astype(str) == "True").sum() == 10


def test_every_real_molecule_of_every_slice_is_there(element, dataset):
    """One element, all the slices, blanks and control probes left out."""
    _, attrs = pz.open_element(element)
    tx = pd.concat([pd.read_csv(s["dir"] / "cellpose-detected_transcripts.csv") for s in dataset["slices"]])
    assert attrs["points:n_points"] == len(tx[~tx["gene"].isin(BLANKS)])


def test_each_molecule_carries_the_slice_it_came_from(element, dataset, combined):
    """The same tag the cells carry in `_sample_id`. Without it a molecule
    cannot be told back to its slice once the slices share a frame."""
    df = pz.read_dataframe(element, level=0, order="tile_major", columns=["sample"])
    counts = df["sample"].value_counts()
    for s in dataset["slices"]:
        tx = pd.read_csv(s["dir"] / "cellpose-detected_transcripts.csv")
        assert counts[s["id"]] == len(tx[~tx["gene"].isin(BLANKS)])


def test_the_molecules_are_moved_onto_the_cells_frame(element, combined, dataset):
    """Each slice's molecules must sit where that slice's cells were put —
    the point of the whole stage. A shift left out shows up as molecules
    piled on the first slice."""
    df = pz.read_dataframe(element, level=0, order="tile_major", columns=["x", "y", "sample"])
    off = offsets(combined)
    for s in dataset["slices"]:
        o = off[s["id"]]
        x0, y0, x1, y1 = s["box"]
        mine = df[df["sample"] == s["id"]]
        assert len(mine) > 0
        # level 0 holds the coordinates as written, in float32
        assert mine["x"].min() >= x0 + o["dx"] - 0.1
        assert mine["x"].max() <= x1 + o["dx"] + 0.1
        assert mine["y"].min() >= y0 + o["dy"] - 0.1
        assert mine["y"].max() <= y1 + o["dy"] + 0.1


def test_molecules_and_cells_of_a_slice_land_together(store, element, dataset):
    """The cells and the molecules of one slice must overlap, which is the
    only thing that makes the two halves one dataset."""
    import anndata as ad

    obs = ad.read_zarr(store / "tables" / "table").obs
    df = pz.read_dataframe(element, level=0, order="tile_major", columns=["x", "y", "sample"])
    for s in dataset["slices"]:
        cells = obs[obs["_sample_id"].astype(str) == s["id"]]
        mine = df[df["sample"] == s["id"]]
        assert len(cells) == s["cells"]
        # the molecules' box contains the cells' centres, give or take a cell
        assert mine["x"].min() - 5 <= cells["center_x"].min()
        assert mine["x"].max() + 5 >= cells["center_x"].max()
        assert mine["y"].min() - 5 <= cells["center_y"].min()
        assert mine["y"].max() + 5 >= cells["center_y"].max()


def test_molecules_in_no_cell_are_marked(element, dataset):
    """MERSCOPE writes `-1` for a molecule in no cell; the shared rule
    (meyes_common.unassigned_mask) is what decides, not this stage."""
    df = pz.read_dataframe(element, level=0, order="tile_major", columns=["unassigned"])
    tx = pd.concat([pd.read_csv(s["dir"] / "cellpose-detected_transcripts.csv") for s in dataset["slices"]])
    real = tx[~tx["gene"].isin(BLANKS)]
    assert int(np.asarray(df["unassigned"]).sum()) == int((real["cell_id"].astype(str) == "-1").sum())


def test_a_slice_with_no_molecules_is_skipped_not_fatal(combined, tmp_path):
    """A cells-only slice must not stop the others being written."""
    data = json.loads((combined / "sample_offsets.json").read_text())
    data["samples"][0]["transcripts"] = None
    hacked = tmp_path / "combined_output"
    hacked.mkdir()
    (hacked / "sample_offsets.json").write_text(json.dumps(data))

    store = tmp_path / "partial.zarr"
    import zarr

    zarr.open_group(store, mode="w", zarr_format=3).attrs["spatialdata_attrs"] = {"version": "0.2"}
    run("slices_to_points.py", hacked, store)
    _, attrs = pz.open_element(store / "points" / "molecules")
    assert attrs["points:n_points"] > 0
