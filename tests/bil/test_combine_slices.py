"""scripts/combine_slices_v3.py, judged from outside.

The cell combine is the step the whole BIL pipeline hangs off: it is what
puts many slices in one frame, and — since the molecules are far too large to
concatenate the way the cells are — it is also the only place that records
where each slice's molecules went. These tests pin both halves.

  python -m pytest tests/bil/test_combine_slices.py -q
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

# sample ids are digits on the cluster, so every reader must ask for text
IDS = {"id": str, "_sample_id": str}

from conftest import ROOT, SCRIPTS


def run_combine(input_dir: Path, out_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / "combine_slices_v3.py"), str(input_dir), str(out_dir)],
        capture_output=True, text=True, cwd=ROOT,
    )


def offsets(combined: Path) -> dict:
    return json.loads((combined / "sample_offsets.json").read_text())


def test_every_slice_is_found_and_combined(combined, dataset):
    """Non-canonical names and a missing id column do not hide a slice."""
    meta = pd.read_csv(combined / "cell_metadata.csv", dtype=IDS)
    assert len(meta) == dataset["cells"]
    assert set(meta["_sample_id"]) == {s["id"] for s in dataset["slices"]}


def test_the_offsets_name_each_slices_molecules(combined, dataset):
    """The molecules are not combined, so the combine must say where they are.

    Without this the single-molecule half has to re-discover the slices, and
    the two halves can disagree about which slice a file belongs to.
    """
    data = offsets(combined)
    by_id = {s["sample_id"]: s for s in data["samples"]}
    assert set(by_id) == {s["id"] for s in dataset["slices"]}
    for s in dataset["slices"]:
        tx = by_id[s["id"]]["transcripts"]
        assert tx is not None, f"{s['id']}: molecules not recorded"
        assert Path(tx).exists()
        assert Path(tx).name == "cellpose-detected_transcripts.csv"


def test_the_offsets_invert_the_combine(combined, dataset):
    """`original = combined - (dx, dy)` must hold for every cell, or the
    molecules cannot be placed in the cells' frame."""
    meta = pd.read_csv(combined / "cell_metadata.csv", dtype=IDS)
    by_id = {s["sample_id"]: s for s in offsets(combined)["samples"]}
    for sample_id, rows in meta.groupby("_sample_id"):
        o = by_id[str(sample_id)]
        spec = next(s for s in dataset["slices"] if s["id"] == str(sample_id))
        x0, y0, x1, y1 = spec["box"]
        back_x = rows["center_x"] - o["dx"]
        back_y = rows["center_y"] - o["dy"]
        assert back_x.between(x0 - 1e-6, x1 + 1e-6).all()
        assert back_y.between(y0 - 1e-6, y1 + 1e-6).all()


def test_the_slices_do_not_overlap_once_combined(combined):
    """The point of the grid: two slices must not land on top of each other."""
    meta = pd.read_csv(combined / "cell_metadata.csv", dtype=IDS)
    boxes = [
        (g["center_x"].min(), g["center_x"].max(), g["center_y"].min(), g["center_y"].max())
        for _, g in meta.groupby("_sample_id")
    ]
    for i, a in enumerate(boxes):
        for b in boxes[i + 1:]:
            apart = a[1] < b[0] or b[1] < a[0] or a[3] < b[2] or b[3] < a[2]
            assert apart, f"{a} and {b} overlap"


def test_a_rerun_does_not_read_its_own_output(dataset, tmp_path):
    """Run the combine with its output inside the folder it reads — what
    happens whenever a finished sample is processed again, or when a
    downloaded copy keeps the output beside the raw slices. The second run
    must see the same slices, not the slices plus itself."""
    out = tmp_path / "in_place"
    assert run_combine(dataset["input"], out).returncode == 0
    first = pd.read_csv(out / "combined_output" / "cell_metadata.csv")

    # put the output where a rerun would walk over it
    nested = tmp_path / "nested"
    nested.mkdir()
    for s in dataset["slices"]:
        shutil.copytree(s["dir"], nested / s["id"])
    shutil.copytree(out / "combined_output", nested / "combined_output")

    again = tmp_path / "again"
    assert run_combine(nested, again).returncode == 0
    second = pd.read_csv(again / "combined_output" / "cell_metadata.csv")
    assert len(second) == len(first) == dataset["cells"]
