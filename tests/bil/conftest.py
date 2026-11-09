"""A tiny multi-slice MERSCOPE dataset, shaped like the BIL ones.

The real slices on the cluster are named by the segmentation that produced
them (`cellpose_metadata.csv`, `cellpose-cell-by-gene.csv`,
`cellpose-detected_transcripts.csv`), sit several directories below the
sample folder, and do not agree with each other about the cell id column —
one slice names it `id`, another leaves it unnamed. The fixture reproduces
all of that at a size that runs in a second, so the tests exercise the
matching the pipeline actually needs rather than a tidy case.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

SCRIPTS = ROOT / "scripts"
GENES = ["Slc17a7", "Gad2", "Drd1", "Aqp4", "Apoe", "Pdgfra"]
BLANKS = ["Blank-1", "NegControlProbe_3"]

# Each slice: how many cells, how many molecules, the box its cells sit in,
# and whether its metadata names the id column. Different boxes so the grid
# has to move them apart, and so a shift that was not applied shows up.
SLICES = [
    {"id": "1232698492", "cells": 120, "molecules": 3000, "box": (0.0, 0.0, 100.0, 80.0), "named_id": True},
    {"id": "1232698496", "cells": 90, "molecules": 2200, "box": (20.0, 10.0, 140.0, 70.0), "named_id": False},
]


def _slice(folder: Path, spec: dict, rng: np.random.Generator) -> dict:
    """One slice's three CSVs, under the directory depth BIL uses."""
    d = folder / spec["id"] / "merfish_output" / "cellpose_cyto2_nuclei"
    d.mkdir(parents=True)
    x0, y0, x1, y1 = spec["box"]
    n = spec["cells"]

    # cell ids are the long integers the instrument writes, unique per slice
    # 21-digit ids, as the instrument writes them (too wide for an int64, so
    # built as text — which is also how every reader has to handle them)
    ids = [f"{rng.integers(10**17, 10**18)}{i:04d}" for i in range(n)]
    cx = rng.uniform(x0, x1, n).round(3)
    cy = rng.uniform(y0, y1, n).round(3)
    meta = pd.DataFrame({
        "id": ids, "fov": rng.integers(0, 4, n), "volume": rng.uniform(100, 900, n).round(3),
        "center_x": cx, "center_y": cy,
        "min_x": cx - 3, "min_y": cy - 3, "max_x": cx + 3, "max_y": cy + 3,
        "barcodeCount": rng.integers(10, 400, n),
    })
    if spec["named_id"]:
        meta.to_csv(d / "cellpose_metadata.csv", index=False)
    else:
        # the other shape: no name on the id column at all
        meta.drop(columns=["id"]).set_index(pd.Index(ids)).to_csv(d / "cellpose_metadata.csv")

    counts = pd.DataFrame(rng.integers(0, 12, (n, len(GENES) + len(BLANKS))), columns=GENES + BLANKS)
    counts.insert(0, "cell", ids)
    counts.to_csv(d / "cellpose-cell-by-gene.csv", index=False)

    m = spec["molecules"]
    # a quarter of the molecules are in no cell, written as MERSCOPE writes it
    cell_of = np.where(rng.random(m) < 0.25, "-1", rng.choice(ids, m))
    pd.DataFrame({
        "barcode_id": rng.integers(0, 300, m),
        "global_x": rng.uniform(x0, x1, m).round(4),
        "global_y": rng.uniform(y0, y1, m).round(4),
        "global_z": rng.uniform(0, 10, m).round(4),
        "fov": rng.integers(0, 4, m),
        "gene": rng.choice(GENES + BLANKS, m),
        "transcript_id": [f"ENSMUST{i:011d}" for i in range(m)],
        "cell_id": cell_of,
    }).to_csv(d / "cellpose-detected_transcripts.csv")  # unnamed index, as BIL has

    return {**spec, "dir": d, "ids": ids}


@pytest.fixture(scope="session")
def dataset(tmp_path_factory) -> dict:
    """A sample folder of slices, untouched by the tests that read it."""
    folder = tmp_path_factory.mktemp("bil_input") / "ace-test-set"
    folder.mkdir()
    rng = np.random.default_rng(0)
    slices = [_slice(folder, s, rng) for s in SLICES]
    return {
        "input": folder,
        "slices": slices,
        "cells": sum(s["cells"] for s in slices),
        "molecules": sum(s["molecules"] for s in slices),
    }


@pytest.fixture(scope="session")
def combined(dataset, tmp_path_factory) -> Path:
    """`dataset` put through scripts/combine_slices_v3.py: one combined cell
    table, the artifact masks, and sample_offsets.json."""
    import subprocess

    out = tmp_path_factory.mktemp("combined")
    r = subprocess.run(
        [sys.executable, "-I", str(SCRIPTS / "combine_slices_v3.py"), str(dataset["input"]), str(out)],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    return out / "combined_output"
