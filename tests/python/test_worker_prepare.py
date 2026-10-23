"""The worker's prepare pipeline (worker/entrypoint.py, docs/plans/ONE-PIPELINE.md
step 2): the scripts/prepare.py call it makes for an upload, its reading of
stats.json, and the reason it gives when prepare.py stops — and, run for
real, that the store it writes is the one prepare.py writes on its own."""
from __future__ import annotations

import json
import subprocess
import sys
import types
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "worker"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.modules.setdefault("boto3", types.ModuleType("boto3"))

import entrypoint  # noqa: E402


def molecules(path: Path) -> Path:
    pq.write_table(
        pa.table({"Gene": ["A", "B", "A"], "X": [0.0, 1.0, 2.0], "Y": [0.0, 1.0, 2.0]}),
        path,
    )
    return path


def test_cells_are_prepared_from_the_raw_folder_with_their_labels(tmp_path):
    raw = tmp_path / "raw"
    (raw / "__annotations__").mkdir(parents=True)
    (raw / "set1.h5ad").write_bytes(b"")
    (raw / "__annotations__" / "mapping.csv").write_text("cell_id,label\n")

    cmd = entrypoint.prepare_command(
        "single_cell", raw, tmp_path, tmp_path / "out", "ds_x", {},
        {"mmcCsv": "__annotations__/mapping.csv"},
    )

    assert cmd[1:5] == ["scripts/prepare.py", str(raw), str(tmp_path / "out"), "--name"]
    assert cmd[5] == "ds_x"
    # The labels were moved out of the folder prepare.py reads
    assert cmd[6] == "--labels" and not (raw / "__annotations__").exists()


def test_molecules_are_prepared_from_the_one_table_with_the_chosen_columns(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    table = molecules(raw / "spots.parquet")

    cmd = entrypoint.prepare_command(
        "single_molecule", raw, tmp_path, tmp_path / "out", "sm_x",
        {"columnMapping": {"gene": "Gene", "x": "X", "y": "Y", "z": ""}}, {},
    )

    assert cmd[2] == str(table)
    assert cmd[-2:] == ["--molecule-columns", "Gene,X,Y"]


def test_stats_json_gives_the_rows_counts(tmp_path):
    (tmp_path / "stats.json").write_text(json.dumps({"cells": 0, "genes": 0, "molecules": 3, "molecule_genes": 2, "dims": 2}))
    assert entrypoint.read_prepare_stats(tmp_path) == {"numCells": 3, "numGenes": 2, "spatialDimensions": 2}

    (tmp_path / "stats.json").write_text(json.dumps({"cells": 10, "genes": 5, "molecules": 3, "dims": 3}))
    assert entrypoint.read_prepare_stats(tmp_path)["numCells"] == 10


def test_the_reason_prepare_py_gives_is_the_one_reported(tmp_path):
    returncode, reason = entrypoint.run_processor([
        sys.executable, "-c",
        "print('reading ...'); print('error: raw: not an AnnData table, a Xenium export'); raise SystemExit(2)",
    ])
    assert returncode == 2
    assert reason == "raw: not an AnnData table, a Xenium export"

    returncode, reason = entrypoint.run_processor([
        sys.executable, "-c", "print('FAIL  tables/table: X has no rows'); raise SystemExit(4)",
    ])
    assert reason == "the store failed validation: tables/table: X has no rows"


def test_the_worker_writes_what_prepare_py_writes(tmp_path, monkeypatch):
    """The same command on the same input: the worker's store and a person's
    are the same, file for file."""
    raw = tmp_path / "raw"
    raw.mkdir()
    molecules(raw / "spots.parquet")
    monkeypatch.chdir(ROOT)

    worker_cmd = entrypoint.prepare_command("single_molecule", raw, tmp_path, tmp_path / "worker", "sm_x", {}, {})
    assert entrypoint.run_processor(worker_cmd)[0] == 0
    person = subprocess.run(
        [sys.executable, "scripts/prepare.py", str(raw / "spots.parquet"), str(tmp_path / "person"), "--name", "sm_x"],
        capture_output=True, text=True,
    )
    assert person.returncode == 0, person.stdout + person.stderr

    def files(root: Path) -> dict[str, bytes]:
        return {
            str(p.relative_to(root)): p.read_bytes()
            for p in sorted(root.rglob("*"))
            if p.is_file() and p.name != "stats.json"
        }

    assert files(tmp_path / "worker") == files(tmp_path / "person")
    assert json.loads((tmp_path / "worker" / "mapping.json").read_text()) == {"cells": "sm_x.zarr"}


def test_the_molecules_z_is_kept_only_when_asked(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    molecules(raw / "spots.parquet")

    asked = entrypoint.prepare_command("single_molecule", raw, tmp_path, tmp_path / "out", "sm_x", {"keepZ": True}, {})
    assert asked[-1] == "--keep-z"
    assert "--keep-z" not in entrypoint.prepare_command("single_molecule", raw, tmp_path, tmp_path / "out", "sm_x", {}, {})


def points_axes(out: Path, name: str) -> list[str]:
    attrs = json.loads((out / f"{name}.zarr" / "points" / "molecules" / "zarr.json").read_text())["attributes"]
    return attrs["axes"]


def test_keep_z_writes_3d_molecules_beside_cells_and_alone(tmp_path, monkeypatch):
    """--keep-z: the molecules keep their z, over the cells too (whose 2-D
    transforms then pass z through), and the store still validates."""
    tiny = ROOT / "tests" / "data" / "tiny"
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "cells.h5ad").write_bytes((tiny / "h5ad" / "tiny.h5ad").read_bytes())
    (raw / "transcripts.parquet").write_bytes((tiny / "single-molecule" / "tiny.parquet").read_bytes())
    monkeypatch.chdir(ROOT)

    def run(source: Path, out: Path, *flags: str) -> None:
        r = subprocess.run(
            [sys.executable, "scripts/prepare.py", str(source), str(out), "--name", "s", *flags],
            capture_output=True, text=True,
        )
        assert r.returncode == 0, r.stdout + r.stderr

    run(raw, tmp_path / "flat")
    assert points_axes(tmp_path / "flat", "s") == ["x", "y"]
    run(raw, tmp_path / "deep", "--keep-z")
    assert points_axes(tmp_path / "deep", "s") == ["x", "y", "z"]
    run(raw / "transcripts.parquet", tmp_path / "alone", "--keep-z")
    assert points_axes(tmp_path / "alone", "s") == ["x", "y", "z"]
