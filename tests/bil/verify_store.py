#!/usr/bin/env python3
"""Check a finished store the way a sceptic would, from outside.

  python tests/bil/verify_store.py <sample_dir> [--sample 20000]

<sample_dir> is what the pipeline wrote for one dataset — the folder holding
`combined_output/`, `mmc_output*/` and `store/` — or the store itself, in
which case the cross-checks against the inputs are skipped.

It answers the three questions that matter, and nothing else:

  1. does the store pass the gate the viewer applies?
  2. is every molecule of every slice there, and no blanks?
  3. did each slice's molecules land on that slice's cells, and did the
     cell-type labels land on the right cells?

Exit 0 when every check passes, 1 when any fails. Each line is OK or FAIL,
so this is usable in a loop over a fleet.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

FAILED: list[str] = []


def check(ok: bool, line: str) -> None:
    print(f"{'OK  ' if ok else 'FAIL'} {line}")
    if not ok:
        FAILED.append(line)


def find_store(folder: Path) -> Path:
    if (folder / "zarr.json").exists():
        return folder
    for d in (folder, folder / "store"):
        stores = sorted(d.glob("*.zarr")) if d.is_dir() else []
        if stores:
            return stores[0]
    raise SystemExit(f"{folder}: no .zarr here, nor in store/")


def source_rows(path: Path) -> int:
    """Rows of a transcripts CSV, without reading it into memory."""
    out = subprocess.run(["wc", "-l", str(path)], capture_output=True, text=True)
    return int(out.stdout.split()[0]) - 1 if out.returncode == 0 else -1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("--sample", type=int, default=20_000, help="rows the validator checks per level")
    ap.add_argument("--skip-source", action="store_true", help="do not re-count the source CSVs")
    args = ap.parse_args(argv)

    import numpy as np
    import pandas as pd
    import points_zarr as pz
    import validate_store
    import anndata as ad

    folder = args.folder.resolve()
    store = find_store(folder)
    combined = folder / "combined_output"
    print(f"store: {store}\n")

    # ── 1. the gate ──────────────────────────────────────────
    mapping = store.parent / "mapping.json"
    report = validate_store.validate(store, sample=args.sample)
    check(not report.failed, f"validator: {'passed' if not report.failed else 'FAILED'}")
    for level, _, message in report.lines:
        if level == "FAIL":
            print(f"       {message}")
    if mapping.exists():
        check(True, f"mapping.json: {json.loads(mapping.read_text())}")

    stats_path = store.parent / "stats.json"
    stats = json.loads(stats_path.read_text()) if stats_path.exists() else {}
    if stats:
        check(stats.get("molecules", 0) > 0 or not (store / "points").exists(),
              f"stats.json: {stats['cells']:,} cells, {stats['genes']} genes, "
              f"{stats.get('molecules', 0):,} molecules")

    table = ad.read_zarr(store / "tables" / "table")
    print(f"\ncells: {table.n_obs:,} x {table.n_vars} genes")

    element = store / "points" / "molecules"
    if not element.exists():
        print("\n(no points element — cells only)")
        return 1 if FAILED else 0

    _, attrs = pz.open_element(element)
    total = int(attrs["points:n_points"])
    levels = attrs["multiscales"]["layout"]
    print(f"molecules: {total:,} in {len(levels)} levels, "
          f"level 0 tile {levels[0]['points:grid']['tile_size']:g} um\n")

    # ── 2. every molecule of every slice ─────────────────────
    frame = pz.read_dataframe(element, level=0, order="tile_major",
                              columns=[c for c in ["x", "y", "sample", "unassigned"]
                                       if c in list(pz.open_element(element)[0]["0"]["tile_major"]["columns"][:])])
    per_slice = frame["sample"].value_counts() if "sample" in frame else None

    offsets = {}
    if (combined / "sample_offsets.json").exists():
        offsets = {s["sample_id"]: s
                   for s in json.loads((combined / "sample_offsets.json").read_text())["samples"]}
        check(len(offsets) > 0, f"sample_offsets.json: {len(offsets)} slice(s)")

    if offsets and not args.skip_source:
        counted = 0
        for sid, off in offsets.items():
            if off.get("transcripts"):
                n = source_rows(Path(off["transcripts"]))
                counted += max(n, 0)
        if counted > 0:
            dropped = counted - total
            check(0 <= dropped < counted * 0.2,
                  f"molecules: {total:,} written of {counted:,} source rows "
                  f"({dropped:,} blanks/controls dropped)")

    # ── 3. molecules on their own cells, labels on their own cells ──
    if per_slice is not None and "_sample_id" in table.obs:
        print()
        for sid in sorted(per_slice.index.astype(str)):
            mine = frame[frame["sample"].astype(str) == sid]
            cells = table.obs[table.obs["_sample_id"].astype(str) == sid]
            if len(cells) == 0:
                check(False, f"{sid}: {len(mine):,} molecules but no cells with that _sample_id")
                continue
            inside = (
                (mine["x"] >= cells["center_x"].min() - 100) & (mine["x"] <= cells["center_x"].max() + 100) &
                (mine["y"] >= cells["center_y"].min() - 100) & (mine["y"] <= cells["center_y"].max() + 100)
            ).mean()
            check(inside > 0.99,
                  f"{sid}: {len(mine):>12,} molecules, {len(cells):>8,} cells, "
                  f"{inside * 100:5.1f}% of molecules on those cells")

    if "unassigned" in frame:
        u = np.asarray(frame["unassigned"])
        print(f"\nin no cell: {u.sum():,} of {len(u):,} ({u.mean() * 100:.1f}%)")

    # the labels, joined the only way that is safe
    mmc = next((p / "mapping_output.csv" for p in sorted(folder.glob("mmc_output*"))
                if (p / "mapping_output.csv").exists()), None)
    if mmc is not None:
        labels = pd.read_csv(mmc, comment="#", dtype=str)
        column = next((c for c in ("class_name", "cluster_name", "supercluster_name")
                       if c in labels.columns), None)
        # The table's index is what identifies a cell: <_sample_id>_<id> for a
        # combined dataset, the id itself otherwise. MapMyCells keys by the
        # same thing, so the index is the join — never a rebuilt string.
        keys = [np.asarray(table.obs_names, dtype=str)]
        if "_sample_id" in table.obs and "slice_cell_id" in table.obs:
            keys.append((table.obs["_sample_id"].astype(str) + "_"
                         + table.obs["slice_cell_id"].astype(str)).to_numpy())
        if column:
            by_id = labels.set_index("cell_id")[column]
            key = next((k for k in keys if pd.Index(k).isin(by_id.index).all()), None)
            check(key is not None, f"{mmc.parent.name}: every cell found by <sample>_<cell_id>")
            if key is not None and column in table.obs:
                truth = by_id.reindex(key)
                check(bool((truth.to_numpy() == table.obs[column].astype(str).to_numpy()).all()),
                      f"{mmc.parent.name}: obs/{column} equals the mapping for every cell")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("every check passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
