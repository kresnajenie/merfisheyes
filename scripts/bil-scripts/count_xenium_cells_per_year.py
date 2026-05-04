#!/usr/bin/env python3
"""
count_xenium_cells_per_year.py

Like count_cells_per_year.py but only for Xenium datasets (rows in the
spreadsheet whose Notes column == "xenium", case-insensitive).

For each Xenium dataset:
  - publication_date = mtime of the BIL data path (same approach as the
    MERFISH script).
  - n_cells = number of UNIQUE cell_id values in the Xenium
    `output-XETG*selected_transcripts.csv` file inside the BIL path.
    (The transcripts CSV has one row per detected molecule, so cell_ids
    repeat — we want the deduplicated count.)

Outputs (in --out-dir, default /bil/data/meyes/_cells_inventory/):
  cells_per_year_xenium.csv          year, n_cells
  cells_inventory_xenium.csv         per-dataset detail (transcripts file path, etc.)
  cells_per_year_xenium_manifest.csv per-dataset included/excluded with reasons

Usage:
  python count_xenium_cells_per_year.py
"""

import argparse
import csv
import gzip
import logging
import os
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("count_xenium")

DEFAULT_DATASETS = "scripts/bil-scripts/all_datasets.csv"
DEFAULT_OUT = "/bil/data/meyes/_cells_inventory"

# How deep to descend looking for the transcripts CSV. Xenium folders
# usually keep it at the top level or one subdir down.
MAX_WALK_DEPTH = 3
DEFAULT_WALK_TIMEOUT = 60.0
DEFAULT_MAX_WALK_FILES = 5000

# Cell IDs Xenium uses to mean "no cell" — exclude from unique count.
NO_CELL_TOKENS = {"", "unassigned", "-1", "nan", "none"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datasets-csv", default=DEFAULT_DATASETS,
                   help=f"Spreadsheet with Notes + 'bil data path' columns "
                        f"(default: {DEFAULT_DATASETS}).")
    p.add_argument("--out-dir", default=DEFAULT_OUT,
                   help=f"Output directory (default: {DEFAULT_OUT}).")
    p.add_argument("--max-walk-depth", type=int, default=MAX_WALK_DEPTH,
                   help=f"Subdir levels to descend looking for the transcripts "
                        f"CSV (default: {MAX_WALK_DEPTH}).")
    p.add_argument("--walk-timeout", type=float, default=DEFAULT_WALK_TIMEOUT,
                   help=f"Per-dataset cap on the walk, in seconds "
                        f"(default: {DEFAULT_WALK_TIMEOUT}).")
    p.add_argument("--max-walk-files", type=int, default=DEFAULT_MAX_WALK_FILES,
                   help=f"Per-dataset cap on files inspected during the walk "
                        f"(default: {DEFAULT_MAX_WALK_FILES}).")
    return p.parse_args()


def read_dataset_rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return [
            {(k or "").strip().lower(): (v or "").strip() for k, v in r.items()}
            for r in reader
        ]


def is_xenium(row: dict) -> bool:
    return row.get("notes", "").strip().lower() == "xenium"


def path_publication_date(path_str: str) -> str:
    if not path_str:
        return ""
    p = Path(path_str.rstrip("/"))
    try:
        st = os.stat(p)
    except (PermissionError, FileNotFoundError):
        return ""
    except OSError as e:
        logger.warning("stat failed for %s: %s", p, e)
        return ""
    return datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d")


def is_xenium_transcripts(name: str) -> bool:
    """Match output-XETG*selected_transcripts.csv (.gz optional)."""
    low = name.lower()
    if not low.startswith("output-xetg"):
        return False
    return low.endswith("selected_transcripts.csv") or \
        low.endswith("selected_transcripts.csv.gz")


def find_transcripts_csv(root: Path, max_depth: int,
                          deadline: float, max_files: int) -> tuple[Path | None, str]:
    """BFS the BIL path looking for the Xenium transcripts CSV. Returns
    (path, bail_reason). bail_reason is '' on success."""
    try:
        if not root.exists():
            return None, "missing_path"
    except PermissionError:
        return None, "permission_denied"
    except OSError as e:
        return None, f"stat_error:{e.errno}"

    queue: list[tuple[Path, int]] = [(root, 0)]
    files_seen = 0
    while queue:
        if time.monotonic() > deadline:
            return None, "walk_timeout"
        d, depth = queue.pop(0)
        try:
            entries = list(d.iterdir())
        except (PermissionError, OSError):
            continue
        for entry in entries:
            try:
                if entry.is_file():
                    files_seen += 1
                    if is_xenium_transcripts(entry.name):
                        return entry, ""
                    if files_seen >= max_files:
                        return None, "walk_max_files"
                    if time.monotonic() > deadline:
                        return None, "walk_timeout"
                elif entry.is_dir() and depth < max_depth:
                    queue.append((entry, depth + 1))
            except OSError:
                continue
    return None, "no_transcripts_csv"


def count_unique_cell_ids(transcripts_path: Path) -> int:
    """Stream the transcripts CSV and return len(set(cell_id))."""
    opener = gzip.open if transcripts_path.suffix == ".gz" else open
    seen: set[str] = set()
    try:
        with opener(transcripts_path, "rt", errors="replace") as fh:
            reader = csv.reader(fh)
            header = next(reader, None)
            if header is None:
                return 0
            try:
                idx = header.index("cell_id")
            except ValueError:
                # Try common variants.
                col_lookup = {h.strip().lower(): i for i, h in enumerate(header)}
                idx = col_lookup.get("cell_id")
                if idx is None:
                    idx = col_lookup.get("cellid")
                if idx is None:
                    logger.warning("No cell_id column in %s (header: %s)",
                                   transcripts_path, header[:5])
                    return -1
            for row in reader:
                if idx >= len(row):
                    continue
                cid = row[idx].strip()
                if cid.lower() in NO_CELL_TOKENS:
                    continue
                seen.add(cid)
    except OSError as e:
        logger.warning("Failed to read %s: %s", transcripts_path, e)
        return -1
    return len(seen)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        stream=sys.stdout)

    datasets_csv = Path(args.datasets_csv)
    if not datasets_csv.exists():
        raise SystemExit(f"Datasets CSV not found: {datasets_csv}")

    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_dataset_rows(datasets_csv)
    xenium_rows = [r for r in rows if is_xenium(r)]
    logger.info("Loaded %d datasets (%d xenium) from %s",
                len(rows), len(xenium_rows), datasets_csv)

    inventory: list[dict] = []
    for r in xenium_rows:
        sample = r.get("dataset", "")
        bil_path = r.get("bil data path", "")
        out: dict = {
            "dataset": sample,
            "status": r.get("status", ""),
            "notes": r.get("notes", ""),
            "bil_path": bil_path,
            "publication_date": "",
            "year": "",
            "n_cells": "",
            "transcripts_path": "",
            "source": "",
        }

        date = path_publication_date(bil_path)
        out["publication_date"] = date
        out["year"] = date[:4] if date else ""

        if not bil_path:
            out["source"] = "no_bil_path"
            inventory.append(out)
            logger.info("[%s] date=— n_cells=— source=no_bil_path",
                        sample)
            continue

        deadline = time.monotonic() + args.walk_timeout
        transcripts_path, bail = find_transcripts_csv(
            Path(bil_path.rstrip("/")), args.max_walk_depth,
            deadline, args.max_walk_files,
        )
        if transcripts_path is None:
            out["source"] = bail or "no_transcripts_csv"
            inventory.append(out)
            logger.info("[%s] date=%s n_cells=— source=%s",
                        sample, out["publication_date"] or "—", out["source"])
            continue

        n = count_unique_cell_ids(transcripts_path)
        out["transcripts_path"] = str(transcripts_path)
        if n >= 0:
            out["n_cells"] = n
            out["source"] = "xenium_transcripts_unique_cell_id"
        else:
            out["source"] = "transcripts_read_failed"

        inventory.append(out)
        logger.info("[%s] date=%s n_cells=%s source=%s",
                    sample, out["publication_date"] or "—",
                    out["n_cells"] if out["n_cells"] != "" else "—",
                    out["source"])

    # ── Inventory ─────────────────────────────────────────────────
    inv_path = out_dir / "cells_inventory_xenium.csv"
    with open(inv_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["dataset", "status", "notes", "bil_path",
                        "publication_date", "year", "n_cells",
                        "transcripts_path", "source"],
        )
        writer.writeheader()
        writer.writerows(inventory)
    logger.info("Wrote %s (%d rows)", inv_path, len(inventory))

    # ── Aggregate by year + manifest ──────────────────────────────
    cells_by_year: dict[str, int] = defaultdict(int)
    manifest_rows: list[dict] = []
    n_included = n_excluded = 0
    for r in inventory:
        year = r["year"]
        n_cells = r["n_cells"] if isinstance(r["n_cells"], int) else None

        reasons = []
        if not year:
            reasons.append("no_publication_date")
        if n_cells is None:
            reasons.append(f"no_cell_count:{r['source'] or 'unknown'}")
        included = not reasons

        if included:
            cells_by_year[year] += n_cells  # type: ignore[arg-type]
            n_included += 1
        else:
            n_excluded += 1

        manifest_rows.append({
            "dataset": r["dataset"],
            "included": "true" if included else "false",
            "year": year,
            "n_cells": n_cells if n_cells is not None else "",
            "reason": "; ".join(reasons),
        })

    yearly_path = out_dir / "cells_per_year_xenium.csv"
    with open(yearly_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["year", "n_cells"])
        for y in sorted(cells_by_year):
            writer.writerow([y, cells_by_year[y]])
    logger.info("Wrote %s (%d years)", yearly_path, len(cells_by_year))

    manifest_path = out_dir / "cells_per_year_xenium_manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["dataset", "included", "year", "n_cells", "reason"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    logger.info("Wrote %s (%d included, %d excluded)",
                manifest_path, n_included, n_excluded)


if __name__ == "__main__":
    main()
