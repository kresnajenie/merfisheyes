#!/usr/bin/env python3
"""
count_cells_per_year.py

Walk every dataset listed in the all-datasets spreadsheet, recover its
publication date from the BIL data path's mtime, and count the total
number of cells in the dataset.

Cell count source priority (per dataset):
  1. /bil/data/meyes/<dataset>/combined_output/cell_metadata.csv (post-pipeline)
  2. cell_metadata.csv anywhere under the BIL data path        (MERSCOPE-style)
  3. cells.csv / cells.csv.gz / cells.parquet under BIL path   (Xenium-style)
  4. cell_by_gene.h5ad / *.h5ad under BIL path                  (h5ad pipeline)

Outputs (in --out-dir):
  cells_inventory.csv     one row per dataset with date, year, n_cells, source
  cells_per_year.csv      yearly aggregation: year, n_datasets, n_cells

Usage:
  python count_cells_per_year.py                                # defaults
  python count_cells_per_year.py --datasets-csv path.csv \
      --out-dir /bil/data/meyes/_cells_inventory
"""

import argparse
import csv
import gzip
import logging
import os
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("count_cells")

DEFAULT_DATASETS = "scripts/bil-scripts/all_datasets.csv"
DEFAULT_MEYES_BASE = "/bil/data/meyes"
DEFAULT_OUT = "/bil/data/meyes/_cells_inventory"

# Limit recursion depth when scanning a BIL data dir.
MAX_WALK_DEPTH = 4


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datasets-csv", default=DEFAULT_DATASETS,
                   help=f"Spreadsheet with Dataset + 'bil data path' columns "
                        f"(default: {DEFAULT_DATASETS}).")
    p.add_argument("--meyes-base", default=DEFAULT_MEYES_BASE,
                   help=f"Where meyes processing outputs live "
                        f"(default: {DEFAULT_MEYES_BASE}).")
    p.add_argument("--out-dir", default=DEFAULT_OUT,
                   help=f"Output directory (default: {DEFAULT_OUT}).")
    p.add_argument("--max-walk-depth", type=int, default=MAX_WALK_DEPTH,
                   help=f"How many subdir levels to scan under each BIL "
                        f"path when looking for cell files "
                        f"(default: {MAX_WALK_DEPTH}).")
    return p.parse_args()


def read_dataset_rows(path: Path) -> list[dict]:
    """Read the spreadsheet and return list of dicts with normalized keys."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for r in reader:
            # Normalize keys (strip and lowercase) so 'bil data path' becomes
            # easy to address.
            norm = {(k or "").strip().lower(): (v or "").strip()
                    for k, v in r.items()}
            rows.append(norm)
    return rows


def path_publication_date(path_str: str) -> str:
    if not path_str:
        return ""
    p = Path(path_str.rstrip("/"))
    try:
        st = os.stat(p)
    except (OSError, FileNotFoundError) as e:
        logger.warning("stat failed for %s: %s", p, e)
        return ""
    return datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d")


def count_csv_rows(path: Path) -> int:
    """Count data rows (subtract one for header). Handles .gz transparently."""
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt", errors="replace") as fh:
        for _ in fh:
            n += 1
    return max(n - 1, 0)


def count_parquet_rows(path: Path) -> int:
    try:
        import pyarrow.parquet as pq
    except ImportError:
        logger.warning("pyarrow not available — skipping parquet at %s", path)
        return -1
    return int(pq.ParquetFile(str(path)).metadata.num_rows)


def count_h5ad_obs(path: Path) -> int:
    """Read only obs index dim from h5ad (cheap; uses h5py)."""
    try:
        import h5py
    except ImportError:
        logger.warning("h5py not available — skipping h5ad at %s", path)
        return -1
    try:
        with h5py.File(str(path), "r") as f:
            # AnnData stores obs as a group; row count == len of any obs col,
            # or the _index dataset if present.
            obs = f.get("obs")
            if obs is None:
                return -1
            if "_index" in obs:
                return int(obs["_index"].shape[0])
            for k in obs.keys():
                v = obs[k]
                if hasattr(v, "shape") and v.shape:
                    return int(v.shape[0])
            return -1
    except Exception as e:
        logger.warning("Failed to read h5ad %s: %s", path, e)
        return -1


def walk_files(root: Path, max_depth: int):
    """BFS yielding files up to max_depth, skipping permission errors."""
    queue: list[tuple[Path, int]] = [(root, 0)]
    while queue:
        d, depth = queue.pop(0)
        try:
            entries = list(d.iterdir())
        except (PermissionError, OSError):
            continue
        for entry in entries:
            try:
                if entry.is_file():
                    yield entry
                elif entry.is_dir() and depth < max_depth:
                    queue.append((entry, depth + 1))
            except OSError:
                continue


# Cell-file lookup: ordered list of (matcher, counter, label_template).
CELL_FILE_PATTERNS: list[tuple] = [
    (lambda n: n == "cell_metadata.csv",
        count_csv_rows, "merscope_cell_metadata"),
    (lambda n: n == "cells.csv" or n == "cells.csv.gz",
        count_csv_rows, "xenium_cells_csv"),
    (lambda n: n == "cells.parquet",
        count_parquet_rows, "xenium_cells_parquet"),
    (lambda n: n == "cell_by_gene.h5ad",
        count_h5ad_obs, "h5ad_cell_by_gene"),
    (lambda n: n.endswith(".h5ad"),
        count_h5ad_obs, "h5ad_other"),
]


def count_cells_in_path(root: Path, max_depth: int) -> tuple[int, str, str]:
    """Walk root looking for the first recognized cell file. Returns
    (n_cells, source_label, relative_path)."""
    if not root.exists():
        return -1, "missing_path", ""
    found_by_pattern: dict[int, tuple[Path, int]] = {}
    for f in walk_files(root, max_depth):
        for idx, (matcher, _, _) in enumerate(CELL_FILE_PATTERNS):
            if matcher(f.name):
                # Keep the first hit per pattern; preserve discovery order.
                if idx not in found_by_pattern:
                    found_by_pattern[idx] = (f, idx)
                break
    for idx in range(len(CELL_FILE_PATTERNS)):
        if idx not in found_by_pattern:
            continue
        path = found_by_pattern[idx][0]
        _, counter, label = CELL_FILE_PATTERNS[idx]
        n = counter(path)
        if n >= 0:
            try:
                rel = str(path.relative_to(root))
            except ValueError:
                rel = str(path)
            return n, label, rel
    return -1, "no_recognized_files", ""


def count_cells_in_meyes(meyes_base: Path, sample: str) -> tuple[int, str, str]:
    """Try the canonical post-pipeline file first."""
    candidate = meyes_base / sample / "combined_output" / "cell_metadata.csv"
    if not candidate.exists():
        return -1, "no_meyes_metadata", ""
    n = count_csv_rows(candidate)
    return n, "meyes_combined_output", str(candidate)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        stream=sys.stdout)

    datasets_csv = Path(args.datasets_csv)
    if not datasets_csv.exists():
        raise SystemExit(f"Datasets CSV not found: {datasets_csv}")

    meyes_base = Path(args.meyes_base).expanduser()
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = read_dataset_rows(datasets_csv)
    logger.info("Loaded %d datasets from %s", len(rows), datasets_csv)

    inventory: list[dict] = []
    for r in rows:
        sample = r.get("dataset", "")
        if not sample:
            continue
        bil_path = r.get("bil data path", "")
        status = r.get("status", "")
        notes = r.get("notes", "")

        out: dict = {
            "dataset": sample,
            "status": status,
            "notes": notes,
            "bil_path": bil_path,
            "publication_date": "",
            "year": "",
            "n_cells": "",
            "cell_source": "",
            "cell_source_path": "",
            "in_meyes": "",
        }

        date = path_publication_date(bil_path) if bil_path else ""
        out["publication_date"] = date
        out["year"] = date[:4] if date else ""

        # 1. Try meyes first (most reliable post-pipeline count).
        n, source, source_path = count_cells_in_meyes(meyes_base, sample)
        out["in_meyes"] = "true" if n >= 0 else "false"

        # 2. Fall back to the BIL path itself.
        if n < 0 and bil_path:
            n, source, source_path = count_cells_in_path(
                Path(bil_path.rstrip("/")), args.max_walk_depth,
            )

        if n >= 0:
            out["n_cells"] = n
            out["cell_source"] = source
            out["cell_source_path"] = source_path
        else:
            out["cell_source"] = source  # explanation: missing_path / no_recognized_files / no_meyes_metadata

        inventory.append(out)
        logger.info("[%s] date=%s n_cells=%s source=%s",
                    sample, out["publication_date"] or "—",
                    out["n_cells"] if out["n_cells"] != "" else "—",
                    out["cell_source"] or "—")

    # ── Write per-dataset inventory ───────────────────────────────
    inv_path = out_dir / "cells_inventory.csv"
    fieldnames = ["dataset", "status", "notes", "bil_path",
                  "publication_date", "year", "n_cells",
                  "cell_source", "cell_source_path", "in_meyes"]
    with open(inv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(inventory)
    logger.info("Wrote %s (%d rows)", inv_path, len(inventory))

    # ── Aggregate by year ─────────────────────────────────────────
    cells_by_year: dict[str, int] = defaultdict(int)
    for r in inventory:
        year = r["year"]
        if not year or not isinstance(r["n_cells"], int):
            continue
        cells_by_year[year] += r["n_cells"]

    yearly_path = out_dir / "cells_per_year.csv"
    with open(yearly_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["year", "n_cells"])
        for y in sorted(cells_by_year):
            writer.writerow([y, cells_by_year[y]])
    logger.info("Wrote %s (%d years)", yearly_path, len(cells_by_year))


if __name__ == "__main__":
    main()
