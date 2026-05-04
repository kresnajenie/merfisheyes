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
import re
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("count_cells")

DEFAULT_DATASETS = "scripts/bil-scripts/all_datasets.csv"
DEFAULT_MEYES_BASE = "/bil/data/meyes"
DEFAULT_OUT = "/bil/data/meyes/_cells_inventory"

# Limit recursion depth when scanning a BIL data dir.
MAX_WALK_DEPTH = 4
# Recognized artifact-mask filenames inside combined_output/.
MASK_NAME_RE = re.compile(r"^artifact_mask_p(\d+)(?:\.csv)?$", re.IGNORECASE)
# Column-name patterns for the artifact flag inside a mask CSV.
ARTIFACT_COL_NAMES = ("is_artifact", "artifact", "artifact_mask",
                      "mask", "p25", "is_mask")
# Hard caps so a giant BIL tree (e.g. raw transcript tiles, image stacks)
# can't stall the whole scrape on a single dataset.
DEFAULT_WALK_TIMEOUT = 60.0      # seconds per dataset
DEFAULT_MAX_WALK_FILES = 20000   # files inspected per dataset


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
    p.add_argument("--walk-timeout", type=float, default=DEFAULT_WALK_TIMEOUT,
                   help=f"Per-dataset wall-clock cap on the BIL fallback "
                        f"walk, in seconds (default: {DEFAULT_WALK_TIMEOUT}).")
    p.add_argument("--max-walk-files", type=int, default=DEFAULT_MAX_WALK_FILES,
                   help=f"Per-dataset cap on files inspected during the BIL "
                        f"fallback walk (default: {DEFAULT_MAX_WALK_FILES}).")
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


class WalkBudget:
    """Tracks per-dataset walk limits and why we stopped (if we stopped)."""
    def __init__(self, max_files: int, deadline: float):
        self.max_files = max_files
        self.deadline = deadline
        self.files_seen = 0
        self.bail_reason: str | None = None

    def file_seen(self) -> bool:
        self.files_seen += 1
        if self.max_files is not None and self.files_seen >= self.max_files:
            self.bail_reason = "walk_max_files"
            return False
        if time.monotonic() > self.deadline:
            self.bail_reason = "walk_timeout"
            return False
        return True

    def time_left(self) -> bool:
        if time.monotonic() > self.deadline:
            self.bail_reason = "walk_timeout"
            return False
        return True


def walk_files(root: Path, max_depth: int, budget: WalkBudget):
    """BFS yielding files up to max_depth, skipping permission errors and
    bailing once the budget is exhausted."""
    queue: list[tuple[Path, int]] = [(root, 0)]
    while queue:
        if not budget.time_left():
            return
        d, depth = queue.pop(0)
        try:
            entries = list(d.iterdir())
        except (PermissionError, OSError):
            continue
        for entry in entries:
            try:
                if entry.is_file():
                    yield entry
                    if not budget.file_seen():
                        return
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


def count_cells_in_path(root: Path, max_depth: int,
                         walk_timeout: float,
                         max_walk_files: int) -> tuple[int, str, str]:
    """Walk root looking for the first recognized cell file. Returns
    (n_cells, source_label, relative_path). Bails out and labels the source
    as walk_timeout / walk_max_files if the budget is exhausted before any
    recognized file is found."""
    try:
        if not root.exists():
            return -1, "missing_path", ""
    except PermissionError:
        return -1, "permission_denied", ""
    except OSError as e:
        return -1, f"stat_error:{e.errno}", ""
    budget = WalkBudget(max_files=max_walk_files,
                         deadline=time.monotonic() + walk_timeout)
    found_by_pattern: dict[int, Path] = {}
    for f in walk_files(root, max_depth, budget):
        for idx, (matcher, _, _) in enumerate(CELL_FILE_PATTERNS):
            if matcher(f.name):
                if idx not in found_by_pattern:
                    found_by_pattern[idx] = f
                # Highest-priority match found — skip the rest of the walk.
                if idx == 0:
                    break
                break
        if 0 in found_by_pattern:
            break
    for idx in range(len(CELL_FILE_PATTERNS)):
        if idx not in found_by_pattern:
            continue
        path = found_by_pattern[idx]
        _, counter, label = CELL_FILE_PATTERNS[idx]
        n = counter(path)
        if n >= 0:
            try:
                rel = str(path.relative_to(root))
            except ValueError:
                rel = str(path)
            return n, label, rel
    if budget.bail_reason:
        return -1, budget.bail_reason, ""
    return -1, "no_recognized_files", ""


def count_cells_in_meyes(meyes_base: Path, sample: str) -> tuple[int, str, str]:
    """Try the canonical post-pipeline file first."""
    candidate = meyes_base / sample / "combined_output" / "cell_metadata.csv"
    if not candidate.exists():
        return -1, "no_meyes_metadata", ""
    n = count_csv_rows(candidate)
    return n, "meyes_combined_output", str(candidate)


def find_highest_percentile_mask(combined_dir: Path) -> Path | None:
    """Return the path to the artifact_mask_p<N> entry with the largest N,
    resolved to the actual CSV file (with or without .csv suffix). Returns
    None if no mask is found."""
    if not combined_dir.is_dir():
        return None
    best: tuple[int, Path] | None = None
    try:
        entries = list(combined_dir.iterdir())
    except (PermissionError, OSError):
        return None
    for entry in entries:
        m = MASK_NAME_RE.match(entry.name)
        if not m:
            continue
        n = int(m.group(1))
        if best is None or n > best[0]:
            best = (n, entry)
    if best is None:
        return None
    candidate = best[1]
    # Resolve to a readable CSV file.
    if candidate.is_file():
        return candidate
    if candidate.is_dir():
        for child in candidate.glob("*.csv"):
            if child.is_file():
                return child
        return None
    with_csv = combined_dir / f"{candidate.name}.csv"
    if with_csv.is_file():
        return with_csv
    return None


def count_mask_artifacts(mask_path: Path) -> int:
    """Count rows in the mask CSV where the artifact flag is truthy.
    Returns -1 on failure."""
    try:
        with open(mask_path, newline="") as fh:
            reader = csv.reader(fh)
            header = next(reader, None)
            if header is None:
                return 0
            # Pick the artifact column.
            col_idx = None
            for i, name in enumerate(header):
                low = name.strip().lower()
                if low in ARTIFACT_COL_NAMES or low.startswith("artifact_mask"):
                    col_idx = i
                    break
            if col_idx is None and len(header) == 1:
                col_idx = 0
            if col_idx is None:
                # Default to the last column — fine for narrow mask CSVs that
                # store IDs first and the boolean flag last.
                col_idx = len(header) - 1

            n = 0
            for row in reader:
                if col_idx >= len(row):
                    continue
                v = row[col_idx].strip().lower()
                if v in ("1", "true", "t", "yes", "y"):
                    n += 1
                elif v in ("0", "false", "f", "no", "n", ""):
                    continue
                else:
                    try:
                        if float(v) > 0:
                            n += 1
                    except ValueError:
                        pass
            return n
    except OSError as e:
        logger.warning("Failed to read mask %s: %s", mask_path, e)
        return -1


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
            "n_cells_raw": "",
            "n_cells_artifact": "",
            "mask_name": "",
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
                walk_timeout=args.walk_timeout,
                max_walk_files=args.max_walk_files,
            )

        # 3. Mask filter: if combined_output has an artifact_mask_p<N> file,
        #    record the artifact count separately so the final TSV can
        #    report both raw and filtered totals. Skip silently if absent.
        n_filtered = n
        n_artifact = 0
        if n >= 0 and out["in_meyes"] == "true":
            mask_path = find_highest_percentile_mask(
                meyes_base / sample / "combined_output"
            )
            if mask_path is not None:
                got = count_mask_artifacts(mask_path)
                if got >= 0:
                    n_artifact = got
                    out["mask_name"] = mask_path.name
                    n_filtered = max(n - n_artifact, 0)

        if n >= 0:
            out["n_cells_raw"] = n
            out["n_cells_artifact"] = n_artifact if out["mask_name"] else ""
            out["n_cells"] = n_filtered  # raw when no mask, filtered otherwise
            out["cell_source"] = source
            out["cell_source_path"] = source_path
        else:
            out["cell_source"] = source

        logger.info("[%s] date=%s raw=%s filtered=%s mask=%s source=%s",
                    sample, out["publication_date"] or "—",
                    out["n_cells_raw"] if out["n_cells_raw"] != "" else "—",
                    out["n_cells"] if out["n_cells"] != "" else "—",
                    out["mask_name"] or "—",
                    out["cell_source"] or "—")
        inventory.append(out)

    # ── Write per-dataset inventory ───────────────────────────────
    inv_path = out_dir / "cells_inventory.csv"
    fieldnames = ["dataset", "status", "notes", "bil_path",
                  "publication_date", "year", "n_cells",
                  "n_cells_raw", "n_cells_artifact", "mask_name",
                  "cell_source", "cell_source_path", "in_meyes"]
    with open(inv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(inventory)
    logger.info("Wrote %s (%d rows)", inv_path, len(inventory))

    # ── Aggregate by year + write manifest ───────────────────────
    raw_by_year: dict[str, int] = defaultdict(int)
    filtered_by_year: dict[str, int] = defaultdict(int)
    manifest_rows: list[dict] = []
    n_included = 0
    n_excluded = 0
    for r in inventory:
        year = r["year"]
        n_raw = r["n_cells_raw"] if isinstance(r["n_cells_raw"], int) else None
        n_filtered = r["n_cells"] if isinstance(r["n_cells"], int) else None

        reasons = []
        if not year:
            reasons.append("no_publication_date")
        if n_raw is None:
            reasons.append(f"no_cell_count:{r['cell_source'] or 'unknown'}")
        included = not reasons

        if included:
            raw_by_year[year] += n_raw  # type: ignore[arg-type]
            filtered_by_year[year] += n_filtered  # type: ignore[arg-type]
            n_included += 1
        else:
            n_excluded += 1

        manifest_rows.append({
            "dataset": r["dataset"],
            "included": "true" if included else "false",
            "year": year,
            "n_cells_raw": n_raw if n_raw is not None else "",
            "n_cells_filtered": n_filtered if n_filtered is not None else "",
            "mask_name": r["mask_name"] or "",
            "reason": "; ".join(reasons),
        })

    yearly_path = out_dir / "cells_per_year.csv"
    with open(yearly_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["year", "n_cells_raw", "n_cells_filtered"])
        for y in sorted(raw_by_year):
            writer.writerow([y, raw_by_year[y], filtered_by_year[y]])
    logger.info("Wrote %s (%d years)", yearly_path, len(raw_by_year))

    manifest_path = out_dir / "cells_per_year_manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["dataset", "included", "year", "n_cells_raw",
                        "n_cells_filtered", "mask_name", "reason"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    logger.info("Wrote %s (%d included, %d excluded)",
                manifest_path, n_included, n_excluded)


if __name__ == "__main__":
    main()
