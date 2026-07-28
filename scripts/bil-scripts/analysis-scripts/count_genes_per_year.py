#!/usr/bin/env python3
"""
count_genes_per_year.py

Counts the TOTAL NUMBER OF GENES IMAGED per year.

"Genes imaged" here means individual detected-transcript molecules — every
RNA spot the microscope called, NOT the size of the gene panel. A 300-gene
MERFISH run that detects 67 million transcripts contributes 67,000,000, not
300. Per-year totals therefore land in the hundreds of millions / billions.

For each dataset in the all-datasets spreadsheet:
  - publication_date = mtime of the BIL data path.
  - n_genes_imaged = total detected molecules, found from (in priority order):
      1. The single-molecule processed manifest the meyes pipeline writes
         (`manifest.json[.gz]` whose statistics carry `total_molecules`),
         searched first under <meyes-base>/<dataset>/ then the BIL path.
      2. Fallback: a streamed row count of a raw molecule table in the BIL
         path — MERSCOPE `detected_transcripts.csv[.gz]` or Xenium
         `transcripts.csv[.gz]` (one row = one detected molecule). Parquet
         molecule tables are counted via pyarrow metadata when available.

Outputs (in --out-dir):
  genes_inventory.csv             one row per dataset: date, year,
                                  n_genes_imaged, source, source path
  genes_per_year.csv              year, n_datasets, total_genes_imaged
                                  (plus a final TOTAL row across all years)
  genes_per_year_manifest.csv     per-dataset included/excluded with reasons

total_genes_imaged = sum of per-dataset detected-molecule counts for that
year. There is no "unique genes" column anymore — a molecule total is a
plain sum, the union of panels is a different question and is not what this
script answers.

Note (multi-region datasets): like count_xenium_cells_per_year.py, a single
spreadsheet row resolves to a single best source. If a dataset has several
independent single-molecule manifests, the largest is used (see
--sum-regions to add them instead).

Usage:
  python count_genes_per_year.py                          # defaults
  python count_genes_per_year.py --datasets-csv path.csv \
      --out-dir /bil/data/meyes/_genes_inventory
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import logging
import os
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

logger = logging.getLogger("count_genes")

DEFAULT_DATASETS = "scripts/bil-scripts/analysis-scripts/all_datasets.csv"
DEFAULT_MEYES_BASE = "/bil/data/meyes"
DEFAULT_OUT = "/bil/data/meyes/_genes_inventory"

MAX_WALK_DEPTH = 4
DEFAULT_WALK_TIMEOUT = 60.0
DEFAULT_MAX_WALK_FILES = 20000

# Raw molecule tables can have huge headers / fields.
csv.field_size_limit(sys.maxsize)


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
                   help=f"How many subdir levels to scan under each search "
                        f"root (default: {MAX_WALK_DEPTH}).")
    p.add_argument("--walk-timeout", type=float, default=DEFAULT_WALK_TIMEOUT,
                   help=f"Per-search wall-clock cap, in seconds "
                        f"(default: {DEFAULT_WALK_TIMEOUT}).")
    p.add_argument("--max-walk-files", type=int, default=DEFAULT_MAX_WALK_FILES,
                   help=f"Per-search cap on files inspected "
                        f"(default: {DEFAULT_MAX_WALK_FILES}).")
    p.add_argument("--sum-regions", action="store_true",
                   help="If a dataset has multiple single-molecule manifests, "
                        "sum their total_molecules instead of taking the "
                        "largest. Off by default to avoid double-counting "
                        "overlapping exports of the same molecules.")
    return p.parse_args()


def read_dataset_rows(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return [
            {(k or "").strip().lower(): (v or "").strip() for k, v in r.items()}
            for r in reader
        ]


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


# ── Budgeted directory walk ───────────────────────────────────────

class WalkBudget:
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


# ── Source 1: single-molecule manifest ────────────────────────────

def read_manifest_total_molecules(path: Path) -> int | None:
    """Return statistics.total_molecules from a (optionally gzipped) manifest
    JSON, or None if the file is not a molecule manifest / unreadable."""
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            data = json.load(fh)
    except (OSError, ValueError) as e:
        logger.debug("manifest unreadable %s: %s", path, e)
        return None
    if not isinstance(data, dict):
        return None
    stats = data.get("statistics")
    if not isinstance(stats, dict):
        return None
    tm = stats.get("total_molecules")
    if isinstance(tm, bool) or not isinstance(tm, (int, float)):
        return None
    tm = int(tm)
    return tm if tm >= 0 else None


def find_molecule_manifests(root: Path, max_depth: int,
                            walk_timeout: float,
                            max_walk_files: int,
                            ) -> tuple[list[tuple[Path, int]], str]:
    """Walk `root` collecting (path, total_molecules) for every
    single-molecule manifest. Returns (hits, bail_reason)."""
    try:
        if not root.exists():
            return [], "missing_path"
    except PermissionError:
        return [], "permission_denied"
    except OSError as e:
        return [], f"stat_error:{e.errno}"

    budget = WalkBudget(max_files=max_walk_files,
                        deadline=time.monotonic() + walk_timeout)
    hits: list[tuple[Path, int]] = []
    for f in walk_files(root, max_depth, budget):
        name = f.name.lower()
        if name not in ("manifest.json", "manifest.json.gz"):
            continue
        tm = read_manifest_total_molecules(f)
        if tm is not None:
            hits.append((f, tm))
    return hits, (budget.bail_reason or "")


def pick_manifest_total(hits: list[tuple[Path, int]],
                        sum_regions: bool) -> tuple[int, str]:
    """Collapse manifest hits to one number. De-duplicates exact-equal
    totals (same molecules exported under two paths) before summing."""
    if not sum_regions:
        path, tm = max(hits, key=lambda h: h[1])
        return tm, str(path)
    seen: set[int] = set()
    total = 0
    used: list[str] = []
    for path, tm in sorted(hits, key=lambda h: -h[1]):
        if tm in seen:
            continue
        seen.add(tm)
        total += tm
        used.append(str(path))
    return total, "; ".join(used)


# ── Source 2: raw molecule table row count ────────────────────────

RAW_MOLECULE_TABLES = (
    "detected_transcripts.csv",
    "detected_transcripts.csv.gz",
    "transcripts.csv",
    "transcripts.csv.gz",
)


def count_csv_data_rows(path: Path) -> int:
    """Stream a molecule CSV and count data rows (one row = one molecule)."""
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            reader = csv.reader(fh)
            if next(reader, None) is None:        # header
                return 0
            n = 0
            for _ in reader:
                n += 1
            return n
    except OSError as e:
        logger.warning("Failed to read %s: %s", path, e)
        return -1


def count_parquet_rows(path: Path) -> int:
    """Row count from parquet footer metadata (no full read). -1 if pyarrow
    is unavailable or the file is unreadable."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return -1
    try:
        return pq.ParquetFile(str(path)).metadata.num_rows
    except Exception as e:
        logger.warning("Failed to read parquet %s: %s", path, e)
        return -1


def molecules_from_raw_tables(root: Path, max_depth: int,
                              walk_timeout: float,
                              max_walk_files: int,
                              ) -> tuple[int, str, str]:
    """Walk `root` for the first raw molecule table and count its rows.
    Returns (n_molecules, source_label, source_path); n_molecules < 0 means
    nothing usable was found (source_label carries the reason)."""
    try:
        if not root.exists():
            return -1, "missing_path", ""
    except PermissionError:
        return -1, "permission_denied", ""
    except OSError as e:
        return -1, f"stat_error:{e.errno}", ""

    budget = WalkBudget(max_files=max_walk_files,
                        deadline=time.monotonic() + walk_timeout)
    csv_hit: Path | None = None
    parquet_hit: Path | None = None
    for f in walk_files(root, max_depth, budget):
        low = f.name.lower()
        if csv_hit is None and low in RAW_MOLECULE_TABLES:
            csv_hit = f
            break
        if parquet_hit is None and low in ("detected_transcripts.parquet",
                                           "transcripts.parquet"):
            parquet_hit = f

    if csv_hit is not None:
        n = count_csv_data_rows(csv_hit)
        if n >= 0:
            return n, "raw_molecule_csv_rows", str(csv_hit)
        return -1, "raw_molecule_csv_unreadable", str(csv_hit)
    if parquet_hit is not None:
        n = count_parquet_rows(parquet_hit)
        if n >= 0:
            return n, "raw_molecule_parquet_rows", str(parquet_hit)
        return -1, "raw_molecule_parquet_no_pyarrow", str(parquet_hit)
    if budget.bail_reason:
        return -1, budget.bail_reason, ""
    return -1, "no_molecule_table", ""


# ── Per-dataset resolution ────────────────────────────────────────

def resolve_molecules(sample: str, bil_path: str, meyes_base: Path,
                       args) -> tuple[int | None, str, str]:
    """Return (n_genes_imaged, source, source_path). n is None on failure."""
    # 1. meyes processed single-molecule manifest.
    meyes_dir = meyes_base / sample
    hits, _ = find_molecule_manifests(
        meyes_dir, args.max_walk_depth,
        args.walk_timeout, args.max_walk_files)
    if hits:
        total, used = pick_manifest_total(hits, args.sum_regions)
        return total, "meyes_sm_manifest", used

    if not bil_path:
        return None, "no_bil_path", ""

    bil_root = Path(bil_path.rstrip("/"))

    # 2a. single-molecule manifest somewhere in the BIL path.
    hits, bail = find_molecule_manifests(
        bil_root, args.max_walk_depth,
        args.walk_timeout, args.max_walk_files)
    if hits:
        total, used = pick_manifest_total(hits, args.sum_regions)
        return total, "bil_sm_manifest", used

    # 2b. raw molecule table row count.
    n, label, src = molecules_from_raw_tables(
        bil_root, args.max_walk_depth,
        args.walk_timeout, args.max_walk_files)
    if n >= 0:
        return n, label, src
    return None, (label or bail or "no_molecule_source"), src


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

        out: dict = {
            "dataset": sample,
            "status": r.get("status", ""),
            "notes": r.get("notes", ""),
            "bil_path": bil_path,
            "publication_date": "",
            "year": "",
            "n_genes_imaged": "",
            "source": "",
            "source_path": "",
        }

        date = path_publication_date(bil_path) if bil_path else ""
        out["publication_date"] = date
        out["year"] = date[:4] if date else ""

        n, source, source_path = resolve_molecules(
            sample, bil_path, meyes_base, args)
        out["source"] = source
        out["source_path"] = source_path
        if n is not None:
            out["n_genes_imaged"] = n

        logger.info("[%s] date=%s n_genes_imaged=%s source=%s",
                    sample, out["publication_date"] or "—",
                    out["n_genes_imaged"] if out["n_genes_imaged"] != "" else "—",
                    out["source"] or "—")
        inventory.append(out)

    # ── Per-dataset inventory ─────────────────────────────────────
    inv_path = out_dir / "genes_inventory.csv"
    fieldnames = ["dataset", "status", "notes", "bil_path",
                  "publication_date", "year", "n_genes_imaged",
                  "source", "source_path"]
    with open(inv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(inventory)
    logger.info("Wrote %s (%d rows)", inv_path, len(inventory))

    # ── Aggregate by year + manifest ──────────────────────────────
    total_by_year: dict[str, int] = defaultdict(int)
    datasets_by_year: dict[str, int] = defaultdict(int)
    manifest_rows: list[dict] = []
    n_included = n_excluded = 0
    for r in inventory:
        year = r["year"]
        n_genes = r["n_genes_imaged"] if isinstance(r["n_genes_imaged"], int) else None

        reasons = []
        if not year:
            reasons.append("no_publication_date")
        if n_genes is None:
            reasons.append(f"no_gene_count:{r['source'] or 'unknown'}")
        included = not reasons

        if included:
            total_by_year[year] += n_genes  # type: ignore[arg-type]
            datasets_by_year[year] += 1
            n_included += 1
        else:
            n_excluded += 1

        manifest_rows.append({
            "dataset": r["dataset"],
            "included": "true" if included else "false",
            "year": year,
            "n_genes_imaged": n_genes if n_genes is not None else "",
            "reason": "; ".join(reasons),
        })

    yearly_path = out_dir / "genes_per_year.csv"
    with open(yearly_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["year", "n_datasets", "total_genes_imaged"])
        grand_total = 0
        grand_datasets = 0
        for y in sorted(total_by_year):
            writer.writerow([y, datasets_by_year[y], total_by_year[y]])
            grand_total += total_by_year[y]
            grand_datasets += datasets_by_year[y]
        writer.writerow(["TOTAL", grand_datasets, grand_total])
    logger.info("Wrote %s (%d years, %d genes imaged total)",
                yearly_path, len(total_by_year), grand_total)

    manifest_path = out_dir / "genes_per_year_manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["dataset", "included", "year",
                        "n_genes_imaged", "reason"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    logger.info("Wrote %s (%d included, %d excluded)",
                manifest_path, n_included, n_excluded)


if __name__ == "__main__":
    main()
