#!/usr/bin/env python3
"""
count_genes_per_year.py

Like count_cells_per_year.py, but counts genes instead of cells.

For each dataset in the all-datasets spreadsheet:
  - publication_date = mtime of the BIL data path.
  - gene panel = column names from cell_by_gene (the post-pipeline file in
    /bil/data/meyes/<dataset>/combined_output/, with a fallback walk of the
    BIL data path for cell_by_gene.csv / cell_by_gene.h5ad).
  - n_genes = len(panel).

Outputs (in --out-dir):
  genes_inventory.csv             one row per dataset with date, year, n_genes, source, genes (semicolon-joined)
  genes_per_year.csv              year, n_datasets, sum_genes, unique_genes
  genes_per_year_manifest.csv     per-dataset included/excluded with reasons

sum_genes  = sum of per-dataset panel sizes (shared genes double-counted)
unique_genes = size of the union of gene panels across all datasets in the year

Usage:
  python count_genes_per_year.py                          # defaults
  python count_genes_per_year.py --datasets-csv path.csv \
      --out-dir /bil/data/meyes/_genes_inventory
"""

from __future__ import annotations

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

logger = logging.getLogger("count_genes")

DEFAULT_DATASETS = "scripts/bil-scripts/all_datasets.csv"
DEFAULT_MEYES_BASE = "/bil/data/meyes"
DEFAULT_OUT = "/bil/data/meyes/_genes_inventory"

MAX_WALK_DEPTH = 4
DEFAULT_WALK_TIMEOUT = 60.0
DEFAULT_MAX_WALK_FILES = 20000

# How many gene names to embed in the inventory CSV (full list kept for the
# year-level union; this is just a peek so the inventory stays readable).
INVENTORY_GENE_PREVIEW = 25

# Raise the csv field-size limit — cell_by_gene.csv headers can be huge.
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
                   help=f"How many subdir levels to scan under each BIL "
                        f"path when looking for cell_by_gene "
                        f"(default: {MAX_WALK_DEPTH}).")
    p.add_argument("--walk-timeout", type=float, default=DEFAULT_WALK_TIMEOUT,
                   help=f"Per-dataset wall-clock cap on the BIL fallback "
                        f"walk, in seconds (default: {DEFAULT_WALK_TIMEOUT}).")
    p.add_argument("--max-walk-files", type=int, default=DEFAULT_MAX_WALK_FILES,
                   help=f"Per-dataset cap on files inspected during the BIL "
                        f"fallback walk (default: {DEFAULT_MAX_WALK_FILES}).")
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


def read_cell_by_gene_header(path: Path) -> list[str]:
    """Read the first row of cell_by_gene.csv (.gz ok) and return the gene
    names — every column except the first (cell ID). Returns [] on failure."""
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", errors="replace") as fh:
            reader = csv.reader(fh)
            header = next(reader, None)
    except OSError as e:
        logger.warning("Failed to read %s: %s", path, e)
        return []
    if not header:
        return []
    return [c.strip() for c in header[1:] if c.strip()]


def read_h5ad_var_names(path: Path) -> list[str]:
    """Return var (gene) names from an h5ad file, or [] on failure."""
    try:
        import h5py
    except ImportError:
        logger.warning("h5py not available — skipping h5ad at %s", path)
        return []
    try:
        with h5py.File(str(path), "r") as f:
            var = f.get("var")
            if var is None:
                return []
            for key in ("_index", "index", "feature_name", "gene_name",
                        "gene_symbol", "gene_ids"):
                if key in var:
                    ds = var[key]
                    if hasattr(ds, "shape") and ds.shape:
                        return [
                            v.decode("utf-8") if isinstance(v, bytes) else str(v)
                            for v in ds[:]
                        ]
            # Last resort: first string-typed dataset under var.
            for k in var.keys():
                ds = var[k]
                if hasattr(ds, "shape") and ds.shape:
                    try:
                        return [
                            v.decode("utf-8") if isinstance(v, bytes) else str(v)
                            for v in ds[:]
                        ]
                    except Exception:
                        continue
            return []
    except Exception as e:
        logger.warning("Failed to read h5ad %s: %s", path, e)
        return []


def genes_from_meyes(meyes_base: Path, sample: str) -> tuple[list[str], str, str]:
    """Try the canonical post-pipeline cell_by_gene first."""
    combined = meyes_base / sample / "combined_output"
    for name, label in (
        ("cell_by_gene.csv", "meyes_cell_by_gene_csv"),
        ("cell_by_gene.csv.gz", "meyes_cell_by_gene_csv"),
        ("cell_by_gene.h5ad", "meyes_cell_by_gene_h5ad"),
    ):
        candidate = combined / name
        try:
            exists = candidate.exists()
        except PermissionError:
            return [], "permission_denied", ""
        if not exists:
            continue
        genes = (read_h5ad_var_names(candidate)
                 if candidate.suffix == ".h5ad"
                 else read_cell_by_gene_header(candidate))
        if genes:
            return genes, label, str(candidate)
        return [], f"{label}_unreadable", str(candidate)
    return [], "no_meyes_cell_by_gene", ""


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


# Ordered list of (matcher, reader, label).
GENE_FILE_PATTERNS: list[tuple] = [
    (lambda n: n == "cell_by_gene.csv" or n == "cell_by_gene.csv.gz",
        read_cell_by_gene_header, "bil_cell_by_gene_csv"),
    (lambda n: n == "cell_by_gene.h5ad",
        read_h5ad_var_names, "bil_cell_by_gene_h5ad"),
    (lambda n: n.endswith(".h5ad"),
        read_h5ad_var_names, "bil_h5ad_other"),
]


def genes_from_bil_path(root: Path, max_depth: int,
                        walk_timeout: float,
                        max_walk_files: int) -> tuple[list[str], str, str]:
    """Walk root looking for the first recognized gene-bearing file."""
    try:
        if not root.exists():
            return [], "missing_path", ""
    except PermissionError:
        return [], "permission_denied", ""
    except OSError as e:
        return [], f"stat_error:{e.errno}", ""
    budget = WalkBudget(max_files=max_walk_files,
                         deadline=time.monotonic() + walk_timeout)
    found_by_pattern: dict[int, Path] = {}
    for f in walk_files(root, max_depth, budget):
        for idx, (matcher, _, _) in enumerate(GENE_FILE_PATTERNS):
            if matcher(f.name):
                if idx not in found_by_pattern:
                    found_by_pattern[idx] = f
                break
        if 0 in found_by_pattern:
            break
    for idx in range(len(GENE_FILE_PATTERNS)):
        if idx not in found_by_pattern:
            continue
        path = found_by_pattern[idx]
        _, reader, label = GENE_FILE_PATTERNS[idx]
        genes = reader(path)
        if genes:
            try:
                rel = str(path.relative_to(root))
            except ValueError:
                rel = str(path)
            return genes, label, rel
    if budget.bail_reason:
        return [], budget.bail_reason, ""
    return [], "no_recognized_files", ""


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
    # Keep the full per-dataset gene list around so we can compute yearly
    # unions without re-reading the inventory CSV.
    panels: dict[str, list[str]] = {}

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
            "n_genes": "",
            "gene_source": "",
            "gene_source_path": "",
            "in_meyes": "",
            "genes_preview": "",
        }

        date = path_publication_date(bil_path) if bil_path else ""
        out["publication_date"] = date
        out["year"] = date[:4] if date else ""

        # 1. Try meyes combined_output first.
        genes, source, source_path = genes_from_meyes(meyes_base, sample)
        out["in_meyes"] = "true" if genes else "false"

        # 2. Fall back to the BIL path itself.
        if not genes and bil_path:
            genes, source, source_path = genes_from_bil_path(
                Path(bil_path.rstrip("/")), args.max_walk_depth,
                walk_timeout=args.walk_timeout,
                max_walk_files=args.max_walk_files,
            )

        if genes:
            out["n_genes"] = len(genes)
            out["gene_source"] = source
            out["gene_source_path"] = source_path
            out["genes_preview"] = ";".join(genes[:INVENTORY_GENE_PREVIEW])
            panels[sample] = genes
        else:
            out["gene_source"] = source

        logger.info("[%s] date=%s n_genes=%s source=%s",
                    sample, out["publication_date"] or "—",
                    out["n_genes"] if out["n_genes"] != "" else "—",
                    out["gene_source"] or "—")
        inventory.append(out)

    # ── Write per-dataset inventory ───────────────────────────────
    inv_path = out_dir / "genes_inventory.csv"
    fieldnames = ["dataset", "status", "notes", "bil_path",
                  "publication_date", "year", "n_genes",
                  "gene_source", "gene_source_path", "in_meyes",
                  "genes_preview"]
    with open(inv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(inventory)
    logger.info("Wrote %s (%d rows)", inv_path, len(inventory))

    # ── Aggregate by year + manifest ──────────────────────────────
    sum_by_year: dict[str, int] = defaultdict(int)
    union_by_year: dict[str, set[str]] = defaultdict(set)
    datasets_by_year: dict[str, int] = defaultdict(int)
    manifest_rows: list[dict] = []
    n_included = 0
    n_excluded = 0
    for r in inventory:
        year = r["year"]
        n_genes = r["n_genes"] if isinstance(r["n_genes"], int) else None

        reasons = []
        if not year:
            reasons.append("no_publication_date")
        if n_genes is None:
            reasons.append(f"no_gene_count:{r['gene_source'] or 'unknown'}")
        included = not reasons

        if included:
            sum_by_year[year] += n_genes  # type: ignore[arg-type]
            datasets_by_year[year] += 1
            union_by_year[year].update(panels.get(r["dataset"], []))
            n_included += 1
        else:
            n_excluded += 1

        manifest_rows.append({
            "dataset": r["dataset"],
            "included": "true" if included else "false",
            "year": year,
            "n_genes": n_genes if n_genes is not None else "",
            "reason": "; ".join(reasons),
        })

    yearly_path = out_dir / "genes_per_year.csv"
    with open(yearly_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["year", "n_datasets", "sum_genes", "unique_genes"])
        for y in sorted(sum_by_year):
            writer.writerow([y, datasets_by_year[y],
                             sum_by_year[y], len(union_by_year[y])])
    logger.info("Wrote %s (%d years)", yearly_path, len(sum_by_year))

    manifest_path = out_dir / "genes_per_year_manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["dataset", "included", "year", "n_genes", "reason"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    logger.info("Wrote %s (%d included, %d excluded)",
                manifest_path, n_included, n_excluded)


if __name__ == "__main__":
    main()
