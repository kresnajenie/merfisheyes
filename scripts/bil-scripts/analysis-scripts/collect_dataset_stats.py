#!/usr/bin/env python3
"""
collect_dataset_stats.py

Collect per-dataset preprocessing statistics into a single growing CSV, one row
per dataset. Designed to run at the tail of stage 1 (after combine + map_my_cell)
so it can read both the combined matrices and the MapMyCells output.

For each dataset folder under ${MEYES_BASE} (default /bil/data/meyes), reads:
  <meyes_base>/<dataset>/combined_output/cell_metadata.csv   (samples, cells)
  <meyes_base>/<dataset>/combined_output/cell_by_gene.csv    (genes, counts)
  <meyes_base>/<dataset>/combined_output/artifact_mask_p<P>.csv  (artifact flags)
  <meyes_base>/<dataset>/mmc_output/mapping_output.csv        (cell types)

Metrics collected:
  n_samples            distinct sections combined (unique _sample_id)
  n_cells              total cells
  n_genes              panel width (cell_by_gene columns minus the cell-id column)
  n_genes_detected     genes with >=1 count anywhere in the dataset
  total_counts         summed transcript counts across all cells
  mean/median_counts_per_cell   per-cell total-count distribution
  n_artifact / pct_artifact     cells flagged by the p<P> artifact mask
  n_clusters           distinct leaf clusters assigned (cluster_name)
  n_classes            distinct broad classes assigned (class_name/supercluster_name)
  mapping_method       'correlation' or 'hierarchical' (from CSV header / default)

Behavior (mirrors collect_correlation_sections.py):
  - If --out exists, it is loaded and merged with this run.
  - Any dataset scraped this run REPLACES its existing row.
  - Datasets not mentioned this run are left untouched.

Usage:
  # Single dataset (how the stage-1 launcher calls it, after MMC)
  python collect_dataset_stats.py --datasets ace-irk-sag --species human \
      --out /bil/data/meyes/_dataset_stats/dataset_stats.csv

  # Several at once
  python collect_dataset_stats.py --datasets ace-den-fix ace-dip-use \
      --out /bil/data/meyes/_dataset_stats/dataset_stats.csv
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import fcntl
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


COLUMNS = [
    "dataset", "species", "n_samples", "n_cells", "n_genes", "n_genes_detected",
    "total_counts", "mean_counts_per_cell", "median_counts_per_cell",
    "mask_percentile", "n_artifact", "pct_artifact",
    "n_clusters", "n_classes", "mapping_method", "scraped_at",
]

META_ID_CANDIDATES = ("EntityID", "id", "cell_id")
CLUSTER_COLS = ("cluster_name", "cluster_label")
CLASS_COLS = ("class_name", "supercluster_name", "subclass_name")
CHUNK = 100_000


def _meta_paths(meyes_base: Path, sample: str):
    combined = meyes_base / sample / "combined_output"
    mmc = meyes_base / sample / "mmc_output" / "mapping_output.csv"
    return combined, mmc


def _samples_and_cells(meta_path: Path):
    """Return (n_samples, n_cells) from cell_metadata.csv."""
    header = pd.read_csv(meta_path, nrows=0).columns.tolist()
    if "_sample_id" in header:
        sid = pd.read_csv(meta_path, usecols=["_sample_id"])
        return int(sid["_sample_id"].nunique()), int(len(sid))
    # No section column (e.g. a single-section combine) — count rows only.
    n_cells = sum(len(c) for c in pd.read_csv(meta_path, usecols=[0], chunksize=CHUNK))
    return pd.NA, int(n_cells)


def _gene_stats(cbg_path: Path):
    """Single chunked pass over cell_by_gene.csv.

    Returns (n_genes, n_genes_detected, total_counts, mean_per_cell,
    median_per_cell).
    """
    gene_cols = None
    detected = None          # bool per gene: has any nonzero count
    per_cell_sums = []       # one total per cell

    for chunk in pd.read_csv(cbg_path, chunksize=CHUNK):
        cols = [c for c in chunk.columns if c != "cell"]
        if gene_cols is None:
            gene_cols = cols
            detected = np.zeros(len(gene_cols), dtype=bool)
        mat = chunk[gene_cols].to_numpy()
        per_cell_sums.append(mat.sum(axis=1))
        detected |= (mat > 0).any(axis=0)

    if gene_cols is None:
        raise ValueError(f"{cbg_path} has no data rows")

    per_cell = np.concatenate(per_cell_sums)
    return (
        len(gene_cols),
        int(detected.sum()),
        float(per_cell.sum()),
        float(per_cell.mean()),
        float(np.median(per_cell)),
    )


def _artifact_stats(combined: Path, percentile: int, n_cells: int):
    """Return (n_artifact, pct_artifact) from artifact_mask_p<P>.csv, or (NA, NA)."""
    mask_path = combined / f"artifact_mask_p{percentile}.csv"
    if not mask_path.exists():
        return pd.NA, pd.NA
    flag_col = "is_artifact"
    cols = pd.read_csv(mask_path, nrows=0).columns.tolist()
    if flag_col not in cols:
        # Fall back to the last column if the schema ever changes.
        flag_col = cols[-1]
    flags = pd.read_csv(mask_path, usecols=[flag_col])[flag_col]
    is_art = flags.astype(str).str.strip().str.lower().isin(("true", "1"))
    n_art = int(is_art.sum())
    pct = round(n_art / n_cells * 100, 3) if n_cells else pd.NA
    return n_art, pct


def _celltype_stats(mmc_csv: Path):
    """Return (n_clusters, n_classes, mapping_method) from mapping_output.csv."""
    if not mmc_csv.exists():
        return pd.NA, pd.NA, pd.NA

    # mapping_method is recorded as a leading "# mapping method: ..." comment by
    # the correlation runner; hierarchical writes no comment.
    method = "hierarchical"
    with open(mmc_csv) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "mapping method:" in line:
                method = line.split("mapping method:", 1)[1].strip().split()[0]
                break

    df = pd.read_csv(mmc_csv, comment="#")

    def _nunique(candidates):
        for col in candidates:
            if col in df.columns and not df[col].isna().all():
                return int(df[col].nunique())
        return pd.NA

    return _nunique(CLUSTER_COLS), _nunique(CLASS_COLS), method


def scrape_dataset(name: str, meyes_base: Path, percentile: int, species: str) -> dict:
    combined, mmc_csv = _meta_paths(meyes_base, name)
    meta_path = combined / "cell_metadata.csv"
    cbg_path = combined / "cell_by_gene.csv"
    if not meta_path.exists() or not cbg_path.exists():
        raise FileNotFoundError(f"combined_output CSVs not found in {combined}")

    n_samples, n_cells = _samples_and_cells(meta_path)
    n_genes, n_detected, total, mean_pc, median_pc = _gene_stats(cbg_path)
    n_art, pct_art = _artifact_stats(combined, percentile, n_cells)
    n_clusters, n_classes, method = _celltype_stats(mmc_csv)

    return {
        "dataset": name,
        "species": species or pd.NA,
        "n_samples": n_samples,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_genes_detected": n_detected,
        "total_counts": total,
        "mean_counts_per_cell": round(mean_pc, 3),
        "median_counts_per_cell": round(median_pc, 3),
        "mask_percentile": percentile,
        "n_artifact": n_art,
        "pct_artifact": pct_art,
        "n_clusters": n_clusters,
        "n_classes": n_classes,
        "mapping_method": method,
        "scraped_at": dt.datetime.now().isoformat(timespec="seconds"),
    }


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--datasets", nargs="+", required=True,
                   help="Dataset folder name(s) to scrape.")
    p.add_argument("--out", required=True,
                   help="Output accumulator CSV. Created or merged into.")
    p.add_argument("--species", default="",
                   help="Species label applied to the scraped dataset(s).")
    p.add_argument("--mask-percentile", type=int, default=25,
                   help="Artifact-mask percentile to read (default: 25).")
    p.add_argument("--meyes-base",
                   default=os.environ.get("MEYES_BASE", "/bil/data/meyes"),
                   help="Root containing one folder per dataset (default: %(default)s).")
    return p.parse_args()


@contextlib.contextmanager
def _accumulator_lock(out_path: Path):
    """Exclusive lock so concurrent per-sample jobs can't clobber the CSV.

    The heavy scraping happens before this is taken, so the lock is held only
    for the brief load → merge → write critical section.
    """
    lock_path = out_path.with_name(out_path.name + ".lock")
    with open(lock_path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def main():
    args = parse_args()
    meyes_base = Path(args.meyes_base)
    out_path = Path(args.out).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Heavy work first, without the lock held.
    new_rows = []
    for name in args.datasets:
        try:
            row = scrape_dataset(name, meyes_base, args.mask_percentile, args.species)
        except (FileNotFoundError, ValueError) as e:
            print(f"[skip ] {name}: {e}", file=sys.stderr)
            continue
        new_rows.append(row)
        print(f"[ok   ] {name}: samples={row['n_samples']} cells={row['n_cells']:,} "
              f"genes={row['n_genes']} ({row['n_genes_detected']} detected) "
              f"clusters={row['n_clusters']} classes={row['n_classes']}")

    if not new_rows:
        sys.exit("Nothing to write — no datasets produced usable rows.")

    new_df = pd.DataFrame(new_rows, columns=COLUMNS)
    scraped = set(new_df["dataset"])

    # Critical section: re-read the latest file under lock, merge, write.
    with _accumulator_lock(out_path):
        if out_path.exists():
            existing = pd.read_csv(out_path)
            for col in COLUMNS:  # tolerate older files missing new columns
                if col not in existing.columns:
                    existing[col] = pd.NA
            existing = existing[COLUMNS]
            print(f"[load] {out_path}: {len(existing)} dataset row(s)")
        else:
            existing = pd.DataFrame(columns=COLUMNS)
            print(f"[load] {out_path}: not present — creating new file")

        kept = existing[~existing["dataset"].isin(scraped)]
        combined = pd.concat([kept, new_df], ignore_index=True)
        combined = combined.sort_values("dataset").reset_index(drop=True)
        combined.to_csv(out_path, index=False)
        replaced = sorted(set(existing["dataset"]) & scraped)

    print(f"\nWrote {out_path}: {len(combined)} dataset row(s)")
    if replaced:
        print(f"       Replaced rows for: {', '.join(replaced)}")


if __name__ == "__main__":
    main()
