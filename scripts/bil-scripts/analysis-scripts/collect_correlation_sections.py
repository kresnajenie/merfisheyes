#!/usr/bin/env python3
"""
collect_correlation_sections.py

Scrape per-section mean correlation coefficients from mapping_output.csv files
produced by `map_my_cell.py --method correlation`, and accumulate them in a
single TSV that can be grown over time.

For each dataset folder under ${MEYES_BASE} (default /bil/data/meyes), reads:
  <meyes_base>/<dataset>/mmc_output_corr/mapping_output.csv

Section ID is the prefix of cell_id before the first underscore — build_h5ad
writes cell_id as "{_sample_id}_{original_cell_id}", so the prefix is the
slice/sample id.

Output TSV columns:
  dataset, section_id, mean_r, median_r, n_cells, source_csv, scraped_at

Behavior:
  - If --out file exists, rows are loaded and merged with the new scrape.
  - Any dataset being scraped this run REPLACES its existing rows in the TSV
    (so re-running with the same name picks up fresh mapping_output.csv data).
  - Datasets not mentioned this run are left untouched.

Usage:
  # First run
  python collect_correlation_sections.py \
      --datasets ace-den-fix ace-dip-use \
      --out ~/correlation_sections.tsv

  # Later, add more — existing rows for ace-den-fix/ace-dip-use are preserved
  python collect_correlation_sections.py \
      --datasets new-dataset-1 new-dataset-2 \
      --out ~/correlation_sections.tsv
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
from pathlib import Path

import pandas as pd


REQUIRED_COLS = ["dataset", "section_id", "mean_r", "median_r",
                 "n_cells", "source_csv", "scraped_at"]


def mapping_csv(meyes_base: Path, sample: str) -> Path:
    return meyes_base / sample / "mmc_output_corr" / "mapping_output.csv"


def scrape_dataset(name: str, csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, comment="#")
    for col in ("cell_id", "correlation_coefficient"):
        if col not in df.columns:
            raise ValueError(f"{csv_path} missing required column '{col}'")

    df["correlation_coefficient"] = pd.to_numeric(
        df["correlation_coefficient"], errors="coerce"
    )
    df = df.dropna(subset=["correlation_coefficient"])
    if df.empty:
        raise ValueError(f"{csv_path} has no numeric correlation_coefficient values")

    df["section_id"] = df["cell_id"].astype(str).str.split("_", n=1).str[0]
    grouped = (
        df.groupby("section_id")["correlation_coefficient"]
        .agg(mean_r="mean", median_r="median", n_cells="size")
        .reset_index()
    )
    grouped["dataset"] = name
    grouped["source_csv"] = str(csv_path)
    grouped["scraped_at"] = dt.datetime.now().isoformat(timespec="seconds")
    return grouped[REQUIRED_COLS]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--datasets", nargs="+", required=True,
                   help="Dataset folder names to scrape (one or more).")
    p.add_argument("--out", required=True,
                   help="Output TSV path. Will be created or appended to.")
    p.add_argument("--meyes-base",
                   default=os.environ.get("MEYES_BASE", "/bil/data/meyes"),
                   help="Root containing one folder per dataset (default: %(default)s).")
    return p.parse_args()


def main():
    args = parse_args()
    meyes_base = Path(args.meyes_base)
    out_path = Path(args.out).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Load existing accumulator if present
    if out_path.exists():
        existing = pd.read_csv(out_path, sep="\t")
        missing = [c for c in REQUIRED_COLS if c not in existing.columns]
        if missing:
            sys.exit(f"Existing {out_path} is missing expected columns: {missing}")
        print(f"[load] {out_path}: "
              f"{existing['dataset'].nunique()} datasets, {len(existing)} rows")
    else:
        existing = pd.DataFrame(columns=REQUIRED_COLS)
        print(f"[load] {out_path}: not present — creating new file")

    new_rows = []
    missing_csvs = []
    for name in args.datasets:
        csv = mapping_csv(meyes_base, name)
        if not csv.exists():
            missing_csvs.append((name, csv))
            continue
        try:
            scraped = scrape_dataset(name, csv)
        except ValueError as e:
            print(f"[skip ] {name}: {e}", file=sys.stderr)
            continue
        new_rows.append(scraped)
        print(f"[ok   ] {name}: sections={len(scraped)} "
              f"cells={int(scraped['n_cells'].sum())} "
              f"mean_r={scraped['mean_r'].mean():.4f}")

    if missing_csvs:
        print("[warn] missing mapping_output.csv:", file=sys.stderr)
        for n, c in missing_csvs:
            print(f"         {n}: {c}", file=sys.stderr)

    if not new_rows:
        sys.exit("Nothing to write — no datasets produced usable rows.")

    new_df = pd.concat(new_rows, ignore_index=True)
    scraped_names = set(new_df["dataset"].unique())

    # Drop rows for datasets we just re-scraped, then append fresh rows
    kept = existing[~existing["dataset"].isin(scraped_names)]
    combined = pd.concat([kept, new_df], ignore_index=True)
    combined = combined.sort_values(["dataset", "section_id"]).reset_index(drop=True)

    combined.to_csv(out_path, sep="\t", index=False)
    print(f"\nWrote {out_path}: "
          f"{combined['dataset'].nunique()} datasets, {len(combined)} rows")
    if not kept.empty:
        replaced = sorted(set(existing['dataset']) & scraped_names)
        if replaced:
            print(f"       Replaced rows for: {', '.join(replaced)}")


if __name__ == "__main__":
    main()
