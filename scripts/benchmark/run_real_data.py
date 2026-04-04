#!/usr/bin/env python3
"""Run browser benchmarks on real datasets.

Expects a folder organized as:
    test-data/
    ├── single-cell/
    │   ├── h5ad/        *.h5ad files
    │   ├── xenium/      folders with cells.csv + features.tsv + cell_feature_matrix/
    │   └── merscope/    folders with cell_metadata.csv + cell_by_gene.csv
    └── single-molecule/
        ├── xenium/      .parquet or .csv files (feature_name columns)
        └── merscope/    .parquet or .csv files (gene columns)

Output: CSV with the same schema as run_browser_benchmarks.py

Usage:
    python run_real_data.py --data-dir /path/to/test-data
    python run_real_data.py --data-dir /path/to/test-data --output real_results.csv
    python run_real_data.py --data-dir /path/to/test-data --types single_cell
"""

import argparse
import asyncio
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    BenchmarkResult, DEV_SERVER_URL, BENCHMARK_TIMEOUT_S, DEFAULT_RESULTS_DIR,
)
from run_browser_benchmarks import (
    BrowserBenchmarkRunner, append_result_csv, get_file_size_mb, load_existing_results,
)


# ---------------------------------------------------------------------------
# Detect real datasets
# ---------------------------------------------------------------------------
def _count_cells_h5ad(path: str) -> tuple[int, int]:
    """Quick read of n_obs and n_vars from an h5ad file."""
    try:
        import h5py
        with h5py.File(path, "r") as f:
            if "X" in f:
                shape = f["X"].shape if hasattr(f["X"], "shape") else None
                if shape:
                    return shape[0], shape[1]
                # Sparse CSR
                if "data" in f["X"] and "indptr" in f["X"]:
                    n_cells = len(f["X"]["indptr"]) - 1
                    n_genes = f["X"].attrs.get("shape", [0, 0])[1] if "shape" in f["X"].attrs else 0
                    return n_cells, n_genes
            # Try obs/_index length
            if "obs" in f and "_index" in f["obs"]:
                n_cells = len(f["obs"]["_index"])
            elif "obs" in f:
                n_cells = f["obs"].attrs.get("_index_length", 0)
            else:
                n_cells = 0
            if "var" in f and "_index" in f["var"]:
                n_genes = len(f["var"]["_index"])
            else:
                n_genes = 0
            return n_cells, n_genes
    except Exception:
        return 0, 0


def _count_cells_csv(path: str) -> int:
    """Count lines in a CSV (minus header)."""
    count = 0
    with open(path) as f:
        for _ in f:
            count += 1
    return max(0, count - 1)


def _count_molecules_parquet(path: str) -> tuple[int, int]:
    """Read row count and unique genes from a parquet file."""
    try:
        import pyarrow.parquet as pq
        meta = pq.read_metadata(path)
        n_rows = meta.num_rows
        # Read just the gene column to count unique genes
        table = pq.read_table(path, columns=["feature_name"])
        n_genes = len(table.column("feature_name").unique())
        return n_rows, n_genes
    except Exception:
        try:
            table = pq.read_table(path, columns=["gene"])
            n_genes = len(table.column("gene").unique())
            return meta.num_rows, n_genes
        except Exception:
            return 0, 0


def scan_real_single_cell(data_dir: str) -> list[dict]:
    """Discover real single-cell datasets from directory structure."""
    sc_dir = os.path.join(data_dir, "single-cell")
    if not os.path.exists(sc_dir):
        sc_dir = os.path.join(data_dir, "single_cell")
    if not os.path.exists(sc_dir):
        return []

    datasets = []

    # H5AD files
    h5ad_dir = os.path.join(sc_dir, "h5ad")
    if os.path.isdir(h5ad_dir):
        for f in sorted(os.listdir(h5ad_dir)):
            if f.endswith(".h5ad"):
                path = os.path.join(h5ad_dir, f)
                n_cells, n_genes = _count_cells_h5ad(path)
                datasets.append({
                    "n_cells": n_cells, "n_genes": n_genes,
                    "filetype": "h5ad", "path": path,
                    "file_size_mb": round(get_file_size_mb(path), 2),
                })

    # Xenium folders
    xenium_dir = os.path.join(sc_dir, "xenium")
    if os.path.isdir(xenium_dir):
        for folder in sorted(os.listdir(xenium_dir)):
            folder_path = os.path.join(xenium_dir, folder)
            if os.path.isdir(folder_path) and os.path.exists(os.path.join(folder_path, "cells.csv")):
                n_cells = _count_cells_csv(os.path.join(folder_path, "cells.csv"))
                # Count genes from features.tsv if present
                n_genes = 0
                for feat_name in ["features.tsv", "features.csv"]:
                    feat_path = os.path.join(folder_path, feat_name)
                    if os.path.exists(feat_path):
                        n_genes = max(0, _count_cells_csv(feat_path))
                        break
                datasets.append({
                    "n_cells": n_cells, "n_genes": n_genes,
                    "filetype": "xenium", "path": folder_path,
                    "file_size_mb": round(get_file_size_mb(folder_path), 2),
                })

    # MERSCOPE folders
    merscope_dir = os.path.join(sc_dir, "merscope")
    if os.path.isdir(merscope_dir):
        for folder in sorted(os.listdir(merscope_dir)):
            folder_path = os.path.join(merscope_dir, folder)
            meta_path = os.path.join(folder_path, "cell_metadata.csv")
            if os.path.isdir(folder_path) and os.path.exists(meta_path):
                n_cells = _count_cells_csv(meta_path)
                # Count genes from cell_by_gene.csv header
                n_genes = 0
                cbg_path = os.path.join(folder_path, "cell_by_gene.csv")
                if os.path.exists(cbg_path):
                    with open(cbg_path) as f:
                        header = f.readline().strip().split(",")
                        n_genes = len(header) - 1  # first col is cell id
                datasets.append({
                    "n_cells": n_cells, "n_genes": n_genes,
                    "filetype": "merscope", "path": folder_path,
                    "file_size_mb": round(get_file_size_mb(folder_path), 2),
                })

    return datasets


def scan_real_single_molecule(data_dir: str) -> list[dict]:
    """Discover real single-molecule datasets from directory structure."""
    sm_dir = os.path.join(data_dir, "single-molecule")
    if not os.path.exists(sm_dir):
        sm_dir = os.path.join(data_dir, "single_molecule")
    if not os.path.exists(sm_dir):
        return []

    datasets = []
    for subfolder in ["xenium", "merscope"]:
        sub_dir = os.path.join(sm_dir, subfolder)
        if not os.path.isdir(sub_dir):
            continue
        for f in sorted(os.listdir(sub_dir)):
            path = os.path.join(sub_dir, f)
            if f.endswith(".parquet"):
                n_mol, n_genes = _count_molecules_parquet(path)
                datasets.append({
                    "n_molecules": n_mol, "n_genes": n_genes,
                    "filetype": "parquet", "path": path,
                    "file_size_mb": round(get_file_size_mb(path), 2),
                    "source_format": subfolder,
                })
            elif f.endswith(".csv"):
                n_mol = _count_cells_csv(path)
                datasets.append({
                    "n_molecules": n_mol, "n_genes": 0,  # unknown without reading full file
                    "filetype": "csv", "path": path,
                    "file_size_mb": round(get_file_size_mb(path), 2),
                    "source_format": subfolder,
                })

    return datasets


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
async def run_real_benchmarks(args):
    runner = BrowserBenchmarkRunner(
        base_url=args.url,
        timeout_s=args.timeout,
        headless=args.headless,
    )

    if not runner.is_server_up():
        print(f"ERROR: Dev server not running at {args.url}")
        sys.exit(1)

    print(f"Real data benchmark")
    print(f"  Data dir: {args.data_dir}")
    print(f"  Output:   {args.output}")
    print()

    tested = load_existing_results(args.output) if args.resume else set()
    completed = 0
    errors = 0

    if "single_cell" in args.types:
        sc_datasets = scan_real_single_cell(args.data_dir)
        print(f"Found {len(sc_datasets)} single cell datasets")
        for i, ds in enumerate(sc_datasets):
            if ds["path"] in tested:
                print(f"  [{i+1}] SKIP: {ds['path']}")
                continue
            print(f"  [{i+1}/{len(sc_datasets)}] {ds['filetype']} — "
                  f"{ds['n_cells']:,} cells × {ds['n_genes']:,} genes "
                  f"({ds['file_size_mb']:.1f} MB)")
            result = await runner.run_single_cell(
                path=ds["path"], filetype=ds["filetype"],
                n_cells=ds["n_cells"], n_genes=ds["n_genes"],
            )
            result.notes = "real_data"
            append_result_csv(args.output, result)
            if result.crashed:
                errors += 1
                print(f"    CRASHED: {result.error}")
            else:
                completed += 1
                print(f"    load={result.load_time_s}s  mem={result.memory_peak_mb}MB")

    if "single_molecule" in args.types:
        sm_datasets = scan_real_single_molecule(args.data_dir)
        print(f"\nFound {len(sm_datasets)} single molecule datasets")
        for i, ds in enumerate(sm_datasets):
            if ds["path"] in tested:
                print(f"  [{i+1}] SKIP: {ds['path']}")
                continue
            print(f"  [{i+1}/{len(sm_datasets)}] {ds['filetype']} — "
                  f"{ds.get('n_molecules', '?'):,} molecules "
                  f"({ds['file_size_mb']:.1f} MB)")
            result = await runner.run_single_molecule(
                path=ds["path"], filetype=ds["filetype"],
                n_molecules=ds.get("n_molecules", 0),
                n_genes=ds.get("n_genes", 0),
            )
            result.notes = f"real_data:{ds.get('source_format', '')}"
            append_result_csv(args.output, result)
            if result.crashed:
                errors += 1
                print(f"    CRASHED: {result.error}")
            else:
                completed += 1
                print(f"    load={result.load_time_s}s  mem={result.memory_peak_mb}MB")

    print(f"\nDone. {completed} passed, {errors} crashed.")
    print(f"Results: {args.output}")


def main():
    parser = argparse.ArgumentParser(description="Run browser benchmarks on real datasets")
    parser.add_argument("--data-dir", required=True,
                        help="Root folder with single-cell/ and single-molecule/ subdirectories")
    parser.add_argument("--output", default=os.path.join(DEFAULT_RESULTS_DIR, "real_data_benchmarks.csv"),
                        help="Output CSV path")
    parser.add_argument("--url", default=DEV_SERVER_URL)
    parser.add_argument("--timeout", type=int, default=BENCHMARK_TIMEOUT_S)
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--types", nargs="+", default=["single_cell", "single_molecule"],
                        choices=["single_cell", "single_molecule"])
    args = parser.parse_args()
    asyncio.run(run_real_benchmarks(args))


if __name__ == "__main__":
    DEFAULT_RESULTS_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "benchmark_results",
    )
    main()
