#!/usr/bin/env python3
"""Generate synthetic single-cell spatial transcriptomics datasets.

Produces h5ad, Xenium, and MERSCOPE format files at configurable scales.
288 combinations by default (12 cell counts × 8 gene counts × 3 filetypes).

Usage:
    python generate_single_cell.py
    python generate_single_cell.py --cells 1000 10000 --genes 100 500 --filetypes h5ad
    python generate_single_cell.py --dry-run
    python generate_single_cell.py --max-file-size-gb 5
"""

import argparse
import gzip
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Add parent to path so config is importable when run directly
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    SINGLE_CELL_SIZES, SINGLE_CELL_GENES, SINGLE_CELL_FILETYPES,
    SINGLE_CELL_OBS_COLUMNS, CELL_TYPES, LEIDEN_CLUSTERS,
    MAX_FILE_SIZE_GB, get_sparsity,
    estimate_h5ad_size_mb, estimate_xenium_size_mb, estimate_merscope_size_mb,
)

CHUNK_SIZE = 100_000  # rows per chunk for streaming writes


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def make_gene_names(n_genes: int) -> list[str]:
    return [f"gene_{i}" for i in range(n_genes)]


def make_obs_columns(n_cells: int, n_cols: int, rng: np.random.Generator) -> dict:
    """Generate categorical observation columns."""
    cols: dict[str, np.ndarray] = {}
    # Fixed columns
    cols["leiden"] = rng.choice(LEIDEN_CLUSTERS, size=n_cells)
    cols["celltype"] = rng.choice(CELL_TYPES, size=n_cells)
    # Additional categorical columns
    for i in range(2, n_cols):
        n_cats = rng.integers(5, 30)
        cats = [f"cat{i}_{j}" for j in range(n_cats)]
        cols[f"obs_col_{i}"] = rng.choice(cats, size=n_cells)
    return cols


def make_spatial(n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Random 2-D spatial coordinates in [-1000, 1000]."""
    return rng.uniform(-1000, 1000, size=(n_cells, 2)).astype(np.float32)


def make_umap(n_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Random 2-D UMAP embedding in [-10, 10]."""
    return rng.uniform(-10, 10, size=(n_cells, 2)).astype(np.float32)


def estimate_size(filetype: str, n_cells: int, n_genes: int) -> float:
    """Return estimated file size in MB."""
    if filetype == "h5ad":
        return estimate_h5ad_size_mb(n_cells, n_genes)
    if filetype == "xenium":
        return estimate_xenium_size_mb(n_cells, n_genes)
    return estimate_merscope_size_mb(n_cells, n_genes)


def fmt_size(mb: float) -> str:
    if mb < 1024:
        return f"{mb:.1f} MB"
    return f"{mb / 1024:.2f} GB"


def output_path_for(base_dir: str, filetype: str, n_cells: int, n_genes: int) -> str:
    tag = f"sc_{n_cells}c_{n_genes}g"
    if filetype == "h5ad":
        return os.path.join(base_dir, "h5ad", f"{tag}.h5ad")
    return os.path.join(base_dir, filetype, tag)


# ---------------------------------------------------------------------------
# Sparse matrix generation (chunked for large datasets)
# ---------------------------------------------------------------------------
def generate_sparse_csr_chunked(n_cells, n_genes, density, rng):
    """Build a CSR sparse matrix in chunks to limit peak memory."""
    from scipy.sparse import csr_matrix, vstack as sparse_vstack

    chunks = []
    for start in range(0, n_cells, CHUNK_SIZE):
        end = min(start + CHUNK_SIZE, n_cells)
        chunk_rows = end - start
        nnz = max(1, int(chunk_rows * n_genes * density))
        rows = rng.integers(0, chunk_rows, size=nnz)
        cols = rng.integers(0, n_genes, size=nnz)
        vals = rng.exponential(2.0, size=nnz).astype(np.float32)
        chunk = csr_matrix((vals, (rows, cols)), shape=(chunk_rows, n_genes))
        chunks.append(chunk)
    return sparse_vstack(chunks, format="csr")


# ---------------------------------------------------------------------------
# H5AD generator
# ---------------------------------------------------------------------------
def generate_h5ad(n_cells: int, n_genes: int, n_obs_cols: int, path: str):
    import anndata
    import scipy.sparse

    rng = np.random.default_rng(42)
    print(f"  Generating sparse matrix ({n_cells}×{n_genes}) ...")
    density = get_sparsity(n_cells, n_genes)
    X = generate_sparse_csr_chunked(n_cells, n_genes, density, rng)

    import pandas as pd
    obs_data = make_obs_columns(n_cells, n_obs_cols, rng)
    obs = pd.DataFrame(obs_data, index=[f"cell_{i}" for i in range(n_cells)])
    for col in obs.columns:
        obs[col] = pd.Categorical(obs[col])

    var = pd.DataFrame(index=make_gene_names(n_genes))

    spatial = make_spatial(n_cells, rng)
    umap = make_umap(n_cells, rng)

    adata = anndata.AnnData(
        X=X,
        obs=obs,
        var=var,
        obsm={"X_spatial": spatial, "X_umap": umap},
    )

    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"  Writing {path} ...")
    adata.write_h5ad(path)


# ---------------------------------------------------------------------------
# Xenium generator
# ---------------------------------------------------------------------------
def generate_xenium(n_cells: int, n_genes: int, n_obs_cols: int, folder: str):
    rng = np.random.default_rng(42)
    os.makedirs(folder, exist_ok=True)
    os.makedirs(os.path.join(folder, "cell_feature_matrix"), exist_ok=True)

    gene_names = make_gene_names(n_genes)
    obs_cols = make_obs_columns(n_cells, n_obs_cols, rng)
    spatial = make_spatial(n_cells, rng)

    # --- cells.csv ---
    cells_path = os.path.join(folder, "cells.csv")
    print(f"  Writing {cells_path} ...")
    col_names = ["cell_id", "x_centroid", "y_centroid"] + list(obs_cols.keys())
    with open(cells_path, "w") as f:
        f.write(",".join(col_names) + "\n")
        for start in range(0, n_cells, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_cells)
            lines = []
            for i in range(start, end):
                parts = [str(i), f"{spatial[i, 0]:.4f}", f"{spatial[i, 1]:.4f}"]
                parts += [obs_cols[c][i] for c in obs_cols]
                lines.append(",".join(parts))
            f.write("\n".join(lines) + "\n")

    # --- features.tsv ---
    features_path = os.path.join(folder, "features.tsv")
    with open(features_path, "w") as f:
        f.write("gene\n")
        for g in gene_names:
            f.write(g + "\n")

    # --- matrix.mtx (gzipped) ---
    density = get_sparsity(n_cells, n_genes)
    mtx_path = os.path.join(folder, "cell_feature_matrix", "matrix.mtx.gz")
    print(f"  Writing {mtx_path} ...")

    total_nnz = 0
    # First pass: count nnz
    for start in range(0, n_cells, CHUNK_SIZE):
        end = min(start + CHUNK_SIZE, n_cells)
        chunk_rows = end - start
        total_nnz += max(1, int(chunk_rows * n_genes * density))

    with gzip.open(mtx_path, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"% Generated by MERFISHeyes benchmark pipeline\n")
        # Xenium matrix: genes × cells (transposed from anndata convention)
        f.write(f"{n_genes} {n_cells} {total_nnz}\n")
        for start in range(0, n_cells, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_cells)
            chunk_rows = end - start
            nnz = max(1, int(chunk_rows * n_genes * density))
            rows_idx = rng.integers(0, n_genes, size=nnz)      # gene index
            cols_idx = rng.integers(0, chunk_rows, size=nnz)    # cell index within chunk
            vals = rng.exponential(2.0, size=nnz).astype(np.float32)
            lines = []
            for j in range(nnz):
                # 1-indexed: gene_row cell_col value
                lines.append(f"{rows_idx[j]+1} {cols_idx[j]+start+1} {vals[j]:.4f}")
            f.write("\n".join(lines) + "\n")

    # --- barcodes.tsv ---
    barcodes_path = os.path.join(folder, "cell_feature_matrix", "barcodes.tsv")
    with open(barcodes_path, "w") as f:
        for i in range(n_cells):
            f.write(f"cell_{i}\n")


# ---------------------------------------------------------------------------
# MERSCOPE generator
# ---------------------------------------------------------------------------
def generate_merscope(n_cells: int, n_genes: int, n_obs_cols: int, folder: str):
    rng = np.random.default_rng(42)
    os.makedirs(folder, exist_ok=True)

    gene_names = make_gene_names(n_genes)
    obs_cols = make_obs_columns(n_cells, n_obs_cols, rng)
    spatial = make_spatial(n_cells, rng)

    # --- cell_metadata.csv ---
    meta_path = os.path.join(folder, "cell_metadata.csv")
    print(f"  Writing {meta_path} ...")
    with open(meta_path, "w") as f:
        f.write("EntityID,center_x,center_y\n")
        for start in range(0, n_cells, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_cells)
            lines = []
            for i in range(start, end):
                lines.append(f"{i},{spatial[i,0]:.4f},{spatial[i,1]:.4f}")
            f.write("\n".join(lines) + "\n")

    # --- cell_categories.csv ---
    cat_path = os.path.join(folder, "cell_categories.csv")
    print(f"  Writing {cat_path} ...")
    cat_col_names = ["EntityID"] + list(obs_cols.keys())
    with open(cat_path, "w") as f:
        f.write(",".join(cat_col_names) + "\n")
        for start in range(0, n_cells, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_cells)
            lines = []
            for i in range(start, end):
                parts = [str(i)] + [obs_cols[c][i] for c in obs_cols]
                lines.append(",".join(parts))
            f.write("\n".join(lines) + "\n")

    # --- cell_by_gene.csv (wide format) ---
    cbg_path = os.path.join(folder, "cell_by_gene.csv")
    density = get_sparsity(n_cells, n_genes)
    print(f"  Writing {cbg_path} (wide format, density={density:.4f}) ...")
    with open(cbg_path, "w") as f:
        f.write("cell," + ",".join(gene_names) + "\n")
        for start in range(0, n_cells, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_cells)
            chunk_rows = end - start
            # Generate sparse chunk as dense for CSV writing
            nnz = max(1, int(chunk_rows * n_genes * density))
            rows_idx = rng.integers(0, chunk_rows, size=nnz)
            cols_idx = rng.integers(0, n_genes, size=nnz)
            vals = rng.exponential(2.0, size=nnz).astype(np.float32)

            # Build dense chunk (zeros by default)
            chunk = np.zeros((chunk_rows, n_genes), dtype=np.float32)
            chunk[rows_idx, cols_idx] = vals

            lines = []
            for r in range(chunk_rows):
                row_vals = ",".join(
                    f"{v:.2f}" if v > 0 else "0" for v in chunk[r]
                )
                lines.append(f"{start + r},{row_vals}")
            f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
GENERATORS = {
    "h5ad": generate_h5ad,
    "xenium": generate_xenium,
    "merscope": generate_merscope,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate synthetic single-cell data")
    parser.add_argument("--output-dir", default=os.path.join(DEFAULT_SYNTH_DIR, "single_cell"),
                        help="Output directory")
    parser.add_argument("--cells", type=int, nargs="+", default=SINGLE_CELL_SIZES,
                        help="Cell counts to generate")
    parser.add_argument("--genes", type=int, nargs="+", default=SINGLE_CELL_GENES,
                        help="Gene counts to generate")
    parser.add_argument("--filetypes", nargs="+", default=SINGLE_CELL_FILETYPES,
                        choices=SINGLE_CELL_FILETYPES, help="File types to generate")
    parser.add_argument("--obs-columns", type=int, default=SINGLE_CELL_OBS_COLUMNS,
                        help="Number of observation columns")
    parser.add_argument("--max-file-size-gb", type=float, default=MAX_FILE_SIZE_GB,
                        help="Skip files estimated larger than this (GB)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show plan without generating files")
    args = parser.parse_args()

    combinations = [
        (nc, ng, ft)
        for nc in args.cells
        for ng in args.genes
        for ft in args.filetypes
    ]
    print(f"Benchmark single-cell data generator")
    print(f"  Combinations: {len(combinations)} "
          f"({len(args.cells)} cell sizes × {len(args.genes)} gene counts × {len(args.filetypes)} filetypes)")
    print(f"  Obs columns:  {args.obs_columns}")
    print(f"  Max file size: {args.max_file_size_gb} GB")
    print()

    # Plan / dry-run table
    total_mb = 0.0
    skipped = 0
    plan = []
    for nc, ng, ft in combinations:
        est_mb = estimate_size(ft, nc, ng)
        path = output_path_for(args.output_dir, ft, nc, ng)
        skip = est_mb > args.max_file_size_gb * 1024
        plan.append((nc, ng, ft, est_mb, path, skip))
        if skip:
            skipped += 1
        else:
            total_mb += est_mb

    # Print summary table
    print(f"{'Cells':>12}  {'Genes':>8}  {'Type':>10}  {'Est Size':>12}  {'Status'}")
    print("-" * 70)
    for nc, ng, ft, est_mb, path, skip in plan:
        status = "SKIP (too large)" if skip else "generate"
        print(f"{nc:>12,}  {ng:>8,}  {ft:>10}  {fmt_size(est_mb):>12}  {status}")

    print()
    print(f"Total to generate: {fmt_size(total_mb)}  |  Skipped: {skipped}")

    if args.dry_run:
        print("\n(dry run – no files created)")
        return

    # Generate
    generated = 0
    for nc, ng, ft, est_mb, path, skip in plan:
        if skip:
            continue
        if os.path.exists(path):
            print(f"\n[SKIP] Already exists: {path}")
            generated += 1
            continue

        print(f"\n[{generated+1}/{len(plan)-skipped}] Generating {ft} — "
              f"{nc:,} cells × {ng:,} genes (est {fmt_size(est_mb)})")
        t0 = time.time()
        try:
            GENERATORS[ft](nc, ng, args.obs_columns, path)
            elapsed = time.time() - t0
            actual_mb = 0
            if os.path.isfile(path):
                actual_mb = os.path.getsize(path) / (1024 * 1024)
            elif os.path.isdir(path):
                for root, _, files in os.walk(path):
                    for fname in files:
                        actual_mb += os.path.getsize(os.path.join(root, fname)) / (1024 * 1024)
            print(f"  Done in {elapsed:.1f}s — actual size: {fmt_size(actual_mb)}")
            generated += 1
        except Exception as e:
            print(f"  ERROR: {e}")

    # Write manifest
    manifest = {
        "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "obs_columns": args.obs_columns,
        "datasets": [],
    }
    for nc, ng, ft, est_mb, path, skip in plan:
        if skip or not os.path.exists(path):
            continue
        actual_mb = 0
        if os.path.isfile(path):
            actual_mb = os.path.getsize(path) / (1024 * 1024)
        elif os.path.isdir(path):
            for root, _, files in os.walk(path):
                for fname in files:
                    actual_mb += os.path.getsize(os.path.join(root, fname)) / (1024 * 1024)
        manifest["datasets"].append({
            "n_cells": nc,
            "n_genes": ng,
            "filetype": ft,
            "path": path,
            "file_size_mb": round(actual_mb, 2),
        })

    manifest_path = os.path.join(args.output_dir, "manifest_single_cell.json")
    os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"\nManifest written to {manifest_path}")
    print(f"Generated {generated} files.")


if __name__ == "__main__":
    DEFAULT_SYNTH_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "benchmark_data", "synthetic"
    )
    main()
