#!/usr/bin/env python3
"""
Generate synthetic spatial transcriptomics datasets for benchmarking.

Creates H5AD (single cell) and Parquet (single molecule) files at graduated
sizes for cross-platform performance comparison.

Usage:
    python generate_test_data.py --output-dir ./benchmark_data
    python generate_test_data.py --output-dir ./benchmark_data --sizes small medium
    python generate_test_data.py --output-dir ./benchmark_data --type single-cell
"""

import argparse
import json
import os
import time

import numpy as np

SINGLE_CELL_SIZES = {
    "small":  {"cells": 50_000,    "genes": 500},
    "medium": {"cells": 200_000,   "genes": 500},
    "large":  {"cells": 500_000,   "genes": 500},
    "xl":     {"cells": 1_000_000, "genes": 500},
}

SINGLE_MOLECULE_SIZES = {
    "small":  {"molecules": 1_000_000,  "genes": 300},
    "medium": {"molecules": 10_000_000, "genes": 400},
    "large":  {"molecules": 21_000_000, "genes": 500},
    "xl":     {"molecules": 50_000_000, "genes": 500},
}


def generate_h5ad(n_cells: int, n_genes: int, output_path: str) -> dict | None:
    """Generate a synthetic H5AD file with spatial coordinates and clusters."""
    try:
        import anndata
        import pandas as pd
        from scipy.sparse import random as sparse_random
    except ImportError as e:
        print(f"  ERROR: Missing dependency: {e}")
        print("  Install with: pip install anndata scipy pandas")
        return None

    print(f"  Generating H5AD: {n_cells:,} cells x {n_genes} genes...")
    t0 = time.perf_counter()

    # Sparse expression matrix (~10% density, realistic for scRNA-seq)
    X = sparse_random(n_cells, n_genes, density=0.1, format="csr", dtype=np.float32)

    var = pd.DataFrame(index=[f"Gene_{i}" for i in range(n_genes)])

    n_clusters = min(25, max(5, n_cells // 5000))
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(
                np.random.choice([f"Cluster_{i}" for i in range(n_clusters)], n_cells)
            ),
            "celltype": pd.Categorical(
                np.random.choice([f"CellType_{i}" for i in range(15)], n_cells)
            ),
            "n_counts": np.random.exponential(1000, n_cells).astype(np.float32),
        },
        index=[f"Cell_{i}" for i in range(n_cells)],
    )

    spatial = np.column_stack([
        np.random.uniform(0, 10000, n_cells),
        np.random.uniform(0, 10000, n_cells),
    ]).astype(np.float32)
    umap = np.random.randn(n_cells, 2).astype(np.float32) * 10

    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_spatial"] = spatial
    adata.obsm["X_umap"] = umap

    adata.write_h5ad(output_path)
    elapsed = time.perf_counter() - t0
    size_mb = os.path.getsize(output_path) / (1024 * 1024)

    print(f"    -> {size_mb:.1f} MB, generated in {elapsed:.1f}s")
    return {
        "path": os.path.abspath(output_path),
        "cells": n_cells,
        "genes": n_genes,
        "size_mb": round(size_mb, 1),
        "generation_time_s": round(elapsed, 2),
    }


def generate_parquet(n_molecules: int, n_genes: int, output_path: str) -> dict | None:
    """Generate a synthetic Parquet file with molecule coordinates (Xenium column names)."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as e:
        print(f"  ERROR: Missing dependency: {e}")
        print("  Install with: pip install pyarrow")
        return None

    print(f"  Generating Parquet: {n_molecules:,} molecules x {n_genes} genes...")
    t0 = time.perf_counter()

    genes = [f"Gene_{i}" for i in range(n_genes)]
    schema = pa.schema([
        ("feature_name", pa.string()),
        ("x_location", pa.float32()),
        ("y_location", pa.float32()),
        ("z_location", pa.float32()),
    ])

    # Build in chunks to manage memory, then write as single file.
    # Use snappy + Parquet v1 + no dictionary for hyparquet (pure JS) compatibility.
    CHUNK = 5_000_000
    all_tables = []
    remaining = n_molecules

    while remaining > 0:
        chunk_size = min(CHUNK, remaining)
        all_tables.append(pa.table({
            "feature_name": pa.array(
                np.random.choice(genes, chunk_size).tolist(), type=pa.string()
            ),
            "x_location": np.random.uniform(-5000, 5000, chunk_size).astype(np.float32),
            "y_location": np.random.uniform(-5000, 5000, chunk_size).astype(np.float32),
            "z_location": np.random.uniform(-25, 25, chunk_size).astype(np.float32),
        }, schema=schema))
        remaining -= chunk_size
        done = n_molecules - remaining
        print(f"    {100 * done // n_molecules}% ({done:,} / {n_molecules:,})", end="\r")

    combined = pa.concat_tables(all_tables)
    pq.write_table(combined, output_path, compression="snappy",
                   use_dictionary=False, version="1.0")

    elapsed = time.perf_counter() - t0
    size_mb = os.path.getsize(output_path) / (1024 * 1024)

    print(f"    -> {size_mb:.1f} MB, generated in {elapsed:.1f}s              ")
    return {
        "path": os.path.abspath(output_path),
        "molecules": n_molecules,
        "genes": n_genes,
        "size_mb": round(size_mb, 1),
        "generation_time_s": round(elapsed, 2),
    }


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark test datasets")
    parser.add_argument("--output-dir", default="./benchmark_data")
    parser.add_argument(
        "--sizes", nargs="+", default=["small", "medium", "large", "xl"],
        choices=["small", "medium", "large", "xl"],
    )
    parser.add_argument(
        "--type", default="both",
        choices=["single-cell", "single-molecule", "both"],
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    manifest = {"generated_at": time.strftime("%Y-%m-%d %H:%M:%S"), "datasets": []}

    if args.type in ("single-cell", "both"):
        print("\n=== Single Cell Datasets (H5AD) ===")
        sc_dir = os.path.join(args.output_dir, "single_cell")
        os.makedirs(sc_dir, exist_ok=True)
        for size in args.sizes:
            cfg = SINGLE_CELL_SIZES[size]
            path = os.path.join(sc_dir, f"benchmark_{size}_{cfg['cells']}_cells.h5ad")
            result = generate_h5ad(cfg["cells"], cfg["genes"], path)
            if result:
                result["type"] = "single_cell"
                result["size_tier"] = size
                manifest["datasets"].append(result)

    if args.type in ("single-molecule", "both"):
        print("\n=== Single Molecule Datasets (Parquet) ===")
        sm_dir = os.path.join(args.output_dir, "single_molecule")
        os.makedirs(sm_dir, exist_ok=True)
        for size in args.sizes:
            cfg = SINGLE_MOLECULE_SIZES[size]
            path = os.path.join(sm_dir, f"benchmark_{size}_{cfg['molecules']}_molecules.parquet")
            result = generate_parquet(cfg["molecules"], cfg["genes"], path)
            if result:
                result["type"] = "single_molecule"
                result["size_tier"] = size
                manifest["datasets"].append(result)

    manifest_path = os.path.join(args.output_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"\nDone! Manifest saved to {manifest_path}")
    print(f"Generated {len(manifest['datasets'])} datasets")


if __name__ == "__main__":
    main()
