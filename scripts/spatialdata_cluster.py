#!/usr/bin/env python3
"""
Leiden clustering for a SpatialData store, written as a sidecar

Clusters the cells of a SpatialData zarr store's table and writes the result
next to the store as an `annotations/` folder — the store itself is only
read. The SpatialData viewer picks the folder up as extra cell columns with
precomputed differential expression stats.

Usage:
    python spatialdata_cluster.py store.zarr annotations/
    python spatialdata_cluster.py store.zarr annotations/ --resolution 0.8

Output:
    annotations/
    ├── annotations.json     columns, their categories, which have DE stats
    ├── leiden.bin           uint16 category code per cell, in table order
    └── de/leiden.bin.gz     per-cluster mean / pct expressing per gene

Clustering: normalize_total → log1p → PCA → neighbours → Leiden. DE stats are
means of that log-normalized matrix, which is what the viewer's DEG panel
expects (it un-logs them for fold change) — raw Xenium counts are low enough
that it would mistake their means for log values.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import zarr

from process_spatial_data import write_de_stats_binary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("store", help="SpatialData zarr store")
    parser.add_argument("output_folder", help="annotations folder to write")
    parser.add_argument("--table", default="table")
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    args = parser.parse_args()

    # Opened by path: these stores' consolidated metadata can be stale
    table = zarr.open_group(args.store, mode="r", use_consolidated=False)[
        f"tables/{args.table}"
    ]
    counts = sp.csr_matrix(ad.io.read_elem(table["X"]))
    genes = [str(g) for g in ad.io.read_elem(table["var"]).index]
    print(f"Table: {counts.shape[0]:,} cells x {counts.shape[1]:,} genes")

    adata = ad.AnnData(counts.astype(np.float32))
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    lognorm = adata.X.copy()
    sc.pp.pca(adata, n_comps=args.n_pcs)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors)
    sc.tl.leiden(
        adata, resolution=args.resolution, flavor="igraph", n_iterations=2
    )
    codes = adata.obs["leiden"].cat.codes.to_numpy()
    categories = [str(c) for c in adata.obs["leiden"].cat.categories]
    cell_counts = np.bincount(codes, minlength=len(categories))
    print(f"Leiden (resolution {args.resolution}): {len(categories)} clusters")

    # Per-cluster sums via a cells x clusters membership matrix
    membership = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(categories)),
    )
    sums = np.asarray((lognorm.T @ membership).todense())
    expressing = np.asarray(((counts > 0).astype(np.float32).T @ membership).todense())

    out = Path(args.output_folder)
    (out / "de").mkdir(parents=True, exist_ok=True)
    codes.astype("<u2").tofile(out / "leiden.bin")
    write_de_stats_binary(
        categories,
        cell_counts,
        sums / cell_counts,
        expressing / cell_counts,
        out / "de" / "leiden.bin.gz",
    )
    (out / "annotations.json").write_text(
        json.dumps(
            {
                "num_cells": len(codes),
                "num_genes": len(genes),
                "columns": {"leiden": {"categories": categories}},
                "de_stats": ["leiden"],
            },
            indent=2,
        )
    )
    print(f"Wrote {out}/")


if __name__ == "__main__":
    main()
