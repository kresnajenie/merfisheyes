#!/usr/bin/env python3
"""
Leiden clustering for a SpatialData store

Clusters the cells of a SpatialData zarr store's table. The clusters go into
the table itself, as a categorical `leiden` column of `obs`, and so do the
per-cluster differential expression stats the viewer's DEG panel reads
(`varm/de_mean__leiden`, `varm/de_pct__leiden`, `uns/merfisheyes_de`; see
scripts/de_stats.py), so any reader of the store sees them.

Usage:
    python spatialdata_cluster.py store.zarr annotations/
    python spatialdata_cluster.py store.zarr annotations/ --resolution 0.8
    python spatialdata_cluster.py store.zarr annotations/ --sidecar-only   # leave the table alone

Output:
    store.zarr/tables/<table>/obs/leiden     the cluster of each cell (categorical)
    store.zarr/tables/<table>/varm, uns      its DE stats
    annotations/             only with --sidecar-only, the older layout:
    ├── annotations.json     the columns kept here with their categories
    ├── de/leiden.bin.gz     per-cluster mean / pct expressing per gene
    └── leiden.bin           uint16 category code per cell

Clustering: normalize_total → log1p → PCA → neighbours → Leiden. DE stats are
means of that log-normalized matrix, which is what the viewer's DEG panel
expects (it un-logs them for fold change) — raw Xenium counts are low enough
that it would mistake their means for log values.
"""

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import zarr

sys.path.insert(0, str(Path(__file__).resolve().parent))

from process_spatial_data import write_de_stats_binary  # noqa: E402


def _drop_consolidated(zarr_json: Path) -> None:
    """anndata refuses to write into a group whose zarr.json carries
    consolidated metadata; it is stale once a column is added anyway."""
    try:
        meta = json.loads(zarr_json.read_text())
    except (OSError, ValueError):
        return
    if "consolidated_metadata" in meta:
        del meta["consolidated_metadata"]
        zarr_json.write_text(json.dumps(meta, indent=2))


def write_column(store, table: str, name: str, values: pd.Categorical) -> None:
    """A categorical obs column into `tables/<table>` of `store`, listed in
    the dataframe's column order so anndata and the viewer both see it."""
    store = Path(store)
    for p in (store / "zarr.json", store / "tables" / table / "zarr.json", store / "tables" / table / "obs" / "zarr.json"):
        _drop_consolidated(p)
    group = zarr.open_group(str(store), mode="r+", use_consolidated=False)[f"tables/{table}"]
    obs = group["obs"]
    ad.io.write_elem(obs, name, values)
    order = [c for c in list(obs.attrs.get("column-order", [])) if c != name] + [name]
    obs.attrs["column-order"] = order


def cluster(
    store,
    output_folder,
    *,
    table: str = "table",
    resolution: float = 0.5,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    into_table: bool = True,
    de: bool = True,
    log=print,
) -> Path | None:
    """Cluster `tables/<table>` of `store`: the clusters into the table's obs,
    with their DE stats in the table (`de=False` leaves those to the caller);
    or, with `into_table=False`, both into the sidecar, the DE stats into the
    `annotations/` sidecar."""
    # Opened by path: these stores' consolidated metadata can be stale
    group = zarr.open_group(str(store), mode="r", use_consolidated=False)[f"tables/{table}"]
    counts = sp.csr_matrix(ad.io.read_elem(group["X"]))
    genes = [str(g) for g in ad.io.read_elem(group["var"]).index]
    log(f"Table: {counts.shape[0]:,} cells x {counts.shape[1]:,} genes")

    adata = ad.AnnData(counts.astype(np.float32))
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    lognorm = adata.X.copy()
    # PCA needs fewer components than cells and genes; tiny tables are allowed
    n_comps = max(1, min(n_pcs, counts.shape[0] - 1, counts.shape[1] - 1))
    sc.pp.pca(adata, n_comps=n_comps)
    sc.pp.neighbors(adata, n_neighbors=min(n_neighbors, counts.shape[0] - 1))
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2)
    codes = adata.obs["leiden"].cat.codes.to_numpy()
    categories = [str(c) for c in adata.obs["leiden"].cat.categories]
    cell_counts = np.bincount(codes, minlength=len(categories))
    log(f"Leiden (resolution {resolution}): {len(categories)} clusters")

    # Per-cluster sums via a cells x clusters membership matrix
    membership = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(categories)),
    )
    sums = np.asarray((lognorm.T @ membership).todense())
    expressing = np.asarray(((counts > 0).astype(np.float32).T @ membership).todense())

    if into_table:
        # The clusters into the table, and their DE stats beside them in it
        # (scripts/de_stats.py); no folder beside the store
        write_column(store, table, "leiden", pd.Categorical.from_codes(codes, categories))
        log(f"Wrote tables/{table}/obs/leiden")
        if de:
            import de_stats

            group = zarr.open_group(str(store), mode="r+", use_consolidated=False)[f"tables/{table}"]
            de_stats.write_de_stats(group, log, columns=["leiden"])
            # Writing dropped the root's consolidated metadata; a store opened
            # by URL is listed by it (prepare_folder.py writes it once, last)
            zarr.consolidate_metadata(str(store))
        return None

    out = Path(output_folder)
    (out / "de").mkdir(parents=True, exist_ok=True)
    codes.astype("<u2").tofile(out / "leiden.bin")
    columns = {"leiden": {"categories": categories}}
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
                "columns": columns,
                "de_stats": ["leiden"],
            },
            indent=2,
        )
    )
    log(f"Wrote {out}/")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("store", help="SpatialData zarr store")
    parser.add_argument("output_folder", help="annotations folder to write")
    parser.add_argument("--table", default="table")
    parser.add_argument("--resolution", type=float, default=0.5)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--sidecar-only", action="store_true",
                        help="keep the clusters in annotations/ too, and leave the table untouched")
    args = parser.parse_args()
    cluster(
        args.store,
        args.output_folder,
        table=args.table,
        resolution=args.resolution,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        into_table=not args.sidecar_only,
    )


if __name__ == "__main__":
    main()
