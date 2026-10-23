#!/usr/bin/env python3
"""Convert a combine_slices_v3.py combined_output/ directory into a query .h5ad
for map_my_cell.py.

obs.index is set to ``{_sample_id}_{raw_id}``, the compound key
process_spatial_data.py's --mmc-csv handling already expects and strips
(reorder_df_to_reference(..., strip_compound=True)) when it matches
mapping_output.csv back against cell_by_gene.csv's plain id column.

Streams cell_by_gene.csv in chunks straight into one pre-allocated float32
array instead of round-tripping through full-size pandas copies -- a
naive read + to_numpy() + set_index().loc[] reorder held 2-3 float64/float32
copies of the whole matrix alive at once and OOM-killed a 512GB job on an
8.6M-cell dataset.
"""
import argparse
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd


def convert(combined_dir: Path, output_path: Path):
    metadata_path = combined_dir / "cell_metadata.csv"
    cbg_path = combined_dir / "cell_by_gene.csv"
    if not metadata_path.exists():
        raise FileNotFoundError(f"cell_metadata.csv not found in {combined_dir}")
    if not cbg_path.exists():
        raise FileNotFoundError(f"cell_by_gene.csv not found in {combined_dir}")

    meta_header = pd.read_csv(metadata_path, nrows=0).columns.tolist()
    meta_id_col = next(
        (c for c in ("EntityID", "id", "cell_id") if c in meta_header), None
    )
    if meta_id_col is None:
        raise ValueError(
            "cell_metadata.csv has no recognized ID column (tried EntityID, id, "
            f"cell_id); columns present: {meta_header}"
        )
    if "_sample_id" not in meta_header:
        raise ValueError("cell_metadata.csv is missing _sample_id (expected from combine_slices_v3.py)")

    metadata_df = pd.read_csv(metadata_path, usecols=[meta_id_col, "_sample_id"], dtype=str)
    n_cells = len(metadata_df)
    print(f"Metadata: {n_cells:,} cells")

    cbg_header = pd.read_csv(cbg_path, nrows=0).columns.tolist()
    gene_names = [c for c in cbg_header if c != "cell"]
    n_genes = len(gene_names)
    print(f"Cell-by-gene: {n_genes:,} genes")

    # target row for each cell id, per metadata's order (defines obs order)
    row_of = {cid: i for i, cid in enumerate(metadata_df[meta_id_col].values)}

    X = np.zeros((n_cells, n_genes), dtype=np.float32)
    filled = np.zeros(n_cells, dtype=bool)

    rows_seen = 0
    for chunk in pd.read_csv(cbg_path, chunksize=50_000, dtype={"cell": str}):
        rows_seen += len(chunk)
        target = chunk["cell"].map(row_of)
        matched = target.notna()
        if not matched.all():
            chunk = chunk[matched]
            target = target[matched]
        rows = target.to_numpy(dtype=np.int64)
        X[rows, :] = chunk[gene_names].to_numpy(dtype=np.float32)
        filled[rows] = True
        if rows_seen % 500_000 < 50_000:
            print(f"  ...{rows_seen:,} cell_by_gene rows read")

    missing = int((~filled).sum())
    if missing:
        raise ValueError(
            f"{missing:,} metadata cell IDs not found in cell_by_gene.csv"
        )

    compound_key = metadata_df["_sample_id"].astype(str) + "_" + metadata_df[meta_id_col]
    obs = pd.DataFrame(index=pd.Index(compound_key.values, name=None))
    var = pd.DataFrame(index=pd.Index(gene_names, name=None))
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Wrote {output_path} ({adata.shape[0]:,} cells x {adata.shape[1]:,} genes)")
    return output_path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("combined_dir", type=Path, help="Path to combined_output/ directory")
    p.add_argument("output", type=Path, help="Output .h5ad path")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not args.combined_dir.is_dir():
        print(f"ERROR: not a directory: {args.combined_dir}", file=sys.stderr)
        sys.exit(1)
    convert(args.combined_dir, args.output)
