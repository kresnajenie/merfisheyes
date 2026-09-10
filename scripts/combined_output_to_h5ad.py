#!/usr/bin/env python3
"""Convert a combine_slices_v3.py combined_output/ directory into a query .h5ad
for map_my_cell.py.

obs.index is set to ``{_sample_id}_{raw_id}``, the compound key
process_spatial_data.py's --mmc-csv handling already expects and strips
(reorder_df_to_reference(..., strip_compound=True)) when it matches
mapping_output.csv back against cell_by_gene.csv's plain id column.
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

    metadata_df = pd.read_csv(metadata_path)
    meta_id_col = next(
        (c for c in ("EntityID", "id", "cell_id") if c in metadata_df.columns), None
    )
    if meta_id_col is None:
        raise ValueError(
            "cell_metadata.csv has no recognized ID column (tried EntityID, id, "
            f"cell_id); columns present: {list(metadata_df.columns)}"
        )
    if "_sample_id" not in metadata_df.columns:
        raise ValueError("cell_metadata.csv is missing _sample_id (expected from combine_slices_v3.py)")
    metadata_df[meta_id_col] = metadata_df[meta_id_col].astype(str)

    cbg_df = pd.read_csv(cbg_path)
    cbg_df["cell"] = cbg_df["cell"].astype(str)

    missing = set(metadata_df[meta_id_col]) - set(cbg_df["cell"])
    if missing:
        raise ValueError(
            f"{len(missing):,} metadata cell IDs not found in cell_by_gene.csv "
            f"(sample: {list(missing)[:5]})"
        )

    # Align cbg rows to metadata's order via the plain id -- combine_slices_v3.py
    # doesn't guarantee the two files share row order, only the same id set.
    cbg_df = cbg_df.set_index("cell").loc[metadata_df[meta_id_col].values]
    gene_names = list(cbg_df.columns)
    X = cbg_df.to_numpy(dtype=np.float32)

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
