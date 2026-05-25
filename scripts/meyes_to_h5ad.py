"""Convert a meyes_output/ directory back into an .h5ad file.

Reverses process_spatial_data.py:
  coords/spatial.bin.gz   -> obsm['X_spatial']
  expr/index.json + chunk_*.bin.gz -> adata.X (sparse CSC, float32)
  obs/*.json.gz           -> adata.obs columns (categorical/numerical per metadata.json)
  palettes/*.json         -> adata.uns['{col}_colors']
"""

import argparse
import gzip
import json
import struct
import sys
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_spatial(coords_dir: Path) -> np.ndarray:
    path = coords_dir / "spatial.bin.gz"
    with gzip.open(path, "rb") as f:
        buf = f.read()
    n, d = struct.unpack_from("<II", buf, 0)
    coords = np.frombuffer(buf, dtype=np.float32, offset=8, count=n * d).reshape(n, d)
    return np.array(coords)


def read_chunk(chunk_path: Path):
    """Yield (gene_idx, indices, values) for each gene in a chunk."""
    with gzip.open(chunk_path, "rb") as f:
        buf = f.read()
    _, num_genes, _, _ = struct.unpack_from("<IIII", buf, 0)
    gene_table_start = 16
    for i in range(num_genes):
        entry = gene_table_start + i * 24
        gene_idx, offset, _data_size, _ds2, num_nz, _ = struct.unpack_from("<IIIIII", buf, entry)
        # sparse gene section: 8-byte header (n, n), then indices, then values
        n1, n2 = struct.unpack_from("<II", buf, offset)
        if n1 != num_nz or n2 != num_nz:
            raise ValueError(
                f"chunk {chunk_path.name} gene {gene_idx}: header mismatch "
                f"(n1={n1}, n2={n2}, expected num_nz={num_nz})"
            )
        idx_start = offset + 8
        val_start = idx_start + num_nz * 4
        indices = np.frombuffer(buf, dtype=np.uint32, offset=idx_start, count=num_nz)
        values = np.frombuffer(buf, dtype=np.float32, offset=val_start, count=num_nz)
        yield int(gene_idx), np.array(indices), np.array(values)


def read_expression(expr_dir: Path, num_cells: int) -> tuple[sp.csc_matrix, list[str]]:
    with open(expr_dir / "index.json") as f:
        index = json.load(f)
    total_genes = index["total_genes"]
    gene_names = [None] * total_genes
    for g in index["genes"]:
        gene_names[g["chunk_id"] * index["chunk_size"] + g["position_in_chunk"]] = g["name"]
    # Assemble CSC by collecting indptr, indices, data column-by-column.
    indptr = np.zeros(total_genes + 1, dtype=np.int64)
    cols_data = [None] * total_genes
    cols_idx = [None] * total_genes

    chunk_files = sorted(expr_dir.glob("chunk_*.bin.gz"))
    for cf in chunk_files:
        for gene_idx, indices, values in read_chunk(cf):
            cols_idx[gene_idx] = indices.astype(np.int32, copy=False)
            cols_data[gene_idx] = values

    for j in range(total_genes):
        if cols_idx[j] is None:
            cols_idx[j] = np.empty(0, dtype=np.int32)
            cols_data[j] = np.empty(0, dtype=np.float32)
        indptr[j + 1] = indptr[j] + cols_idx[j].size

    indices_arr = np.concatenate(cols_idx) if total_genes else np.empty(0, dtype=np.int32)
    data_arr = np.concatenate(cols_data) if total_genes else np.empty(0, dtype=np.float32)
    X = sp.csc_matrix((data_arr, indices_arr, indptr), shape=(num_cells, total_genes))
    return X, gene_names


def read_obs(obs_dir: Path, num_cells: int) -> tuple[pd.DataFrame, dict]:
    with open(obs_dir / "metadata.json") as f:
        meta = json.load(f)
    cols = {}
    for col_name, info in meta.items():
        path = obs_dir / f"{col_name}.json.gz"
        if not path.exists():
            print(f"  warn: missing {path.name}, skipping")
            continue
        with gzip.open(path, "rt") as f:
            vals = json.load(f)
        if len(vals) != num_cells:
            print(f"  warn: {col_name} has {len(vals)} values, expected {num_cells}; skipping")
            continue
        series = pd.Series(vals, name=col_name)
        if info["type"] == "categorical":
            series = series.astype("category")
        else:
            series = pd.to_numeric(series, errors="coerce")
        cols[col_name] = series
    df = pd.DataFrame(cols)
    return df, meta


def read_palettes(palettes_dir: Path) -> dict:
    palettes = {}
    if not palettes_dir.exists():
        return palettes
    for p in palettes_dir.glob("*.json"):
        with open(p) as f:
            palettes[p.stem] = json.load(f)
    return palettes


def attach_palettes(adata: anndata.AnnData, palettes: dict):
    for col, palette in palettes.items():
        if col not in adata.obs.columns:
            continue
        if not isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            continue
        categories = list(adata.obs[col].cat.categories)
        if isinstance(palette, dict):
            colors = [palette.get(str(c), "#888888") for c in categories]
        elif isinstance(palette, list):
            colors = [palette[i] if i < len(palette) else "#888888" for i in range(len(categories))]
        else:
            continue
        adata.uns[f"{col}_colors"] = colors


def convert(input_dir: Path, output_path: Path):
    print(f"Reading {input_dir}")
    with open(input_dir / "manifest.json") as f:
        manifest = json.load(f)
    num_cells = manifest["statistics"]["total_cells"]
    total_genes = manifest["statistics"]["total_genes"]
    print(f"  cells={num_cells}, genes={total_genes}")

    print("Reading spatial coordinates...")
    spatial = read_spatial(input_dir / "coords")
    print(f"  spatial shape: {spatial.shape}, dtype: {spatial.dtype}")

    print("Reading expression chunks...")
    X, gene_names = read_expression(input_dir / "expr", num_cells)
    print(f"  X shape: {X.shape}, nnz: {X.nnz}, dtype: {X.dtype}")

    print("Reading obs columns...")
    obs, _meta = read_obs(input_dir / "obs", num_cells)
    print(f"  obs columns: {len(obs.columns)}")

    print("Reading palettes...")
    palettes = read_palettes(input_dir / "palettes")
    print(f"  palettes: {len(palettes)}")

    # Build AnnData. Use integer string indices since the original cell ids are
    # not preserved by process_spatial_data.py (they may live in obs['id']).
    obs.index = obs.index.astype(str)
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))

    adata = anndata.AnnData(X=X.tocsr(), obs=obs, var=var)
    adata.obsm["X_spatial"] = spatial
    attach_palettes(adata, palettes)
    adata.uns["meyes_manifest"] = manifest

    print(f"Writing {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")
    print(f"Done. {output_path} ({output_path.stat().st_size / 1e6:.1f} MB)")
    return adata


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir", type=Path, help="Path to meyes_output/ directory")
    p.add_argument("output", type=Path, help="Output .h5ad path")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if not args.input_dir.is_dir():
        print(f"ERROR: not a directory: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    convert(args.input_dir, args.output)
