"""
Chunked Binary -> H5AD Reverser

Reconstructs an AnnData .h5ad file from the chunked binary format produced by
process_spatial_data.py.

Layout consumed (mirrors process_spatial_data.py output):
    manifest.json
    coords/spatial.bin.gz          [num_points u32, dims u32, float32[...]]
    expr/index.json                {total_genes, num_chunks, chunk_size, genes:[{name,chunk_id,position_in_chunk}]}
    expr/chunk_NNNNN.bin.gz        per-chunk sparse expression (see _read_expr_chunk)
    obs/<col>.json.gz              per-cell list (categorical: strings, "" = missing kept as-is;
                                                  numerical: floats with null)
    obs/metadata.json              {<col>: {type: categorical|numerical, unique_values: int}}
    palettes/<col>.json            (not needed for h5ad reconstruction; ignored)

Reconstruction choices (decided with the dataset owner):
    - Spatial coords  -> adata.obsm['X_spatial']
    - Missing values  -> empty strings "" kept literally (not converted to NaN)
    - obs_names       -> default RangeIndex ('id' kept as a regular obs column)
    - obs columns     -> ALL obs/*.json.gz are included, even ones absent from
                         metadata.json (e.g. barcodeCount, fov), type-inferred.
    - categorical (per metadata.json) -> pandas 'category' dtype
    - numerical (per metadata.json)   -> float (null -> NaN)

Usage:
    python chunked_to_h5ad.py CHUNKED_FOLDER OUTPUT.h5ad
"""

import argparse
import gzip
import json
import struct
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def log(msg, t0=None):
    elapsed = f" [{time.perf_counter() - t0:.1f}s]" if t0 is not None else ""
    print(f"[{time.strftime('%H:%M:%S')}]{elapsed} {msg}", flush=True)


def read_spatial(coords_path: Path) -> np.ndarray:
    """Read coords/spatial.bin.gz -> (num_points, dims) float32 array."""
    with gzip.open(coords_path, "rb") as f:
        raw = f.read()
    num_points, dims = struct.unpack_from("<II", raw, 0)
    coords = np.frombuffer(raw, dtype=np.float32, count=num_points * dims, offset=8)
    return coords.reshape(num_points, dims).copy()


def _read_expr_chunk(chunk_path: Path):
    """
    Parse one expr/chunk_NNNNN.bin.gz.

    Format (from process_spatial_data.write_expression_chunk):
        header        '<IIII' = (1, num_genes, chunk_id, total_cells)
        gene table    num_genes * '<IIIIII'
                          = (gene_idx, offset, data_size, data_size, num_non_zero, 0)
        gene sections concatenated; each section:
            '<II'     = (num_non_zero, num_non_zero)
            uint32[num_non_zero]   row indices
            float32[num_non_zero]  values

    Yields (global_gene_idx, row_indices(int64), values(float32)).
    """
    with gzip.open(chunk_path, "rb") as f:
        buf = f.read()

    _, num_genes, _, _ = struct.unpack_from("<IIII", buf, 0)
    table_start = 16

    for i in range(num_genes):
        (gene_idx, offset, _dsize, _dsize2, num_nz, _reserved) = struct.unpack_from(
            "<IIIIII", buf, table_start + i * 24
        )
        # gene data section: 8-byte header then indices then values
        sec = offset + 8
        # row indices fit in int32 (cell counts << 2^31); halves index memory
        idx = np.frombuffer(buf, dtype=np.uint32, count=num_nz, offset=sec).astype(np.int32)
        vals = np.frombuffer(
            buf, dtype=np.float32, count=num_nz, offset=sec + num_nz * 4
        )
        yield gene_idx, idx, vals.copy()


def build_expression(expr_dir: Path, num_cells: int, t0):
    """Reconstruct the cells x genes sparse matrix (CSR float32) and gene names."""
    index = json.loads((expr_dir / "index.json").read_text())
    total_genes = index["total_genes"]
    # var order = order genes appear in index.json (== global gene_idx)
    gene_names = [g["name"] for g in index["genes"]]
    if len(gene_names) != total_genes:
        log(f"  WARNING: index lists {len(gene_names)} genes, "
            f"total_genes={total_genes}", t0)

    # Collect per-gene columns, then assemble a CSC matrix (one column per gene).
    col_rows = [None] * total_genes
    col_vals = [None] * total_genes

    chunk_files = sorted(expr_dir.glob("chunk_*.bin.gz"))
    log(f"  Reading {len(chunk_files)} expression chunks...", t0)
    for ci, cf in enumerate(chunk_files):
        for gene_idx, idx, vals in _read_expr_chunk(cf):
            col_rows[gene_idx] = idx
            col_vals[gene_idx] = vals
        if (ci + 1) % 50 == 0 or ci == len(chunk_files) - 1:
            log(f"    {ci + 1}/{len(chunk_files)} chunks read", t0)

    # Any gene with no chunk -> empty column
    indptr = np.zeros(total_genes + 1, dtype=np.int64)
    for g in range(total_genes):
        n = 0 if col_rows[g] is None else len(col_rows[g])
        indptr[g + 1] = indptr[g] + n

    total_nnz = int(indptr[-1])
    log(f"  Assembling sparse matrix: {num_cells:,} x {total_genes:,}, "
        f"{total_nnz:,} non-zeros", t0)

    data = np.empty(total_nnz, dtype=np.float32)
    indices = np.empty(total_nnz, dtype=np.int32)
    for g in range(total_genes):
        if col_rows[g] is None:
            continue
        s, e = indptr[g], indptr[g + 1]
        indices[s:e] = col_rows[g]
        data[s:e] = col_vals[g]
        col_rows[g] = None  # free per-gene buffer as we go
        col_vals[g] = None

    mat = sparse.csc_matrix(
        (data, indices, indptr), shape=(num_cells, total_genes)
    )
    return mat.tocsr(), gene_names


def build_obs(obs_dir: Path, num_cells: int, t0) -> pd.DataFrame:
    """Load every obs/*.json.gz into a DataFrame, typed per metadata.json."""
    meta_path = obs_dir / "metadata.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}

    obs = pd.DataFrame(index=pd.RangeIndex(num_cells))
    files = sorted(p for p in obs_dir.glob("*.json.gz"))
    for p in files:
        col = p.name[: -len(".json.gz")]
        with gzip.open(p, "rt", encoding="utf-8") as f:
            values = json.load(f)
        if len(values) != num_cells:
            log(f"  WARNING: obs '{col}' has {len(values):,} values, "
                f"expected {num_cells:,} -- skipping", t0)
            continue

        col_type = meta.get(col, {}).get("type")
        if col_type == "categorical":
            obs[col] = pd.Categorical([str(v) for v in values])
        elif col_type == "numerical":
            obs[col] = pd.to_numeric(pd.Series(values), errors="coerce").astype(
                np.float32
            )
        else:
            # Not in metadata.json (e.g. barcodeCount, fov): infer.
            ser = pd.Series(values)
            num = pd.to_numeric(ser, errors="coerce")
            if num.notna().mean() >= 0.5:
                obs[col] = num.astype(np.float32)
            else:
                obs[col] = pd.Categorical(ser.astype(str))
        log(f"  obs '{col}': {obs[col].dtype} "
            f"({col_type or 'inferred'})", t0)

    return obs


def main():
    parser = argparse.ArgumentParser(
        description="Reconstruct a .h5ad from chunked binary format"
    )
    parser.add_argument("input", type=Path, help="Chunked dataset folder")
    parser.add_argument(
        "output", type=Path, help="Output .h5ad path (must end in .h5ad)"
    )
    args = parser.parse_args()

    if args.output.suffix != ".h5ad":
        print(f"Output path must end in .h5ad: {args.output}")
        sys.exit(1)

    try:
        import anndata as ad
    except ImportError:
        print("anndata is required. Install with: pip install anndata")
        sys.exit(1)

    in_dir: Path = args.input
    if not in_dir.is_dir():
        print(f"Input folder not found: {in_dir}")
        sys.exit(1)

    out_path: Path = args.output

    t0 = time.perf_counter()
    log(f"{'='*60}")
    log(f"Reconstructing h5ad from: {in_dir}")
    log(f"Output: {out_path}")
    log(f"{'='*60}")

    manifest = json.loads((in_dir / "manifest.json").read_text())
    num_cells = manifest["statistics"]["total_cells"]
    log(f"Manifest: {num_cells:,} cells, "
        f"{manifest['statistics']['total_genes']:,} genes, "
        f"type={manifest.get('type')}", t0)

    log("=== Spatial coordinates ===", t0)
    spatial = read_spatial(in_dir / "coords" / "spatial.bin.gz")
    log(f"  Spatial: {spatial.shape}", t0)
    if spatial.shape[0] != num_cells:
        log(f"  WARNING: spatial rows ({spatial.shape[0]:,}) != "
            f"manifest cells ({num_cells:,})", t0)

    log("=== Expression matrix ===", t0)
    X, gene_names = build_expression(in_dir / "expr", num_cells, t0)
    log(f"  Matrix: {X.shape}, nnz={X.nnz:,}", t0)

    log("=== Observation columns ===", t0)
    obs = build_obs(in_dir / "obs", num_cells, t0)
    log(f"  obs: {obs.shape[1]} columns", t0)

    log("=== Assembling AnnData ===", t0)
    var = pd.DataFrame(index=pd.Index(gene_names, name=None))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_spatial"] = spatial
    adata.uns["source_manifest"] = manifest

    log(f"  AnnData: {adata}", t0)
    log("=== Writing h5ad ===", t0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    size_mb = out_path.stat().st_size / 1e6
    log(f"{'='*60}", t0)
    log(f"DONE -> {out_path} ({size_mb:,.1f} MB)", t0)
    log(f"{'='*60}", t0)


if __name__ == "__main__":
    main()
