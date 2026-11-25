#!/usr/bin/env python3
"""mh1: a count-matrix hash that names the same data the same, however it is stored.

    python scripts/metadata/matrix_hash.py STORE.zarr [--table table] [--layer NAME] [--by _sample_id]

The reference for docs/specs/metadata.md §7.7. Equal hashes mean the same cells with the
same counts for the same genes, whatever the chunking, compression, sparse format (CSR,
CSC, dense), dtype width, cell order or gene order. Prints the store's hash, the provenance
block of meyes_provenance.matrix_hash and, with --by, one hash per section.

The algorithm, version mh1:

  1. Genes: the table's feature names, stripped and lower-cased; they must be unique.
  2. Values: each entry rounded half-even to 6 decimals; an integer is written without a
     fraction ("3"), anything else as Python's shortest repr ("0.25"); zeros are left out.
  3. Each cell: its non-zero (gene, value) pairs sorted by gene, written "gene\\tvalue\\n",
     UTF-8; the cell's digest is the first 16 bytes of their SHA-256. A cell with no
     counts hashes the empty string.
  4. The matrix: the cell digests sorted, concatenated after the bytes "mh1\\n", SHA-256.
     The hash is "mh1:" + its hex.
  5. The sketch: the 64 smallest cell digests' first 8 bytes, hex, sorted.

Cell names, obs and var annotations, coordinates and genes with no counts do not enter it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys

import numpy as np
import scipy.sparse as sp

VERSION = "mh1"
SKETCH = 64


def _token(v: float) -> str:
    t = round(float(v), 6)
    return str(int(t)) if t.is_integer() else repr(t)


def cell_digests(matrix, genes) -> list[bytes]:
    """One 16-byte digest per row of `matrix` (cells × genes)."""
    names = [str(g).strip().lower() for g in genes]
    if len(set(names)) != len(names):
        raise ValueError("gene names are not unique once stripped and lower-cased")
    m = sp.csr_matrix(matrix, dtype=np.float64)
    order = np.argsort(np.asarray(names, dtype=object), kind="stable")   # columns by name
    m = m[:, order]
    m.sort_indices()
    names = [names[i] for i in order]
    out = []
    for r in range(m.shape[0]):
        lo, hi = m.indptr[r], m.indptr[r + 1]
        parts = []
        for j, v in zip(m.indices[lo:hi], m.data[lo:hi]):
            tok = _token(v)
            if tok != "0":
                parts.append(f"{names[j]}\t{tok}\n")
        out.append(hashlib.sha256("".join(parts).encode()).digest()[:16])
    return out


def combine(digests: list[bytes]) -> tuple[str, list[str]]:
    """The matrix hash and its sketch, from cell digests in any order."""
    d = sorted(digests)
    h = hashlib.sha256(b"mh1\n" + b"".join(d)).hexdigest()
    return f"{VERSION}:{h}", [x[:8].hex() for x in d[:SKETCH]]


def matrix_hash(matrix, genes) -> str:
    return combine(cell_digests(matrix, genes))[0]


def sketch_overlap(a: list[str], b: list[str]) -> float:
    """Bottom-k estimate of the Jaccard similarity of two stores' cell sets."""
    k = min(len(a), len(b))
    if not k:
        return 0.0
    union = sorted(set(a) | set(b))[:k]
    return len(set(union) & set(a) & set(b)) / k


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("store")
    ap.add_argument("--table", default="table")
    ap.add_argument("--layer", help="hash this layer instead of X (when X is not counts)")
    ap.add_argument("--by", help="obs column naming sections: one hash per section too")
    a = ap.parse_args()

    import anndata as ad
    import time

    t0 = time.time()
    adata = ad.read_zarr(f"{a.store.rstrip('/')}/tables/{a.table}")
    m = adata.layers[a.layer] if a.layer else adata.X
    genes = list(adata.var_names)
    digests = cell_digests(m, genes)
    h, sketch = combine(digests)
    nnz = int(sp.csr_matrix(m).count_nonzero())
    out = {"meyes_matrix_hash": h,
           "matrix_hash": {"version": VERSION, "matrix": f"layers/{a.layer}" if a.layer else "X",
                           "genes": "symbol", "n_cells": len(digests), "n_nonzero": nnz, "sketch": sketch}}
    if a.by:
        col = adata.obs[a.by].astype(str).to_numpy()
        out["sections"] = {s: combine([d for d, c in zip(digests, col) if c == s])[0]
                           for s in sorted(set(col))}
    out["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
