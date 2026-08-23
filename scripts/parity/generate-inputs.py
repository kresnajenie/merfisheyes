#!/usr/bin/env python3
"""Deterministic inputs for the pipeline parity harness (docs/PARITY-TESTING.md).

Small datasets that exercise the rules both pipelines must agree on: control-probe
gene filtering, categorical obs with missing labels, numerical/bool obs, an
all-empty obs column, uns[<col>_colors], obsm embeddings wider than 3 dims,
coordinates falling back to obs columns, sparse X, Xenium and MERSCOPE folders,
and single-molecule parquet / CSV with and without a cell-id column (incl. a gene
whose molecules are all unassigned, NaN coordinates and empty gene names).

Usage: python3 scripts/parity/generate-inputs.py <out-dir>
Writes <out-dir>/cases.json listing every case (read by run.py / run-js.mts).
"""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "testdata"))
from generate import CELLTYPES  # noqa: E402

# 34 real genes + 6 that every gene filter must drop (name-based, both pipelines).
GENES = [f"gene_{i}" for i in range(34)] + [
    "Blank-1",
    "Blank-2",
    "NegControlProbe_00001",
    "NegControlCodeword_0500",
    "UnassignedCodeword_0001",
    "Deprecated_0001",
]
FEATURE_TYPES = {
    "Blank-1": "Negative Control Codeword",
    "Blank-2": "Negative Control Codeword",
    "NegControlProbe_00001": "Negative Control Probe",
    "NegControlCodeword_0500": "Negative Control Codeword",
    "UnassignedCodeword_0001": "Unassigned Codeword",
    "Deprecated_0001": "Deprecated Codeword",
}
N_CELLS = 300
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
           "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


def sc_arrays(seed: int, dims: int):
    rng = np.random.default_rng(seed)
    n, g = N_CELLS, len(GENES)
    coords = np.round((rng.random((n, dims)) - 0.5) * 2000.0, 4)
    expr = rng.random((n, g))
    expr[expr < 0.85] = 0.0
    expr = np.round(expr * 3.0, 2)
    leiden = np.array([f"cluster_{rng.integers(0, 8)}" for _ in range(n)], dtype=object)
    leiden[rng.random(n) < 0.03] = None  # missing labels
    celltype = np.array([CELLTYPES[rng.integers(0, len(CELLTYPES))] for _ in range(n)], dtype=object)
    n_counts = np.round(rng.random(n) * 500, 3)
    n_counts[rng.random(n) < 0.02] = np.nan
    return dict(
        coords=coords,
        expr=expr,
        leiden=leiden,
        celltype=celltype,
        n_counts=n_counts,
        n_genes=rng.integers(5, 40, n).astype(np.int32),
        is_doublet=rng.random(n) < 0.1,
        sample=np.where(rng.random(n) < 0.5, "sample_A", "sample_B"),
        score=np.round(rng.standard_normal(n), 5),
        umap=np.round((rng.random((n, 2)) - 0.5) * 20.0, 4),
        pca=np.round(rng.standard_normal((n, 50)), 4),
    )


def obs_frame(a, with_coords=False):
    obs = pd.DataFrame(
        {
            "leiden": pd.Categorical(a["leiden"]),
            "celltype": pd.Categorical(a["celltype"]),
            "n_counts": a["n_counts"],
            "n_genes": a["n_genes"],
            "is_doublet": a["is_doublet"],
            "sample": a["sample"],
            "score": a["score"],
            "empty_col": np.full(N_CELLS, np.nan),
        },
        index=[f"cell_{i}" for i in range(N_CELLS)],
    )
    if with_coords:
        obs["center_x"] = a["coords"][:, 0]
        obs["center_y"] = a["coords"][:, 1]
    return obs


def write_h5ad(path: Path, seed: int, *, sparse_x=False, dims=2, coords_in_obs=False):
    import anndata as ad
    from scipy import sparse

    a = sc_arrays(seed, dims)
    obs = obs_frame(a, with_coords=coords_in_obs)
    X = sparse.csr_matrix(a["expr"].astype(np.float32)) if sparse_x else a["expr"]
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=GENES))
    if not coords_in_obs:
        adata.obsm["X_spatial"] = a["coords"].astype(np.float32) if sparse_x else a["coords"]
    adata.obsm["X_umap"] = a["umap"].astype(np.float32)
    adata.obsm["X_pca"] = a["pca"].astype(np.float32)
    adata.uns["leiden_colors"] = np.array(PALETTE[: len(obs["leiden"].cat.categories)], dtype=object)
    adata.uns["celltype_colors"] = np.array(PALETTE[: len(obs["celltype"].cat.categories)], dtype=object)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path


def _label(v):
    return "" if v is None else str(v)


def write_xenium(out: Path, seed: int):
    a = sc_arrays(seed, 2)
    out.mkdir(parents=True, exist_ok=True)
    cell_ids = [f"aaaa{i:05d}-1" for i in range(N_CELLS)]
    with open(out / "cells.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell_id", "x_centroid", "y_centroid", "transcript_counts",
                    "cell_area", "nucleus_area", "leiden", "celltype"])
        for i in range(N_CELLS):
            w.writerow([cell_ids[i], a["coords"][i, 0], a["coords"][i, 1],
                        int(a["n_genes"][i]), float(a["n_counts"][i]) if np.isfinite(a["n_counts"][i]) else "",
                        round(float(a["score"][i]) * 10 + 50, 3), _label(a["leiden"][i]), a["celltype"][i]])
    mtx = out / "cell_feature_matrix"
    mtx.mkdir(exist_ok=True)
    # Real Xenium features.tsv: no header, id \t name \t feature type.
    with open(mtx / "features.tsv", "w") as f:
        for i, g in enumerate(GENES):
            f.write(f"ENSG{i:011d}\t{g}\t{FEATURE_TYPES.get(g, 'Gene Expression')}\n")
    with open(mtx / "barcodes.tsv", "w") as f:
        f.write("\n".join(cell_ids) + "\n")
    nz = np.argwhere(a["expr"].T > 0)  # (gene, cell), 1-based in the file
    with gzip.open(mtx / "matrix.mtx.gz", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate real general\n%\n")
        f.write(f"{len(GENES)} {N_CELLS} {len(nz)}\n")
        for g, c in nz:
            f.write(f"{g + 1} {c + 1} {a['expr'][c, g]}\n")
    return out


def write_merscope(out: Path, seed: int):
    a = sc_arrays(seed, 2)
    out.mkdir(parents=True, exist_ok=True)
    ids = [str(10_000_000_000 + i) for i in range(N_CELLS)]
    with open(out / "cell_metadata.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["EntityID", "fov", "volume", "center_x", "center_y",
                    "min_x", "min_y", "max_x", "max_y", "leiden", "celltype"])
        for i in range(N_CELLS):
            x, y = a["coords"][i, 0], a["coords"][i, 1]
            w.writerow([ids[i], int(i % 12), round(float(a["n_counts"][i]) if np.isfinite(a["n_counts"][i]) else 0.0, 3),
                        x, y, x - 5, y - 5, x + 5, y + 5, _label(a["leiden"][i]), a["celltype"][i]])
    with open(out / "cell_by_gene.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell"] + GENES)
        for i in range(N_CELLS):
            w.writerow([ids[i]] + [a["expr"][i, j] if a["expr"][i, j] else 0 for j in range(len(GENES))])
    return out


def _sm_frame(seed: int, dims: int, gene_col: str, xyz, cell_id: bool):
    rng = np.random.default_rng(seed + 7)
    n = 6000
    genes = np.array(GENES, dtype=object)[rng.integers(0, len(GENES), n)]
    x = np.round((rng.random(n) - 0.5) * 2000.0, 4)
    y = np.round((rng.random(n) - 0.5) * 2000.0, 4)
    z = np.round((rng.random(n) - 0.5) * 50.0, 4)
    # Edge rows: non-finite coordinates and empty gene names must be dropped.
    x[:3] = np.nan
    y[3] = np.inf
    genes[4:6] = ""
    cols = {gene_col: genes, xyz[0]: x, xyz[1]: y}
    if dims == 3:
        cols[xyz[2]] = z
    if cell_id:
        cid = rng.integers(1, 400, n)
        cid[rng.random(n) < 0.2] = -1  # unassigned
        cid[genes == "gene_33"] = -1  # a gene whose molecules are ALL unassigned
        cols["cell_id"] = cid
        cols["fov"] = rng.integers(0, 12, n)
    return pd.DataFrame(cols)


def write_sm_parquet(path: Path, seed: int, *, preset: str, dims: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    if preset == "xenium":
        df = _sm_frame(seed, dims, "feature_name", ("x_location", "y_location", "z_location"), cell_id=False)
    else:
        df = _sm_frame(seed, dims, "gene", ("global_x", "global_y", "global_z"), cell_id=True)
    df.to_parquet(path, index=False)
    return path


def write_sm_csv(path: Path, seed: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = _sm_frame(seed, 3, "gene", ("global_x", "global_y", "global_z"), cell_id=True)
    df.to_csv(path, index=False)
    return path


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("out", type=Path)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    out, seed = args.out, args.seed
    cases = []

    def sc(name, fmt, path):
        cases.append({"name": name, "kind": "sc", "format": fmt, "input": str(path.relative_to(out))})

    def sm(name, fmt, path, dataset_type, cell_id_col):
        cases.append({"name": name, "kind": "sm", "format": fmt, "input": str(path.relative_to(out)),
                      "datasetType": dataset_type, "cellIdCol": cell_id_col})

    sc("h5ad_dense", "h5ad", write_h5ad(out / "sc/h5ad_dense/input.h5ad", seed))
    sc("h5ad_csr_3d", "h5ad", write_h5ad(out / "sc/h5ad_csr_3d/input.h5ad", seed + 1, sparse_x=True, dims=3))
    sc("h5ad_obs_coords", "h5ad", write_h5ad(out / "sc/h5ad_obs_coords/input.h5ad", seed + 2, coords_in_obs=True))
    sc("xenium", "xenium", write_xenium(out / "sc/xenium", seed + 3))
    sc("merscope", "merscope", write_merscope(out / "sc/merscope", seed + 4))
    sm("sm_parquet_xenium_3d", "parquet", write_sm_parquet(out / "sm/parquet_xenium_3d/input.parquet", seed + 5, preset="xenium", dims=3), "xenium", None)
    sm("sm_parquet_xenium_2d", "parquet", write_sm_parquet(out / "sm/parquet_xenium_2d/input.parquet", seed + 6, preset="xenium", dims=2), "xenium", None)
    sm("sm_parquet_merscope", "parquet", write_sm_parquet(out / "sm/parquet_merscope/input.parquet", seed + 7, preset="merscope", dims=3), "merscope", "cell_id")
    sm("sm_csv_merscope", "csv", write_sm_csv(out / "sm/csv_merscope/input.csv", seed + 8), "merscope", "cell_id")

    (out / "cases.json").write_text(json.dumps(cases, indent=1) + "\n")
    for c in cases:
        print(f"  {c['kind']}  {c['name']:<24} {c['input']}")


if __name__ == "__main__":
    main()
