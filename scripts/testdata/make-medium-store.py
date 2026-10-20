#!/usr/bin/env python3
"""Build the generated medium SpatialData store the budget suite measures on
(e2e/perf/elements.spec.ts): one section's worth of cells and molecules,
small enough to build in CI in about a minute, large enough that its
heaviest gene (60 % of the molecules) exceeds the Auto point budget, so a
whole-section view is drawn from a coarser level and a zoom reads the finer.

  python scripts/testdata/make-medium-store.py [--molecules 30000000] [--cells 50000]

Outputs into tests/data/generated/medium/points-store/ (gitignored):
  store.zarr/tables/table        cells: CSC X, one categorical column, obsm/spatial
  store.zarr/points/transcripts  molecules in the points format (docs/points-zarr-spec.md)
  store.zarr/images/dapi         a small multiscale image under the cells
  expected.json                  counts to assert on

Deterministic (seeded), so the same store comes out on every machine.
Needs anndata >= 0.12, zarr >= 3, numcodecs >= 0.16.
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import points_zarr as pz  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--molecules", type=int, default=30_000_000)
parser.add_argument("--cells", type=int, default=50_000)
parser.add_argument("--genes", type=int, default=500)
args = parser.parse_args()

OUT = ROOT / "tests" / "data" / "generated" / "medium" / "points-store"
STORE = OUT / "store.zarr"
EXTENT = (8000.0, 6000.0)  # µm, about one MERSCOPE section

shutil.rmtree(OUT, ignore_errors=True)
OUT.mkdir(parents=True)
rng = np.random.default_rng(42)
# Gene names: the three the viewer picks by default (lib/utils/auto-select-genes.ts)
# sit at modest ranks, so opening the store is light; Gene000 is most of the
# molecules, heavier than a third of the Auto budget, so selecting it makes
# the point budget bite and a zoom read new levels
genes = [f"Gene{i:03d}" for i in range(args.genes)]
for name, rank in (("Igfbpl1", 40), ("Drd1", 80), ("Th", 160)):
    genes[rank] = name

# ── cells ────────────────────────────────────────────────────────────────
# Clustered positions: a few hundred blobs, like tissue, so a zoomed view
# is dense and a zoomed-out one is not uniform
centres = rng.uniform([0, 0], EXTENT, size=(300, 2))
blob = rng.integers(0, len(centres), args.cells)
spatial = centres[blob] + rng.normal(0, 120, size=(args.cells, 2))
spatial = np.clip(spatial, [0, 0], EXTENT).astype(np.float32)
# Sparse expression: ~5 % of entries, Poisson counts
nnz_per_cell = max(1, args.genes // 20)
rows = np.repeat(np.arange(args.cells), nnz_per_cell)
cols = rng.integers(0, args.genes, args.cells * nnz_per_cell)
vals = rng.poisson(2.0, args.cells * nnz_per_cell).astype(np.float32) + 1
X = sp.csc_matrix((vals, (rows, cols)), shape=(args.cells, args.genes))
X.sum_duplicates()
leiden = pd.Categorical([str(b % 12) for b in blob])
obs = pd.DataFrame(
    {
        "cell_id": [f"cell-{i}" for i in range(args.cells)],
        "leiden": leiden,
        "region": pd.Categorical(["cell_labels"] * args.cells),
    },
    index=pd.Index([str(i) for i in range(args.cells)]),
)
adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
adata.obsm["spatial"] = spatial.astype(np.float64)

root = zarr.open_group(STORE, mode="w", zarr_format=3)
root.attrs["spatialdata_attrs"] = {"version": "0.2"}
# Nested: the format must be 3 before autosharding turns on, and
# autosharding must be off again before the format is restored
with ad.settings.override(zarr_write_format=3):
    with ad.settings.override(auto_shard_zarr_v3=True):
        ad.io.write_elem(root.require_group("tables"), "table", adata)
root["tables/table"].attrs.update(
    {
        "spatialdata-encoding-type": "ngff:regions_table",
        "region": "cell_labels",
        "region_key": "region",
        "instance_key": "cell_id",
    }
)

# ── image ────────────────────────────────────────────────────────────────
# A soft blob image under the cells, 2 levels, uint16, 1024² chunks
img = root.require_group("images").require_group("dapi")
h, w = 1500, 2000
yy, xx = np.mgrid[0:h, 0:w]
base = np.zeros((h, w), np.float32)
for cx, cy in centres[:60]:
    base += np.exp(-(((xx * 4 - cx) ** 2 + (yy * 4 - cy) ** 2) / (2 * 150**2)))
base = (base / base.max() * 4000 + rng.normal(0, 30, base.shape)).clip(0, 65535)
datasets = []
for level, factor in enumerate([1, 2]):
    arr = base[::factor, ::factor].astype(np.uint16)[None]
    a = img.create_array(
        str(level), shape=arr.shape, chunks=(1, 1024, 1024), dtype="uint16",
    )
    a[:] = arr
    datasets.append(
        {
            "path": str(level),
            "coordinateTransformations": [
                {"type": "scale", "scale": [1.0, 4.0 * factor, 4.0 * factor]},
                {"type": "translation", "translation": [0.0, 0.0, 0.0]},
            ],
        }
    )
img.attrs["ome"] = {
    "version": "0.5",
    "multiscales": [
        {
            "name": "dapi",
            "axes": [
                {"name": "c", "type": "channel"},
                {"name": "y", "type": "space", "unit": "micrometer"},
                {"name": "x", "type": "space", "unit": "micrometer"},
            ],
            "datasets": datasets,
        }
    ],
    "omero": {"channels": [{"label": "DAPI"}]},
}

# ── molecules ────────────────────────────────────────────────────────────
# Each molecule sits near a cell; gene frequencies follow a power law so a
# few genes are abundant (the ones a view budget bites on)
root.require_group("points")
w = pz.PointsWriter(
    STORE / "points" / "transcripts",
    OUT / "points.tmp",
    features=genes,
    dims=2,
    origin=(0.0, 0.0),
    extent=EXTENT,
    tile_size=100,
    band_um=800,
    n_levels=6,
)
freq = 1.0 / np.arange(1, args.genes + 1)
freq[0] = 0  # set below
freq = freq / freq.sum() * 0.4
freq[0] = 0.6  # Gene000: 60 % of all molecules
per_spill = 2_000_000
remaining = args.molecules
part = 0
while remaining > 0:
    n = min(per_spill, remaining)
    cell = rng.integers(0, args.cells, n)
    xy = spatial[cell] + rng.normal(0, 8, size=(n, 2))
    xy = np.clip(xy, [0, 0], [EXTENT[0] - 1e-3, EXTENT[1] - 1e-3])
    w.spill(
        {
            "x": xy[:, 0].astype(np.float32),
            "y": xy[:, 1].astype(np.float32),
            "f": rng.choice(args.genes, n, p=freq).astype(np.uint16),
        },
        seed=part,
        tag=f"p{part}",
    )
    remaining -= n
    part += 1
w.finish()
zarr.consolidate_metadata(str(STORE))

(OUT / "expected.json").write_text(
    json.dumps(
        {
            "cells": args.cells,
            "genes": args.genes,
            "molecules": args.molecules,
            "heavy_gene": "Gene000",
            "default_genes": ["Igfbpl1", "Drd1", "Th"],
        }
    )
)
try:
    print(pz.validate(STORE / "points" / "transcripts"))
except Exception as e:  # noqa: BLE001 — the validator, not the store (see docs/specs/AUDIT-2026-10-17.md)
    print(f"validate skipped: {type(e).__name__}: {e}")
print(f"wrote {args.cells:,} cells and {args.molecules:,} molecules to {OUT.relative_to(ROOT)}")
