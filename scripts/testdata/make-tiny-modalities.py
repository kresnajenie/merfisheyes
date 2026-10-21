#!/usr/bin/env python3
"""Build the committed tiny stores shaped like SpatialData-DB's non-RNA and
CosMx stores, written by the `spatialdata` library itself (0.8).

  python scripts/testdata/make-tiny-modalities.py

Outputs under tests/data/tiny/zarr/, each laid out as the real store is
(roots measured 2026-10-21, docs/plans/SDDB.md):

- codex_v3.zarr — a CODEX / PhenoCycler store: a dense `X` (every image
  channel a feature, `var/marker_role` telling markers from blanks), a
  6-channel uint16 image chunked one channel at a time, `cell_boundaries`
  polygons and the `cell_boundaries_for_vitessce` points (centre + radius),
  everything identity into `global`;
- maldi_v3.zarr — a MALDI-MSI store: the feature index column is
  `formula_adduct` (not `_index`), `X` CSR with a `csc_counts` copy,
  `maldi_spots` circles, and an RGB optical image placed into `global` by an
  affine that scales and rotates it;
- cosmx_v3.zarr — a CosMx store: the table has no `obsm/spatial`; the cells'
  positions are in `cells_polygons` (polygons, identity into `global`) and
  `cells_polygons_for_vitessce` (centre + radius + `global_cell_id`); no image.

Needs spatialdata 0.8 (and its geopandas, shapely, anndata).
"""

import math
import shutil
from pathlib import Path

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.sparse as sp
import spatialdata as sd
import zarr
from shapely.geometry import Point, Polygon
from spatialdata.models import Image2DModel, ShapesModel, TableModel
from spatialdata.transformations import Affine, Identity

OUT_DIR = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "zarr"
rng = np.random.default_rng(1)


def write(sdata: sd.SpatialData, name: str, assay: str, default_table: str = "table") -> Path:
    out = OUT_DIR / name
    if out.exists():
        shutil.rmtree(out)
    sdata.write(out)
    root = zarr.open_group(str(out), mode="r+")
    root.attrs["metadata_general"] = {"Assay": assay, "Default Table": default_table}
    zarr.consolidate_metadata(str(out))
    print(f"wrote {out}")
    return out


def hexagon(x: float, y: float, r: float) -> Polygon:
    return Polygon([(x + r * math.cos(a), y + r * math.sin(a)) for a in np.linspace(0, 2 * math.pi, 7)[:-1]])


def cells_on_grid(n_side: int, pitch: float, radius: float):
    centres = [(pitch * (i + 0.5), pitch * (j + 0.5)) for j in range(n_side) for i in range(n_side)]
    return centres, [hexagon(x, y, radius) for x, y in centres]


# ── CODEX: dense X, many image channels ─────────────────────────────────
SIDE, CHANNELS, CELLS_SIDE = 128, 6, 6
centres, polys = cells_on_grid(CELLS_SIDE, SIDE / CELLS_SIDE, 7.0)
n = len(centres)
ids = [f"cell{i}" for i in range(n)]
image = np.zeros((CHANNELS, SIDE, SIDE), dtype=np.uint16)
yy, xx = np.mgrid[0:SIDE, 0:SIDE]
for c in range(CHANNELS):
    # Each channel a different stripe pattern, so a channel switch can be seen
    image[c] = (2000 * ((xx // (8 * (c + 1)) + yy // 16) % 2) + 500 * c).astype(np.uint16)
markers = ["DAPI", "CD3", "CD20", "blank", "PanCK", "blank (2)"]
roles = ["nuclear", "marker", "marker", "blank", "marker", "blank"]
dense = rng.gamma(2.0, 100.0, size=(n, CHANNELS)).astype(np.float32)
codex_table = TableModel.parse(
    ad.AnnData(
        X=dense,
        obs=pd.DataFrame(
            {"region": pd.Categorical(["cell_boundaries"] * n), "cell_id": ids, "area": rng.uniform(50, 150, n)},
            index=ids,
        ),
        var=pd.DataFrame({"marker_name": markers, "marker_role": roles}, index=markers),
        obsm={"spatial": np.array(centres, dtype=np.float64)},
    ),
    region="cell_boundaries",
    region_key="region",
    instance_key="cell_id",
)
codex = sd.SpatialData(
    images={
        "reg001_image": Image2DModel.parse(
            image,
            dims=("c", "y", "x"),
            scale_factors=[2],
            chunks=(1, 64, 64),
            c_coords=markers,
            transformations={"global": Identity()},
        )
    },
    shapes={
        "cell_boundaries": ShapesModel.parse(
            gpd.GeoDataFrame(geometry=polys, index=pd.Index(ids)),
            transformations={"global": Identity()},
        ),
        "cell_boundaries_for_vitessce": ShapesModel.parse(
            gpd.GeoDataFrame({"radius": np.full(n, 7.0)}, geometry=[Point(x, y) for x, y in centres], index=pd.Index(ids)),
            transformations={"global": Identity()},
        ),
    },
    tables={"table": codex_table},
)
write(codex, "codex_v3.zarr", "CODEX")

# ── MALDI-MSI: a named feature index and an affine image ────────────────
SPOTS_SIDE, PITCH, IONS = 8, 20.0, 5
spots = [(PITCH * (i + 0.5), PITCH * (j + 0.5)) for j in range(SPOTS_SIDE) for i in range(SPOTS_SIDE)]
m = len(spots)
ions = ["C6H12O6+H", "C6H12O6+Na", "C18H34O2-H", "C5H9NO4+K", "C27H46O+H"]
counts = sp.csr_matrix(rng.poisson(3.0, size=(m, IONS)).astype(np.float32))
maldi_table = TableModel.parse(
    ad.AnnData(
        X=counts,
        obs=pd.DataFrame(
            {
                "instance_id": np.arange(m),
                "region": pd.Categorical(["maldi_spots"] * m),
                "total_ion_current": rng.uniform(1e4, 5e4, m),
            },
            index=[str(i) for i in range(m)],
        ),
        var=pd.DataFrame(
            {"formula": [i.split("+")[0].split("-")[0] for i in ions], "mz": rng.uniform(100, 900, IONS)},
            index=pd.Index(ions, name="formula_adduct"),
        ),
        obsm={"spatial": np.array(spots, dtype=np.int64)},
        layers={"csc_counts": counts.tocsc()},
    ),
    region="maldi_spots",
    region_key="region",
    instance_key="instance_id",
)
# The optical image: 40 × 32 pixels, each pixel 4 units, rotated 30° and
# shifted so that it covers the spots' grid (0–160 × 0–160)
OPT_W, OPT_H, OPT_SCALE, ANGLE = 40, 32, 4.0, math.radians(30)
optical = np.zeros((3, OPT_H, OPT_W), dtype=np.uint8)
optical[0] = 200
optical[1] = (np.arange(OPT_W) * 6).astype(np.uint8)[None, :]
optical[2] = (np.arange(OPT_H) * 8).astype(np.uint8)[:, None]
cos_a, sin_a = math.cos(ANGLE), math.sin(ANGLE)
# Homogeneous: rows = output axes (c, y, x) + 1, columns = input axes
# (c, y, x) + translation; the library drops the last row when it writes
affine = np.array(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, OPT_SCALE * cos_a, OPT_SCALE * sin_a, -20.0],
        [0.0, -OPT_SCALE * sin_a, OPT_SCALE * cos_a, 40.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
)
maldi = sd.SpatialData(
    images={
        "optical_image": Image2DModel.parse(
            optical,
            dims=("c", "y", "x"),
            c_coords=["r", "g", "b"],
            transformations={"global": Affine(affine, input_axes=("c", "y", "x"), output_axes=("c", "y", "x"))},
        )
    },
    shapes={
        "maldi_spots": ShapesModel.parse(
            gpd.GeoDataFrame({"radius": np.full(m, PITCH / 2)}, geometry=[Point(x, y) for x, y in spots]),
            transformations={"global": Identity()},
        )
    },
    tables={"table": maldi_table},
)
write(maldi, "maldi_v3.zarr", "MALDI-MSI")

# ── CosMx: no obsm/spatial ──────────────────────────────────────────────
centres, polys = cells_on_grid(5, 30.0, 10.0)
k = len(centres)
gids = [f"c_1_{i}" for i in range(k)]
proteins = ["PanCK", "CD45", "CD68", "Membrane"]
cosmx_counts = sp.csr_matrix(rng.gamma(2.0, 50.0, size=(k, len(proteins))).astype(np.float32))
cosmx_table = TableModel.parse(
    ad.AnnData(
        X=cosmx_counts,
        obs=pd.DataFrame(
            {
                "fov": np.ones(k, dtype=np.int64),
                "Area": rng.uniform(100, 400, k),
                "global_cell_id": gids,
                "region_key": pd.Categorical(["cells_polygons"] * k),
            },
            index=gids,
        ),
        var=pd.DataFrame({"marker_name": proteins}, index=proteins),
        layers={"csc_counts": cosmx_counts.tocsc()},
    ),
    region="cells_polygons",
    region_key="region_key",
    instance_key="global_cell_id",
)
cosmx = sd.SpatialData(
    shapes={
        "cells_polygons": ShapesModel.parse(
            gpd.GeoDataFrame({"fov": np.ones(k, dtype=np.int64)}, geometry=polys, index=pd.Index(gids)),
            transformations={"global": Identity()},
        ),
        # As SpatialData-DB writes it: the cell id is the frame's index
        "cells_polygons_for_vitessce": ShapesModel.parse(
            gpd.GeoDataFrame(
                {"radius": np.full(k, 10.0)},
                geometry=[Point(x, y) for x, y in centres],
                index=pd.Index(gids, name="global_cell_id"),
            ),
            transformations={"global": Identity()},
        ),
    },
    tables={"table": cosmx_table},
)
# As SpatialData-DB writes it: no obsm at all
out = write(cosmx, "cosmx_v3.zarr", "CosMx")
table = zarr.open_group(str(out / "tables" / "table"), mode="r+")
if "obsm" in table:
    del table["obsm"]
zarr.consolidate_metadata(str(out))
print("cosmx table obsm removed")
