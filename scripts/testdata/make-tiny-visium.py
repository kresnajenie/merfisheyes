#!/usr/bin/env python3
"""Build the committed tiny Visium-shaped SpatialData store.

  python scripts/testdata/make-tiny-visium.py

Outputs tests/data/tiny/zarr/visium_v3.zarr, written by the `spatialdata`
library itself (0.8, as SpatialData-DB writes its stores), laid out like a
Visium store there (docs/plans/SDDB.md §1):

- no `global` system: three images in three systems, `tiny_full_image`
  identity into `tiny`, `tiny_hires_image` into `tiny_downscaled_hires`
  (and by a scale of 4 into `tiny`), `tiny_lowres_image` into
  `tiny_downscaled_lowres` (scale 16 into `tiny`);
- the spots as circles (`shapes/tiny`: Point + radius), identity into
  `tiny` and scaled into the two downscaled systems;
- a table annotating the spots, `obsm/spatial` int64 in full-image pixels,
  CSR `X` with a CSC copy in `layers/csc_counts`, and the curated
  `metadata_general["Default Table"]`.

The viewer must draw the full image at one pixel per unit, under the spots,
and not open the downscaled images.

Needs spatialdata 0.8 (and its geopandas, shapely, anndata).
"""

import shutil
from pathlib import Path

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import scipy.sparse as sp
import spatialdata as sd
import zarr
from shapely.geometry import Point
from spatialdata.models import Image2DModel, ShapesModel, TableModel
from spatialdata.transformations import Identity, Scale

OUT = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "zarr" / "visium_v3.zarr"
SAMPLE = "tiny"
FULL, HIRES, LOWRES = 256, 64, 16  # image sides in pixels
SPOTS, GENES = 36, 6

rng = np.random.default_rng(0)


def image(side: int, rgb: tuple[int, int, int]) -> np.ndarray:
    """A flat colour with a darker diagonal, so placement can be seen."""
    a = np.empty((3, side, side), dtype=np.uint8)
    for c in range(3):
        a[c] = rgb[c]
    for i in range(side):
        a[:, i, max(0, i - 1) : i + 2] = 30
    return a


hires_scale = FULL / HIRES  # 4
lowres_scale = FULL / LOWRES  # 16
images = {
    f"{SAMPLE}_full_image": Image2DModel.parse(
        image(FULL, (200, 120, 160)),
        dims=("c", "y", "x"),
        scale_factors=[2],
        chunks=(3, 64, 64),
        transformations={SAMPLE: Identity()},
    ),
    f"{SAMPLE}_hires_image": Image2DModel.parse(
        image(HIRES, (120, 200, 160)),
        dims=("c", "y", "x"),
        transformations={
            f"{SAMPLE}_downscaled_hires": Identity(),
            SAMPLE: Scale([hires_scale, hires_scale], axes=("y", "x")),
        },
    ),
    f"{SAMPLE}_lowres_image": Image2DModel.parse(
        image(LOWRES, (120, 160, 200)),
        dims=("c", "y", "x"),
        transformations={
            f"{SAMPLE}_downscaled_lowres": Identity(),
            SAMPLE: Scale([lowres_scale, lowres_scale], axes=("y", "x")),
        },
    ),
}

# A 6 × 6 grid of spots across the full image
grid = np.arange(6) * 40 + 28
xy = np.array([(x, y) for y in grid for x in grid], dtype=np.int64)
spots = gpd.GeoDataFrame(
    {"radius": np.full(SPOTS, 12.0)},
    geometry=[Point(x, y) for x, y in xy],
    index=pd.Index([f"spot{i}" for i in range(SPOTS)], name="spot_id"),
)
shapes = {
    SAMPLE: ShapesModel.parse(
        spots,
        transformations={
            SAMPLE: Identity(),
            f"{SAMPLE}_downscaled_hires": Scale([1 / hires_scale] * 2, axes=("x", "y")),
            f"{SAMPLE}_downscaled_lowres": Scale([1 / lowres_scale] * 2, axes=("x", "y")),
        },
    )
}

counts = sp.csr_matrix(rng.poisson(1.0, size=(SPOTS, GENES)).astype(np.float32))
adata = ad.AnnData(
    X=counts,
    obs=pd.DataFrame(
        {
            "spot_id": spots.index.to_numpy(),
            "region": pd.Categorical([SAMPLE] * SPOTS),
            "cluster": pd.Categorical([["a", "b", "c"][i % 3] for i in range(SPOTS)]),
        },
        index=[str(i) for i in range(SPOTS)],
    ),
    var=pd.DataFrame(index=[f"gene{i}" for i in range(GENES)]),
    obsm={"spatial": xy},
    layers={"csc_counts": counts.tocsc()},
)
table = TableModel.parse(adata, region=SAMPLE, region_key="region", instance_key="spot_id")

sdata = sd.SpatialData(images=images, shapes=shapes, tables={"table": table})
if OUT.exists():
    shutil.rmtree(OUT)
sdata.write(OUT)

# The curated metadata SpatialData-DB adds, then the metadata consolidated
root = zarr.open_group(str(OUT), mode="r+")
root.attrs["metadata_general"] = {"Assay": "Visium", "Default Table": "table"}
zarr.consolidate_metadata(str(OUT))
print(f"wrote {OUT}")
