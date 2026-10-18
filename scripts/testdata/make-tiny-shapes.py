#!/usr/bin/env python3
"""Build the committed tiny shapes fixture: cell polygons for the tiny
SpatialData table (tests/data/tiny/zarr/spatialdata_v3.zarr), as a
SpatialData shapes element whose GeoParquet the viewer reads as written.

  python scripts/testdata/make-tiny-shapes.py

Outputs into tests/data/tiny/zarr/shapes_v1/cell_boundaries/:
  zarr.json       the shapes group (ngff:shapes, identity transform)
  shapes.parquet  one Polygon (WKB) per cell, indexed by the table's
                  instance ids ("cell-0" .. "cell-29"), with the `geo`
                  metadata GeoParquet readers look for

The E2E fixture server merges this group into the tiny store's consolidated
metadata (e2e/helpers/spatialdata.ts), so a test can open "cells + shapes"
without a store that has both on disk. Needs only numpy and pyarrow: the
WKB is written by hand, no shapely.
"""

import json
import struct
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[2]
TINY = ROOT / "tests" / "data" / "tiny" / "zarr"
OUT = TINY / "shapes_v1" / "cell_boundaries"
CELLS = 30

# The same positions make-tiny-zarr.py wrote (same seed and draw order:
# the X matrix is drawn first, then the coordinates).
rng = np.random.default_rng(0)
rng.poisson(0.6, size=(CELLS, 8))
spatial = rng.uniform(0, 100, size=(CELLS, 2))


def wkb_polygon(points: np.ndarray) -> bytes:
    """A little-endian WKB Polygon with one closed ring."""
    ring = np.vstack([points, points[:1]])
    out = struct.pack("<BIII", 1, 3, 1, len(ring))
    for x, y in ring:
        out += struct.pack("<dd", float(x), float(y))
    return out


polygons = []
for i, (cx, cy) in enumerate(spatial):
    # A hexagon of radius 2.5 units around the cell centre, slightly rotated
    # per cell so neighbours do not line up
    angles = np.linspace(0, 2 * np.pi, 7)[:-1] + i * 0.3
    pts = np.column_stack([cx + 2.5 * np.cos(angles), cy + 2.5 * np.sin(angles)])
    polygons.append(wkb_polygon(pts))

ids = [f"cell-{i}" for i in range(CELLS)]
table = pa.table(
    {"geometry": pa.array(polygons, pa.binary()), "cell_id": pa.array(ids)},
)
geo = {
    "version": "1.0.0",
    "primary_column": "geometry",
    "columns": {
        "geometry": {"encoding": "WKB", "geometry_types": ["Polygon"]},
    },
}
pandas = {
    "index_columns": ["cell_id"],
    "column_indexes": [],
    "columns": [
        {"name": "geometry", "field_name": "geometry", "pandas_type": "bytes", "numpy_type": "object", "metadata": None},
        {"name": "cell_id", "field_name": "cell_id", "pandas_type": "unicode", "numpy_type": "object", "metadata": None},
    ],
    "pandas_version": "2.2.0",
}
table = table.replace_schema_metadata(
    {b"geo": json.dumps(geo).encode(), b"pandas": json.dumps(pandas).encode()},
)

OUT.mkdir(parents=True, exist_ok=True)
pq.write_table(table, OUT / "shapes.parquet", compression="zstd")
(OUT / "zarr.json").write_text(
    json.dumps(
        {
            "zarr_format": 3,
            "node_type": "group",
            "attributes": {
                "encoding-type": "ngff:shapes",
                "spatialdata_attrs": {"version": "0.2"},
                "coordinateTransformations": [
                    {"type": "identity", "input": {"name": "xy", "axes": [{"name": "x", "type": "space"}, {"name": "y", "type": "space"}]}, "output": {"name": "global", "axes": [{"name": "x", "type": "space"}, {"name": "y", "type": "space"}]}}
                ],
            },
        },
        indent=1,
    ),
)
print(f"wrote {CELLS} polygons to {OUT.relative_to(ROOT)}")
