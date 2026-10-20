#!/usr/bin/env python3
"""Build the committed tiny points fixture (docs/points-zarr-spec.md).

  python scripts/testdata/make-tiny-points.py

Outputs into tests/data/tiny/zarr/:
  points_v1.zarr          a store holding one points element,
                          points/transcripts: 3,000 points, 4 genes, 3 levels,
                          with chunks and shards small enough that reads cross
                          both
  points_expected.json    what it must read back as

Needs zarr >= 3 and numcodecs >= 0.16.
"""

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import zarr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import points_zarr as pz  # noqa: E402

OUT = ROOT / "tests" / "data" / "tiny" / "zarr"
store = OUT / "points_v1.zarr"
shutil.rmtree(store, ignore_errors=True)

pz.CHUNK, pz.SHARD = 128, 512
rng = np.random.default_rng(7)
n = 3000
genes = ["Alpha", "Beta", "Gamma", "Delta"]
x = rng.uniform(0, 1000, n).astype(np.float32)
y = rng.uniform(0, 600, n).astype(np.float32)
z = rng.uniform(-5, 5, n).astype(np.float32)
f = rng.choice(4, n, p=[0.6, 0.25, 0.1, 0.05]).astype(np.uint16)

root = zarr.open_group(store, mode="w")
root.attrs["spatialdata_attrs"] = {"version": "0.2"}
root.require_group("points")
w = pz.PointsWriter(store / "points" / "transcripts", OUT / "points_v1.tmp",
                    features=genes, dims=3, origin=(0, 0), extent=(1000, 600),
                    tile_size=100, band_um=200, n_levels=3)
w.spill({"x": x, "y": y, "z": z, "f": f}, seed=1, tag="a")
w.finish()
zarr.consolidate_metadata(str(store))

bbox = [250, 100, 520, 330]
inside = (f == 0) & (x >= bbox[0]) & (x < bbox[2]) & (y >= bbox[1]) & (y < bbox[3])
g = zarr.open_group(store / "points" / "transcripts", mode="r")
expected = {
    "genes": genes,
    "counts": [int((f == i).sum()) for i in range(4)],
    "levels": [e["points:n_points"] for e in g.attrs["multiscales"]["layout"]],
    "beta_xyz_sorted": sorted(
        [round(float(a), 3), round(float(b), 3), round(float(c), 3)]
        for a, b, c in zip(x[f == 1], y[f == 1], z[f == 1])
    ),
    "alpha_bbox": bbox,
    "alpha_in_bbox": int(inside.sum()),
}
(OUT / "points_expected.json").write_text(json.dumps(expected))
print(pz.validate(store / "points" / "transcripts"))
