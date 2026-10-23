#!/usr/bin/env python3
"""Build the tiny SpatialData store whose table carries its DE stats.

  python scripts/testdata/make-tiny-de.py

Outputs into tests/data/tiny/zarr/:
  spatialdata_v3_de.zarr   spatialdata_v3.zarr with scripts/de_stats.py run on
                           its table (varm/de_mean__celltype, de_pct__celltype,
                           uns/merfisheyes_de)
  de_expected.json         what the viewer must read back, computed here by
                           hand from expected.json, not by de_stats.py

Run after make-tiny-zarr.py (which clears the folder).
"""

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import zarr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import de_stats  # noqa: E402

ZARR = ROOT / "tests" / "data" / "tiny" / "zarr"
OUT = ZARR / "spatialdata_v3_de.zarr"

shutil.rmtree(OUT, ignore_errors=True)
shutil.copytree(ZARR / "spatialdata_v3.zarr", OUT)
assert de_stats.write_store_de_stats(OUT, log=print) == ["celltype"]
zarr.consolidate_metadata(OUT)

# By hand: normalize_total to the median total, log1p, mean / fraction per celltype
expected = json.loads((ZARR / "expected.json").read_text())
x = np.array(expected["expression"], dtype=np.float64).T  # cells × genes
celltype = np.array(expected["celltype"])
totals = x.sum(1)
target = np.median(totals[totals > 0])
lognorm = np.log1p(np.divide(x * target, totals[:, None], out=np.zeros_like(x), where=totals[:, None] > 0))
categories = sorted(set(celltype))
(ZARR / "de_expected.json").write_text(json.dumps({
    "celltypes": categories,
    "cellCounts": [int((celltype == c).sum()) for c in categories],
    # genes × categories, row-major
    "means": [float(lognorm[celltype == c, g].mean()) for g in range(x.shape[1]) for c in categories],
    "pctExpressing": [float((x[celltype == c, g] > 0).mean()) for g in range(x.shape[1]) for c in categories],
}))
print(f"wrote {OUT.name} and de_expected.json")
