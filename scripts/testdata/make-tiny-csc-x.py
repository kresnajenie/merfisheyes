#!/usr/bin/env python3
"""Build the tiny SpatialData store whose table keeps X itself as CSC.

  python scripts/testdata/make-tiny-csc-x.py

Outputs tests/data/tiny/zarr/spatialdata_v3_csc_x.zarr: spatialdata_v3.zarr
with tables/table/X rewritten as a csc_matrix and no layers/X_csc, the shape
scripts/prepare_merscope.py writes. Derived from the committed store, so it
reads back as the same expected.json. Run after make-tiny-zarr.py (which
clears the folder).
"""

import shutil
from pathlib import Path

import anndata as ad
import scipy.sparse as sp
import zarr

ZARR = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "zarr"
SOURCE, OUT = ZARR / "spatialdata_v3.zarr", ZARR / "spatialdata_v3_csc_x.zarr"

shutil.rmtree(OUT, ignore_errors=True)
shutil.copytree(SOURCE, OUT)

table = zarr.open_group(OUT / "tables" / "table", mode="r+")
x = sp.csc_matrix(ad.io.read_elem(table["X"]))
del table["X"]
del table["layers"]["X_csc"]
ad.io.write_elem(table, "X", x)
zarr.consolidate_metadata(OUT)
