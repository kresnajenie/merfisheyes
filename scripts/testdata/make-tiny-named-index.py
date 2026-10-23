#!/usr/bin/env python3
"""Build the tiny SpatialData store whose table's indexes have names.

  python scripts/testdata/make-tiny-named-index.py

Outputs tests/data/tiny/zarr/spatialdata_v3_named_index.zarr:
spatialdata_v3.zarr with the table rewritten from an AnnData whose var index
is named `feature_name` and obs index `cell_name`, as many lab h5ad files
are. anndata then stores the indexes as var/feature_name and obs/cell_name
(the group's `_index` attribute names them), not var/_index. Derived from
the committed store, so it reads back as the same expected.json. Run after
make-tiny-zarr.py (which clears the folder).
"""

import shutil
from pathlib import Path

import anndata as ad
import zarr

ZARR = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "zarr"
SOURCE = ZARR / "spatialdata_v3.zarr"
OUT = ZARR / "spatialdata_v3_named_index.zarr"

shutil.rmtree(OUT, ignore_errors=True)
shutil.copytree(SOURCE, OUT)

tables = zarr.open_group(OUT / "tables", mode="r+")
attrs = dict(tables["table"].attrs)
table = ad.io.read_elem(tables["table"])
table.var.index.name = "feature_name"
table.obs.index.name = "cell_name"
del tables["table"]
ad.io.write_elem(tables, "table", table)
tables["table"].attrs.update(attrs)
zarr.consolidate_metadata(OUT)
