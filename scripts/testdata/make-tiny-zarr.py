#!/usr/bin/env python3
"""Build the committed tiny SpatialData fixture for the table adapter tests.

  python scripts/testdata/make-tiny-zarr.py

Outputs into tests/data/tiny/zarr/:
  spatialdata_v3.zarr    SpatialData-style store: tables/table with CSR X
                         plus a CSC copy in layers/X_csc
  expected.json          the values the table must read back as
"""

import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr

OUT = Path(__file__).resolve().parents[2] / "tests" / "data" / "tiny" / "zarr"
CELLS, GENES = 30, 8

rng = np.random.default_rng(0)
dense = rng.poisson(0.6, size=(CELLS, GENES)).astype(np.float32)
dense[:, 3] = 0  # an all-zero gene
genes = [f"gene{i}" for i in range(GENES)]
celltype = [["astro", "neuron", "microglia"][i % 3] for i in range(CELLS)]
spatial = rng.uniform(0, 100, size=(CELLS, 2))


def make(x):
    obs = pd.DataFrame(
        {
            "cell_id": [f"cell-{i}" for i in range(CELLS)],
            "region": pd.Categorical(["cell_labels"] * CELLS),
            "celltype": pd.Categorical(celltype),
            "area": np.linspace(1.5, 9.5, CELLS),
        },
        index=[str(i) for i in range(CELLS)],
    )
    return ad.AnnData(
        X=x, obs=obs, var=pd.DataFrame(index=genes), obsm={"spatial": spatial}
    )


def write(name, adata, group):
    root = zarr.open_group(OUT / name, mode="w", zarr_format=3)
    root.attrs["spatialdata_attrs"] = {"version": "0.2"}
    ad.io.write_elem(root, group, adata)
    root[group].attrs.update(
        {
            "spatialdata-encoding-type": "ngff:regions_table",
            "region": "cell_labels",
            "region_key": "region",
            "instance_key": "cell_id",
        }
    )


shutil.rmtree(OUT, ignore_errors=True)
OUT.mkdir(parents=True)

sdata_table = make(sp.csr_matrix(dense))
sdata_table.layers["X_csc"] = sp.csc_matrix(dense)
write("spatialdata_v3.zarr", sdata_table, "tables/table")

(OUT / "expected.json").write_text(
    json.dumps(
        {
            "genes": genes,
            "expression": dense.T.tolist(),
            "celltype": celltype,
            "spatial": spatial.astype(np.float32).ravel().tolist(),
        }
    )
)
