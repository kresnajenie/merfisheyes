#!/usr/bin/env python3
"""Build the committed tiny SpatialData fixture for the table adapter tests.

  python scripts/testdata/make-tiny-zarr.py

Outputs into tests/data/tiny/zarr/:
  spatialdata_v3.zarr          SpatialData-style store: tables/table with CSR
                               X plus a CSC copy in layers/X_csc, and a
                               multiscale image in images/dapi
  spatialdata_v3_sharded.zarr  the same store with the table as current
                               anndata writes it by default: sharded arrays
                               and nullable strings (plus a `note` column
                               with missing values)
  expected.json                the values the table must read back as

Both stores carry consolidated metadata, which is how a store at a public URL
is listed (the E2E suite opens them that way; see e2e/helpers/spatialdata.ts).

Needs anndata >= 0.12 (for the sharded store).
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


note = [None if i % 5 == 0 else f"note-{i % 3}" for i in range(CELLS)]


def make(x, nullable=False):
    """`nullable` writes strings the way pandas' string dtype is stored: a
    values array and a mask. Otherwise they are plain string arrays."""
    strings = "string" if nullable else object
    obs = pd.DataFrame(
        {
            "cell_id": pd.Series([f"cell-{i}" for i in range(CELLS)], dtype=strings),
            "region": pd.Categorical(["cell_labels"] * CELLS),
            "celltype": pd.Categorical(celltype),
            "area": np.linspace(1.5, 9.5, CELLS),
        },
    ).set_index(pd.Index([str(i) for i in range(CELLS)], dtype=strings))
    if nullable:
        obs["note"] = pd.array(note, dtype="string")
    return ad.AnnData(
        X=x,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, dtype=strings)),
        obsm={"spatial": spatial},
    )


IMAGE_SIZE, IMAGE_LEVELS = 64, 2
# A diagonal ramp with a bright square, so every tile has something to show
yy, xx = np.mgrid[0:IMAGE_SIZE, 0:IMAGE_SIZE]
image = ((yy + xx) * 400).astype(np.uint16)
image[16:32, 16:32] = 60000


def write_image(root, name):
    """An OME-Zarr image (NGFF 0.5, axes c/y/x) over the cells' 0-100 extent."""
    group = root.require_group("images").create_group(name)
    datasets = []
    for level in range(IMAGE_LEVELS):
        step = 2**level
        pixels = image.reshape(
            IMAGE_SIZE // step, step, IMAGE_SIZE // step, step
        ).mean(axis=(1, 3))
        group.create_array(
            str(level),
            data=pixels.astype(np.uint16)[None],
            chunks=(1, 32, 32),
        )
        size = 100 / IMAGE_SIZE * step
        datasets.append(
            {
                "path": str(level),
                "coordinateTransformations": [
                    {"type": "scale", "scale": [1.0, size, size]},
                    # Centre of pixel (0, 0)
                    {"type": "translation", "translation": [0.0, size / 2, size / 2]},
                ],
            }
        )
    group.attrs["ome"] = {
        "version": "0.5",
        "multiscales": [
            {
                "name": name,
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


def write(name, adata, group, shard=False):
    root = zarr.open_group(OUT / name, mode="w", zarr_format=3)
    root.attrs["spatialdata_attrs"] = {"version": "0.2"}
    parent, _, key = group.rpartition("/")
    with ad.settings.override(auto_shard_zarr_v3=shard):
        ad.io.write_elem(root.require_group(parent), key, adata)
    root[group].attrs.update(
        {
            "spatialdata-encoding-type": "ngff:regions_table",
            "region": "cell_labels",
            "region_key": "region",
            "instance_key": "cell_id",
        }
    )
    write_image(root, "dapi")
    zarr.consolidate_metadata(OUT / name)


shutil.rmtree(OUT, ignore_errors=True)
OUT.mkdir(parents=True)

sdata_table = make(sp.csr_matrix(dense))
sdata_table.layers["X_csc"] = sp.csc_matrix(dense)
write("spatialdata_v3.zarr", sdata_table, "tables/table")

sharded_table = make(sp.csr_matrix(dense), nullable=True)
sharded_table.layers["X_csc"] = sp.csc_matrix(dense)
write("spatialdata_v3_sharded.zarr", sharded_table, "tables/table", shard=True)

(OUT / "expected.json").write_text(
    json.dumps(
        {
            "genes": genes,
            "expression": dense.T.tolist(),
            "celltype": celltype,
            "note": ["" if n is None else n for n in note],
            "spatial": spatial.astype(np.float32).ravel().tolist(),
        }
    )
)
