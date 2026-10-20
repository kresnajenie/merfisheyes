"""scripts/prepare_visium.py on a Space Ranger-shaped Visium output made here
(a 10x h5, tissue_positions.csv, scalefactors, the two PNGs), judged from
outside: the validator passes, the spots are the cells' region with their
radius, X has its CSC copy, the images are in their own systems, and the
spatialdata library reads the store back.

  python -m pytest tests/python/test_prepare_visium.py -q

The Visium HD path needs a real output: set MERFISHEYES_VISIUM_HD to one.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

spatialdata = pytest.importorskip("spatialdata")
pytest.importorskip("spatialdata_io")
h5py = pytest.importorskip("h5py")
PIL = pytest.importorskip("PIL")

import prepare_visium  # noqa: E402

SPOTS, GENES = 80, 25


def make_spaceranger(folder: Path) -> Path:
    """A Space Ranger `outs/` for an 8 × 10 spot grid on a 2000 × 1600 px slide."""
    from PIL import Image

    rng = np.random.default_rng(5)
    (folder / "spatial").mkdir(parents=True)
    barcodes = [f"{i:05d}-1" for i in range(SPOTS)]
    genes = [f"Gene{i}" for i in range(GENES)]
    counts = rng.poisson(1.2, size=(SPOTS, GENES)).astype(np.int32)
    # 10x h5: genes × barcodes CSC, i.e. the barcodes × genes CSR arrays
    import scipy.sparse as sp

    csr = sp.csr_matrix(counts)
    with h5py.File(folder / "filtered_feature_bc_matrix.h5", "w") as f:
        f.attrs["library_ids"] = np.array([b"sample"])
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=np.array([b.encode() for b in barcodes]))
        m.create_dataset("data", data=csr.data)
        m.create_dataset("indices", data=csr.indices)
        m.create_dataset("indptr", data=csr.indptr)
        m.create_dataset("shape", data=np.array([GENES, SPOTS], dtype=np.int32))
        feat = m.create_group("features")
        feat.create_dataset("id", data=np.array([f"ENSG{i:05d}".encode() for i in range(GENES)]))
        feat.create_dataset("name", data=np.array([g.encode() for g in genes]))
        feat.create_dataset("feature_type", data=np.array([b"Gene Expression"] * GENES))
        feat.create_dataset("genome", data=np.array([b"GRCh38"] * GENES))
        feat.create_dataset("_all_tag_keys", data=np.array([b"genome"]))

    rows, cols = np.divmod(np.arange(SPOTS), 10)
    positions = "\n".join(
        f"{b},1,{r},{c},{100 + r * 180},{120 + c * 170}" for b, r, c in zip(barcodes, rows, cols)
    )
    (folder / "spatial" / "tissue_positions.csv").write_text(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n" + positions + "\n"
    )
    hires, lowres = 0.3, 0.075
    (folder / "spatial" / "scalefactors_json.json").write_text(
        json.dumps({"tissue_hires_scalef": hires, "tissue_lowres_scalef": lowres,
                    "fiducial_diameter_fullres": 150.0, "spot_diameter_fullres": 110.0})
    )
    for scale, name in ((hires, "tissue_hires_image.png"), (lowres, "tissue_lowres_image.png")):
        h, w = int(1600 * scale), int(2000 * scale)
        img = (rng.uniform(80, 220, size=(h, w, 3))).astype(np.uint8)
        Image.fromarray(img).save(folder / "spatial" / name)
    return folder


def fails(report):
    return [m for lv, _, m in report.lines if lv == "FAIL"]


def test_a_visium_output_becomes_a_store_with_spots_and_images(tmp_path):
    outs = make_spaceranger(tmp_path / "sample_outs")
    folder = tmp_path / "sample"
    r = prepare_visium.prepare_visium(outs, folder, log=lambda *_: None)
    assert not r.failed, fails(r)

    store = folder / "sample.zarr"
    nodes = json.loads((store / "zarr.json").read_text())["consolidated_metadata"]["metadata"]
    table = nodes["tables/table"]["attributes"]
    assert table["region"] == "sample" and table["instance_key"] == "spot_id"
    assert nodes["tables/table/X"]["attributes"]["encoding-type"] == "csr_matrix"
    assert nodes["tables/table/layers/X_csc"]["attributes"]["encoding-type"] == "csc_matrix"
    assert (store / "shapes" / "sample" / "shapes.parquet").is_file()
    images = sorted(k for k in nodes if k.startswith("images/") and k.count("/") == 1)
    assert images == ["images/sample_hires_image", "images/sample_lowres_image"]
    # The viewer's rule (store-layout.md §8): no global, the spots' system is the target
    assert any(w == "coordinates" and "region identity" in m for lv, w, m in r.lines if lv == "OK"), r.lines
    assert json.loads((folder / "mapping.json").read_text()) == {"cells": "sample.zarr"}

    # Read back by the library that wrote it
    sdata = spatialdata.read_zarr(store)
    spots = sdata.shapes["sample"]
    assert len(spots) == SPOTS and set(spots["radius"]) == {55.0}
    assert sdata.tables["table"].n_obs == SPOTS
    assert set(sdata.tables["table"].obs["spot_id"]) == set(spots.index)


@pytest.mark.skipif(not os.environ.get("MERFISHEYES_VISIUM_HD"), reason="MERFISHEYES_VISIUM_HD not set")
def test_a_visium_hd_output(tmp_path):
    folder = tmp_path / "hd"
    r = prepare_visium.prepare_visium(Path(os.environ["MERFISHEYES_VISIUM_HD"]), folder, name="hd", log=lambda *_: None)
    assert not r.failed, fails(r)
    store = folder / "hd.zarr"
    root = json.loads((store / "zarr.json").read_text())
    assert root["attributes"]["metadata_general"]["Default Table"] == "square_008um"
    nodes = root["consolidated_metadata"]["metadata"]
    assert nodes["tables/square_008um/layers/X_csc"]["attributes"]["encoding-type"] == "csc_matrix"
    assert "tables/square_002um" not in nodes
