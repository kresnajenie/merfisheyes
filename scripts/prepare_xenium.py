#!/usr/bin/env python3
"""Turn a Xenium Onboard Analysis export into a folder MERFISHEYES opens as it is.

  python scripts/prepare_xenium.py <xenium_outs> <folder> [--name NAME]
                                   [--labels] [--keep-z] [--workers N]
                                   [--cluster] [--url BASE]

Reads the export with spatialdata-io (`spatialdata_io.xenium`): the cell and
nucleus boundaries, the transcripts, the morphology image(s) and the cell
table, and writes

  <folder>/<name>.zarr/   the SpatialData store (Zarr v3): tables/table with a
                          CSC copy of X, shapes/cell_boundaries and
                          nucleus_boundaries (one GeoParquet each), images/
                          (1024² chunks, a pyramid), points/transcripts
  <folder>/molecules/     the transcripts as one file per gene (what the
                          viewer streams; the store's parquet points are not
                          read), flat on the plane of the cells unless --keep-z
  <folder>/mapping.json   names the store and the molecules
  <folder>/annotations/   with --cluster, Leiden clusters + DEGs when the
                          table has no categorical column

then runs scripts/prepare_folder.py's checks and scripts/validate_store.py.
The raster cell and nucleus labels are not written unless --labels: the
viewer does not draw them, and they are the largest part of the store. The
table then annotates `cell_boundaries` (the polygons) instead of the labels.

Needs: pip install -r scripts/requirements-prepare.txt spatialdata-io
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_folder  # noqa: E402
import spatialdata_points_to_sm  # noqa: E402

CHUNK = 1024


def read_xenium(path: Path, *, labels: bool, log=print):
    try:
        from spatialdata_io import xenium
    except ImportError:
        raise SystemExit("pip install spatialdata-io   (Python 3.12 or newer)")

    log(f"reading {path} with spatialdata-io ...")
    sdata = xenium(
        path,
        cells_as_circles=False,
        cells_labels=labels,
        nucleus_labels=labels,
        morphology_mip=False,  # the focus image is the one drawn; the MIP doubles the store
        # No pyramid unless asked: the viewer draws the coarsest level whole as the backdrop
        image_models_kwargs={"chunks": (1, CHUNK, CHUNK), "scale_factors": [2, 2, 2, 2, 2]},
        labels_models_kwargs={"chunks": (CHUNK, CHUNK), "scale_factors": [2, 2, 2, 2, 2]},
    )
    if not labels and "table" in sdata.tables and "cell_boundaries" in sdata.shapes:
        # The reader ties the table to the raster labels; with those left out,
        # the polygons are the cells' region (the same cell ids). The viewer
        # places the image relative to the region, so it must be an element
        # that exists.
        import pandas as pd

        table = sdata.tables["table"]
        table.obs["region"] = pd.Categorical(["cell_boundaries"] * table.n_obs)
        sdata.set_table_annotates_spatialelement(
            "table", region="cell_boundaries", region_key="region", instance_key="cell_id"
        )
    return sdata


def prepare_xenium(
    path: Path,
    folder: Path,
    *,
    name: str | None = None,
    labels: bool = False,
    keep_z: bool = False,
    workers: int = 1,
    cluster: bool = False,
    url: str | None = None,
    expanded: bool = False,
    log=print,
):
    path, folder = Path(path), Path(folder)
    name = name or path.name.removesuffix("_outs")
    folder.mkdir(parents=True, exist_ok=True)
    store = folder / f"{name}.zarr"
    if store.exists():
        raise SystemExit(f"{store} exists already; remove it or pass --name")

    sdata = read_xenium(path, labels=labels, log=log)
    log(f"  {sdata}")
    log(f"writing {store} ...")
    sdata.write(store)

    molecules = folder / "molecules"
    if "transcripts" in sdata.points:
        log(f"writing {molecules} from points/transcripts ({'3D' if keep_z else 'flat'}) ...")
        spatialdata_points_to_sm.convert(store, molecules, points="transcripts", workers=workers, drop_z=not keep_z)

    log("checking the folder ...")
    return prepare_folder.prepare(
        folder, cells=store, molecules=molecules if molecules.exists() else None, cluster=cluster, url=url,
        expanded=expanded, workers=workers, log=log,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("xenium", type=Path, help="the Xenium output folder (holds experiment.xenium)")
    ap.add_argument("folder", type=Path, help="the folder to write; it is what you drop, save or sync")
    ap.add_argument("--name", help="the store's name (default: the export folder's, without _outs)")
    ap.add_argument("--labels", action="store_true", help="also write the raster cell and nucleus labels")
    ap.add_argument("--keep-z", action="store_true", help="keep the transcripts' z (for a 3D view on their own)")
    ap.add_argument("--workers", type=int, default=4, help="parallel gene-file writers")
    ap.add_argument("--cluster", action="store_true", help="Leiden clusters + DEGs when the table has no cell types")
    ap.add_argument("--url", help="the folder's public URL, for absolute links in mapping.json")
    ap.add_argument("--expanded", action="store_true", help="keep the molecule folder beside the store instead of inside it")
    args = ap.parse_args(argv)
    r = prepare_xenium(
        args.xenium,
        args.folder,
        name=args.name,
        labels=args.labels,
        keep_z=args.keep_z,
        workers=args.workers,
        cluster=args.cluster,
        url=args.url,
        expanded=args.expanded,
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
