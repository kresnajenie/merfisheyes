#!/usr/bin/env python3
"""Turn a Xenium Onboard Analysis export into a folder MERFISHEYES opens as it is.

  python scripts/prepare_xenium.py <xenium_outs> <folder> [--name NAME]
                                   [--labels] [--workers N]
                                   [--cluster] [--url BASE]

Reads the export with spatialdata-io (`spatialdata_io.xenium`): the cell and
nucleus boundaries, the morphology image (when the export has it) and the
cell table, and writes

  <folder>/<name>.zarr/   the SpatialData store (Zarr v3): tables/table with a
                          CSC copy of X, shapes/cell_boundaries and
                          nucleus_boundaries (one GeoParquet each), images/
                          (1024² chunks, a pyramid; when present),
                          points/molecules (tiled, straight from
                          transcripts.parquet, flat on the cells' plane)
  <folder>/mapping.json   names the store
  <folder>/annotations/   with --cluster, Leiden clusters + DEGs when the
                          table has no categorical column

With --expanded the molecules are instead a per-gene folder beside the store.
An export without its morphology files (an upload without images) works the
same; it needs experiment.xenium, cell_feature_matrix.h5, cells.parquet,
cells.zarr.zip, cell_boundaries.parquet, nucleus_boundaries.parquet and
transcripts.parquet. scripts/prepare.py runs this for a Xenium export.

then runs scripts/prepare_folder.py's checks and scripts/validate_store.py.
The raster cell and nucleus labels are not written unless --labels: the
viewer does not draw them, and they are the largest part of the store. The
table then annotates `cell_boundaries` (the polygons) instead of the labels.

Needs: pip install -r scripts/requirements.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_folder  # noqa: E402

CHUNK = 1024


def has_images(path: Path) -> bool:
    """Whether the export holds its morphology image (an upload may not:
    the server sends images only when asked)."""
    return any((path / n).exists() for n in ("morphology_focus", "morphology_focus.ome.tif", "morphology.ome.tif"))


GRAPHCLUST = "analysis/clustering/gene_expression_graphclust/clusters.csv"


def graph_clusters(path: Path):
    """10x's own clustering of the cells — Xenium Onboard Analysis'
    graph-based clusters, `analysis/clustering/gene_expression_graphclust/
    clusters.csv` (Barcode, Cluster) — from analysis.tar.gz or an unpacked
    analysis/ folder, as (column name, labels by cell id); None when the
    export has neither. The browser preview shows the same column
    (XeniumAdapter._augmentClustersFromAnalysis), so it must survive a save."""
    import io
    import tarfile

    import pandas as pd

    unpacked = path / GRAPHCLUST
    if unpacked.exists():
        frame = pd.read_csv(unpacked, dtype=str)
    elif (path / "analysis.tar.gz").exists():
        with tarfile.open(path / "analysis.tar.gz") as archive:
            member = next((m for m in archive.getmembers() if m.name.endswith(GRAPHCLUST)), None)
            if member is None:
                return None
            frame = pd.read_csv(io.BytesIO(archive.extractfile(member).read()), dtype=str)
    else:
        return None
    if frame.shape[1] < 2:
        return None
    ids, column = frame.columns[0], frame.columns[1]
    return column, frame.set_index(frame[ids].str.strip())[column].str.strip()


def read_xenium(path: Path, *, labels: bool, images: bool = True, log=print):
    try:
        from spatialdata_io import xenium
    except ImportError:
        raise SystemExit("pip install spatialdata-io   (Python 3.12 or newer)")

    log(f"reading {path} with spatialdata-io ...")
    images = images and has_images(path)
    sdata = xenium(
        path,
        cells_as_circles=False,
        cells_labels=labels,
        nucleus_labels=labels,
        # The molecules go into the tiled points element from transcripts.parquet;
        # the reader's own parquet copy would only double them
        transcripts=False,
        morphology_focus=images,
        aligned_images=images,
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
    clusters = graph_clusters(path)
    if clusters is not None and "table" in sdata.tables:
        import pandas as pd

        column, labels_by_id = clusters
        table = sdata.tables["table"]
        # A cell the clustering left out gets "", as in the preview
        found = labels_by_id.reindex(table.obs["cell_id"].astype(str).to_numpy())
        table.obs[column] = pd.Categorical(found.fillna("").to_numpy())
        log(f"  {GRAPHCLUST}: {column}, {found.notna().sum():,} of {table.n_obs:,} cells")
    return sdata


def prepare_xenium(
    path: Path,
    folder: Path,
    *,
    name: str | None = None,
    labels: bool = False,
    workers: int = 1,
    cluster: bool = False,
    url: str | None = None,
    expanded: bool = False,
    images: bool = True,
    label_csv: Path | list[Path] | None = None,
    keep_z: bool = False,
    log=print,
):
    path, folder = Path(path), Path(folder)
    name = name or path.name.removesuffix("_outs")
    folder.mkdir(parents=True, exist_ok=True)
    store = folder / f"{name}.zarr"
    if store.exists():
        raise SystemExit(f"{store} exists already; remove it or pass --name")

    sdata = read_xenium(path, labels=labels, images=images, log=log)
    log(f"  {sdata}")
    log(f"writing {store} ...")
    sdata.write(store)

    # The transcripts go straight from the export's table into the store's
    # tiled points element (or, --expanded, a per-gene folder beside it)
    transcripts = path / "transcripts.parquet"
    if not transcripts.exists():
        transcripts = next(path.glob("transcripts.csv*"), None)
    log("checking the folder ...")
    return prepare_folder.prepare(
        folder, cells=store, molecules=transcripts, molecule_columns="xenium", cluster=cluster, url=url,
        expanded=expanded, workers=workers, labels=label_csv, keep_z=keep_z, log=log,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("xenium", type=Path, help="the Xenium output folder (holds experiment.xenium)")
    ap.add_argument("folder", type=Path, help="the folder to write; it is what you drop, save or sync")
    ap.add_argument("--name", help="the store's name (default: the export folder's, without _outs)")
    ap.add_argument("--labels", action="store_true", help="also write the raster cell and nucleus labels")
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
        workers=args.workers,
        cluster=args.cluster,
        url=args.url,
        expanded=args.expanded,
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
