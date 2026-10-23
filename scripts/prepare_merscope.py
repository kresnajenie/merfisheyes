#!/usr/bin/env python3
"""Turn a MERSCOPE (Vizgen) region export into a folder MERFISHEYES opens as it is.

  python scripts/prepare_merscope.py <region_folder> <folder> [--name NAME]
                                     [--z 3] [--workers N]
                                     [--cluster] [--url BASE]

The region folder holds the files the MERSCOPE instrument (or vpt) writes:

  cell_metadata.csv          one row per cell: EntityID, center_x, center_y, ...
  cell_by_gene.csv           counts, one row per cell, one column per gene
  detected_transcripts.csv   (or .parquet) gene, global_x, global_y, global_z, ...
  cell_boundaries.parquet    (optional) vpt polygons: EntityID, ZIndex, Geometry (WKB)
  cell_categories.csv        (optional) cell types / clusters per cell

It writes, with the `spatialdata` library (so the store is one that library
reads back):

  <folder>/<name>.zarr/   tables/table (counts without blanks, CSC; obsm/spatial =
                          the cell centres, microns; the metadata and categories as
                          obs columns), shapes/cell_boundaries (one GeoParquet: each
                          cell's polygon on one z-plane, --z, as microns)
                          and points/molecules (tiled, straight from the
                          transcripts table, flat on the cells' plane)
  <folder>/mapping.json   names the store
  <folder>/annotations/   with --cluster, Leiden clusters + DEGs when the table
                          has no categorical column

then runs scripts/prepare_folder.py's checks and scripts/validate_store.py.

Why not spatialdata-io: its MERSCOPE reader (0.7.1) stops without an
`images/` folder (the mosaic TIFFs and micron_to_mosaic_pixel_transform.csv),
which most shared exports lack, and its fixes are still an open pull request.
The mosaic images themselves are not written here yet: a MERSCOPE mosaic is
tens of gigabytes, more than scripts/image_to_ome_zarr.py holds in memory.

Needs: pip install -r scripts/requirements.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_folder  # noqa: E402
from meyes_common import should_filter_gene  # noqa: E402

ID_COLUMNS = ("EntityID", "cell_id", "cell", "id")


def find(folder: Path, *names: str) -> Path | None:
    for n in names:
        if (folder / n).exists():
            return folder / n
    return None


def id_column(df) -> str:
    """The cell id column: a known name, else the first column."""
    for c in ID_COLUMNS:
        if c in df.columns:
            return c
    return df.columns[0]


def read_table(folder: Path, log):
    """An AnnData from cell_metadata.csv + cell_by_gene.csv (+ cell_categories.csv)."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    meta_path = find(folder, "cell_metadata.csv")
    counts_path = find(folder, "cell_by_gene.csv")
    if meta_path is None or counts_path is None:
        raise SystemExit("cell_metadata.csv and cell_by_gene.csv are needed")

    meta = pd.read_csv(meta_path, dtype={c: str for c in ID_COLUMNS})
    if meta.columns[0].startswith("Unnamed"):
        meta = meta.rename(columns={meta.columns[0]: "EntityID"})
        meta["EntityID"] = meta["EntityID"].astype(str)
    idc = id_column(meta)
    ids = meta[idc].astype(str)
    log(f"  {meta_path.name}: {len(meta):,} cells, id column {idc!r}")

    counts = pd.read_csv(counts_path, index_col=0, dtype={0: str})
    counts.index = counts.index.astype(str)
    genes = [g for g in counts.columns if not should_filter_gene(g)]
    log(f"  {counts_path.name}: {counts.shape[0]:,} cells × {len(genes)} genes ({counts.shape[1] - len(genes)} blanks dropped)")
    if set(counts.index) == set(ids):
        counts = counts.loc[ids]  # the metadata's order
    elif len(counts) == len(meta):
        log("  cell ids differ between the two files; rows taken in file order")
        counts.index = ids.to_numpy()
    else:
        raise SystemExit(f"{meta_path.name} has {len(meta):,} cells, {counts_path.name} {len(counts):,}: not the same cells")

    for c in ("center_x", "center_y"):
        if c not in meta.columns:
            raise SystemExit(f"{meta_path.name} has no {c} column: where are the cells?")
    obs = meta.drop(columns=[idc]).copy()
    obs.index = pd.Index(ids.to_numpy(), name="cell_id")
    obs.insert(0, "cell_id", ids.to_numpy())

    cats_path = find(folder, "cell_categories.csv")
    if cats_path is not None:
        cats = pd.read_csv(cats_path, dtype=str)
        cidc = id_column(cats)
        cats = cats.set_index(cats[cidc].astype(str)).drop(columns=[cidc])
        for c in cats.columns:
            obs[c] = pd.Categorical(cats[c].reindex(obs.index))
        log(f"  {cats_path.name}: {list(cats.columns)}")

    x = sp.csc_matrix(counts[genes].to_numpy(dtype=np.float32))
    return ad.AnnData(
        X=x,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(genes, name="gene")),
        obsm={"spatial": meta[["center_x", "center_y"]].to_numpy(dtype=np.float64)},
    )


def read_boundaries(path: Path, z: int, log):
    """One polygon per cell from a vpt cell_boundaries.parquet: the plane
    ZIndex == z (the largest part of a MultiPolygon), indexed by EntityID."""
    import geopandas as gpd
    import pyarrow.parquet as pq
    import shapely

    pf = pq.ParquetFile(path)
    cols = {c.lower(): c for c in pf.schema_arrow.names}
    geom_col, id_col = cols.get("geometry"), cols.get("entityid")
    if geom_col is None or id_col is None:
        raise SystemExit(f"{path.name}: no Geometry / EntityID columns")
    want = [id_col, geom_col] + ([cols["zindex"]] if "zindex" in cols else [])
    table = pq.read_table(path, columns=want, filters=[(cols["zindex"], "==", z)] if "zindex" in cols else None)
    if table.num_rows == 0:
        planes = pq.read_table(path, columns=[cols["zindex"]]).column(0).unique().to_pylist()
        raise SystemExit(f"{path.name}: no polygons on ZIndex {z}; planes present: {sorted(planes)}")
    geoms = shapely.from_wkb(table.column(geom_col).to_numpy(zero_copy_only=False))
    largest = [
        max(g.geoms, key=lambda p: p.area) if g.geom_type == "MultiPolygon" else g
        for g in geoms
    ]
    ids = [str(v) for v in table.column(id_col).to_pylist()]
    gdf = gpd.GeoDataFrame(index=ids, geometry=largest)
    gdf = gdf[~gdf.index.duplicated()]
    log(f"  {path.name}: {len(gdf):,} polygons on ZIndex {z}")
    return gdf


def write_store(store: Path, table, boundaries, log):
    try:
        import spatialdata
        from spatialdata.models import ShapesModel, TableModel
        from spatialdata.transformations import Identity
    except ImportError:
        raise SystemExit("pip install spatialdata")
    import pandas as pd

    shapes = {}
    if boundaries is not None:
        shapes["cell_boundaries"] = ShapesModel.parse(boundaries, transformations={"global": Identity()})
        table.obs["region"] = pd.Categorical(["cell_boundaries"] * table.n_obs)
        table = TableModel.parse(table, region="cell_boundaries", region_key="region", instance_key="cell_id")
    else:
        table = TableModel.parse(table)
    log(f"writing {store} ...")
    spatialdata.SpatialData(tables={"table": table}, shapes=shapes).write(store)


def prepare_merscope(
    path: Path,
    folder: Path,
    *,
    name: str | None = None,
    z: int = 3,
    workers: int = 1,
    cluster: bool = False,
    url: str | None = None,
    expanded: bool = False,
    label_csv: Path | None = None,
    keep_z: bool = False,
    log=print,
):
    path, folder = Path(path), Path(folder)
    name = name or path.name
    folder.mkdir(parents=True, exist_ok=True)
    store = folder / f"{name}.zarr"
    if store.exists():
        raise SystemExit(f"{store} exists already; remove it or pass --name")

    log(f"reading {path} ...")
    table = read_table(path, log)
    boundaries_path = find(path, "cell_boundaries.parquet")
    boundaries = read_boundaries(boundaries_path, z, log) if boundaries_path else None
    if boundaries is not None:
        missing = len(set(table.obs["cell_id"]) - set(boundaries.index))
        if missing:
            log(f"  {missing:,} cells have no polygon on ZIndex {z}")
    write_store(store, table, boundaries, log)

    # The transcripts go straight from the instrument's table into the store's
    # tiled points element (or, --expanded, a per-gene folder beside it)
    transcripts = find(path, "detected_transcripts.parquet", "detected_transcripts.csv")
    if transcripts is None:
        log("  no detected_transcripts file: no molecules")
    log("checking the folder ...")
    return prepare_folder.prepare(
        folder, cells=store, molecules=transcripts, molecule_columns="merscope", cluster=cluster, url=url,
        expanded=expanded, workers=workers, labels=label_csv, keep_z=keep_z, log=log,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("region", type=Path, help="the MERSCOPE region folder (holds cell_metadata.csv)")
    ap.add_argument("folder", type=Path, help="the folder to write; it is what you drop, save or sync")
    ap.add_argument("--name", help="the store's name (default: the region folder's)")
    ap.add_argument("--z", type=int, default=3, help="the z-plane (ZIndex) whose polygons are kept (default 3)")
    ap.add_argument("--workers", type=int, default=4, help="parallel gene-file writers")
    ap.add_argument("--cluster", action="store_true", help="Leiden clusters + DEGs when the table has no cell types")
    ap.add_argument("--url", help="the folder's public URL, for absolute links in mapping.json")
    ap.add_argument("--expanded", action="store_true", help="keep the molecule folder beside the store instead of inside it")
    args = ap.parse_args(argv)
    r = prepare_merscope(
        args.region,
        args.folder,
        name=args.name,
        z=args.z,
        workers=args.workers,
        cluster=args.cluster,
        url=args.url,
        expanded=args.expanded,
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
