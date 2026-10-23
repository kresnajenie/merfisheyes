#!/usr/bin/env python3
"""Prepare a folder of already-processed parts for MERFISHEYES.

  python scripts/prepare_folder.py <folder> [--name NAME] [--url BASE]
                                   [--cluster] [--resolution 0.5]
                                   [--cells PATH] [--molecules PATH] [--image PATH]

The folder holds any of three parts, each recognised by what is in it:

  cells      an AnnData table — a `.h5ad` file or a Zarr folder (v2 or v3)
             written by `write_zarr` — or a SpatialData store (`zarr.json` + `tables/`)
  molecules  a merfisheyes molecule folder (`manifest.json[.gz]` + `genes/`),
             from scripts/process_single_molecule.py
  image      an OME-Zarr image folder (`zarr.json` naming axes c, y, x)

Afterwards the folder opens as it is — dropped on the homepage's Folder card,
saved with Upload & Save, or synced to a public bucket — and holds:

  <name>.zarr/     one SpatialData store (Zarr v3): tables/table with `layers/X_csc`
                   added when X is not already CSC (the viewer reads one gene as
                   one slice; a CSR X cannot be read per gene) and `obsm/spatial`;
                   images/<image> written by the spatialdata library from the
                   image folder; points/molecules, a tiled points element the
                   viewer streams by view, from the molecule folder; consolidated
                   metadata written last
  sources/         the image and molecule folders the store was made from
  mapping.json     names the store (--url BASE writes an absolute URL for a bucket)
  annotations/     with --cluster, the per-cluster DEG stats; the clusters
                   themselves go into the table as obs/leiden
                   (scripts/spatialdata_cluster.py)

With --expanded the image and the molecule folder stay beside the store and
mapping.json names them, relative to the folder; that layout opens the same way.

A SpatialData store that is there already is repaired where the viewer
cannot read it as written: a shapes element stored as a directory of parquet
parts or in the geoarrow encoding becomes one WKB file (the original kept
as `shapes.parquet.parts`), image levels in chunks over 1024² are rewritten
in 1024² chunks, and with --molecules-from-points a points element becomes
the molecule folder (the store's own parquet points are not read).

The last step runs scripts/validate_store.py on the store and the mapping;
the exit code is 1 when a MUST fails. A Zarr v3 table is moved into the store
(no rewrite); a Zarr v2 table or an .h5ad is read whole and written again, and
the original is left in place — remove it before dropping the folder, or the
viewer sees two candidates for the cells.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import validate_store  # noqa: E402
from meyes_common import should_filter_gene  # noqa: E402

SPATIALDATA_ATTRS = {"version": "0.2"}


# ── what is in the folder ────────────────────────────────────────────────

def read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def transcripts_columns(path: Path) -> bool:
    """True when the .parquet / .csv file has a gene column and x, y columns:
    a platform's names (Xenium, MERSCOPE) or a lab's own (Gene, X, Y), by name
    in any case (meyes_common.match_molecule_columns); anything else needs
    --molecule-columns."""
    try:
        import transcripts_to_points
        from meyes_common import match_molecule_columns

        names = list(transcripts_to_points.column_names(path))
    except Exception:  # noqa: BLE001 — not a table we can read
        return False
    known = any({g, x, y} <= set(names) for g, x, y, _ in transcripts_to_points.KNOWN.values())
    return known or match_molecule_columns(names) is not None


def kind_of(path: Path, molecule_columns: str | None = None) -> str | None:
    """`store`, `table3`, `table2`, `h5ad`, `molecules` (a per-gene folder),
    `transcripts` (a table of molecules), `image`, or None."""
    if path.is_file():
        if path.suffix == ".h5ad":
            return "h5ad"
        if path.suffix.lower() in (".parquet", ".csv") and (molecule_columns or transcripts_columns(path)):
            return "transcripts"
        return None
    if not path.is_dir():
        return None
    root = read_json(path / "zarr.json")
    if root is not None:
        attrs = root.get("attributes", {})
        if "spatialdata_attrs" in attrs or (path / "tables").is_dir():
            return "store"
        if (path / "obs" / "zarr.json").exists() and (path / "var" / "zarr.json").exists():
            return "table3"
        axes = [a.get("name") for a in attrs.get("ome", {}).get("multiscales", [{}])[0].get("axes", [])]
        if axes == ["c", "y", "x"]:
            return "image"
        return None
    if (path / ".zgroup").exists() and (path / "obs" / ".zgroup").exists():
        return "table2"
    if ((path / "manifest.json.gz").exists() or (path / "manifest.json").exists()) and (path / "genes").is_dir():
        return "molecules"
    return None


def find_parts(folder: Path, molecule_columns: str | None = None) -> dict[str, list[Path]]:
    """Every part in the folder, or the folder itself, grouped by what it is.
    Molecules are a per-gene folder or a transcripts table (parquet / csv)."""
    found: dict[str, list[Path]] = {"cells": [], "molecules": [], "image": []}
    candidates = [folder] + sorted(p for p in folder.iterdir() if not p.name.startswith("."))
    for p in candidates:
        kind = kind_of(p, molecule_columns)
        if kind in ("store", "table3", "table2", "h5ad"):
            found["cells"].append(p)
        elif kind in ("molecules", "transcripts"):
            found["molecules"].append(p)
        elif kind == "image":
            found["image"].append(p)
    return found


# ── where each part sits (the same units and frame, or nothing lines up) ──

Box = tuple[float, float, float, float]  # x0, y0, x1, y1


def table_box(table) -> Box | None:
    import anndata as ad

    if "spatial" not in table["obsm"]:
        return None
    xy = ad.io.read_elem(table["obsm/spatial"])
    xy = xy[:, :2]
    return float(xy[:, 0].min()), float(xy[:, 1].min()), float(xy[:, 0].max()), float(xy[:, 1].max())


def molecules_box(folder: Path, genes: int = 3) -> Box | None:
    """The extent of the first few genes' coordinates (raw, as the files hold them)."""
    import gzip

    import numpy as np

    manifest_path = next((folder / n for n in ("manifest.json.gz", "manifest.json") if (folder / n).exists()), None)
    if manifest_path is None:
        return None
    with (gzip.open(manifest_path, "rt") if manifest_path.suffix == ".gz" else manifest_path.open()) as f:
        manifest = json.load(f)
    dims = int(manifest.get("statistics", {}).get("spatial_dimensions", 2))
    scale = float(manifest.get("processing", {}).get("scaling_factor", 1) or 1)
    xs, ys = [], []
    for path in sorted((folder / "genes").glob("*.bin.gz"))[:genes]:
        with gzip.open(path, "rb") as f:
            flat = np.frombuffer(f.read(), dtype="<f4")
        if len(flat) < dims:
            continue
        coords = flat[: len(flat) // dims * dims].reshape(-1, dims)
        xs.append(coords[:, 0])
        ys.append(coords[:, 1])
    if not xs:
        return None
    x, y = np.concatenate(xs) * scale, np.concatenate(ys) * scale
    return float(x.min()), float(y.min()), float(x.max()), float(y.max())


def image_box(image: Path) -> Box | None:
    """Level 0's scale and translation over its shape, as the viewer places it."""
    root = read_json(image / "zarr.json") or {}
    ms = (root.get("attributes", {}).get("ome", {}).get("multiscales") or [{}])[0]
    datasets = ms.get("datasets") or []
    if not datasets:
        return None
    level = read_json(image / str(datasets[0].get("path")) / "zarr.json") or {}
    shape = level.get("shape") or []
    if len(shape) != 3:
        return None
    sx = sy = 1.0
    tx = ty = 0.0
    for t in datasets[0].get("coordinateTransformations") or []:
        if t.get("type") == "scale":
            sy, sx = float(t["scale"][-2]), float(t["scale"][-1])
        elif t.get("type") == "translation":
            ty, tx = float(t["translation"][-2]), float(t["translation"][-1])
    return tx, ty, tx + sx * shape[2], ty + sy * shape[1]


def overlap(a: Box, b: Box) -> float:
    """The fraction of `a` covered by `b`."""
    w = min(a[2], b[2]) - max(a[0], b[0])
    h = min(a[3], b[3]) - max(a[1], b[1])
    area = (a[2] - a[0]) * (a[3] - a[1])
    return 0.0 if w <= 0 or h <= 0 or area <= 0 else (w * h) / area


def fmt(box: Box) -> str:
    return f"x {box[0]:.4g}–{box[2]:.4g}, y {box[1]:.4g}–{box[3]:.4g}"


def check_against_cells(what: str, box: Box | None, cells: Box | None, r) -> None:
    if box is None or cells is None:
        return
    covered = overlap(cells, box)
    span = max(box[2] - box[0], box[3] - box[1]) / max(cells[2] - cells[0], cells[3] - cells[1], 1e-9)
    if covered == 0:
        r.warn(what, f"spans {fmt(box)} while the cells span {fmt(cells)}: they do not overlap — not the same units or frame?")
    elif span > 5 or span < 0.2:
        r.warn(what, f"spans {fmt(box)}, {span:.3g}× the cells' extent ({fmt(cells)}): not the same units?")
    else:
        r.ok(what, f"spans {fmt(box)}, covering {covered:.0%} of the cells' extent")


def check_image_folder(image: Path, r) -> None:
    """The image's own metadata through the validator's image rules
    (docs/specs/image-profile.md), as if it sat inside a store."""
    root = read_json(image / "zarr.json") or {}
    ms = (root.get("attributes", {}).get("ome", {}).get("multiscales") or [{}])[0]
    nodes = {f"images/{image.name}": root}
    for ds in ms.get("datasets") or []:
        level = read_json(image / str(ds.get("path")) / "zarr.json")
        if level:
            nodes[f"images/{image.name}/{ds['path']}"] = level
    validate_store.check_image(validate_store.RootStore({"consolidated_metadata": {"metadata": nodes}}), image.name, r)


# ── repairs of a store the viewer cannot read as written ──────────────────

CHUNK = 1024


def element_dirs(store: Path, kind: str) -> list[Path]:
    d = store / kind
    return sorted(p for p in d.iterdir() if p.is_dir() and (p / "zarr.json").is_file()) if d.is_dir() else []


def repair_shapes(store: Path, log) -> list[str]:
    """One WKB GeoParquet file per shapes element (shapes-profile.md). A
    directory of parts (dask-geopandas, Celldega's tiles) or a file in the
    geoarrow encoding is read whole with geopandas and written again as
    `shapes.parquet`; the parts are kept beside it as `shapes.parquet.parts`."""
    import pyarrow.parquet as pq

    repaired = []
    for el in element_dirs(store, "shapes"):
        path = el / "shapes.parquet"
        if path.is_dir():
            why = f"a directory of {sum(1 for _ in path.iterdir())} parts"
        elif path.is_file():
            meta = pq.read_metadata(path).metadata or {}
            geo = json.loads(meta.get(b"geo", b"{}") or b"{}")
            encoding = (geo.get("columns") or {}).get(geo.get("primary_column", ""), {}).get("encoding", "WKB")
            if str(encoding).upper() == "WKB":
                continue
            why = f"{encoding} geometry encoding"
        else:
            continue
        try:
            import geopandas as gpd
        except ImportError:
            raise SystemExit(f"shapes/{el.name} is {why}: pip install geopandas to rewrite it as one WKB file")
        log(f"  shapes/{el.name}: {why}; writing one WKB file ...")
        gdf = gpd.read_parquet(path)
        tmp = el / "shapes.parquet.tmp"
        gdf.to_parquet(tmp, geometry_encoding="WKB")
        parts = el / "shapes.parquet.parts"
        if parts.exists():
            shutil.rmtree(parts)
        path.rename(parts)
        tmp.rename(path)
        log(f"  shapes/{el.name}: {len(gdf):,} geometries in shapes.parquet (the original kept as shapes.parquet.parts)")
        repaired.append(el.name)
    return repaired


def repair_image_chunks(store: Path, log, chunk: int = CHUNK) -> list[str]:
    """Image levels in chunks of at most `chunk`² (image-profile.md): a level
    stored in larger chunks is rewritten, block by block, with the same
    codecs; the viewer reads one tile at a time and a 4096² uint16 tile is
    32 MB."""
    import zarr

    repaired = []
    for el in element_dirs(store, "images"):
        meta = read_json(el / "zarr.json") or {}
        ms = (meta.get("attributes", {}).get("ome", meta.get("attributes", {})).get("multiscales") or [{}])[0]
        for ds in ms.get("datasets") or []:
            level = el / str(ds.get("path"))
            src = zarr.open_array(level, mode="r")
            if len(src.shape) < 2 or max(src.chunks[-2:]) <= chunk:
                continue
            chunks = tuple([1] * (len(src.shape) - 2) + [min(chunk, src.shape[-2]), min(chunk, src.shape[-1])])
            log(f"  images/{el.name}/{level.name}: chunks {list(src.chunks)} → {list(chunks)} ...")
            tmp = el / f"{level.name}.rechunk"
            if tmp.exists():
                shutil.rmtree(tmp)
            dst = zarr.create_array(
                tmp,
                shape=src.shape,
                dtype=src.dtype,
                chunks=chunks,
                compressors=src.compressors,
                serializer=src.serializer,
                filters=src.filters,
                fill_value=src.fill_value,
                dimension_names=src.metadata.dimension_names,
                attributes=dict(src.attrs),
                zarr_format=3,
            )
            step = src.chunks  # read the source a chunk at a time
            import itertools

            ranges = [range(0, n, s) for n, s in zip(src.shape, step)]
            for start in itertools.product(*ranges):
                sel = tuple(slice(a, min(a + s, n)) for a, s, n in zip(start, step, src.shape))
                dst[sel] = src[sel]
            shutil.rmtree(level)
            tmp.rename(level)
            repaired.append(f"{el.name}/{level.name}")
    return repaired


# ── the wrapper: image and molecules inside the store ─────────────────────

def wrap_image(store: Path, image: Path, log) -> str:
    """The OME-Zarr image folder becomes `images/<name>` of the store, written
    by the spatialdata library from its finest level: the pyramid rebuilt in
    1024² chunks, the level's scale and translation as the element's
    transform into `global` (the cells' system)."""
    try:
        import dask.array as da
        import spatialdata
        from spatialdata.models import Image2DModel
        from spatialdata.transformations import Scale, Sequence, Translation
    except ImportError:
        raise SystemExit("pip install spatialdata   (to put the image inside the store; or pass --expanded)")

    name = image.name.removesuffix(".zarr")
    root = read_json(image / "zarr.json") or {}
    ms = (root.get("attributes", {}).get("ome", {}).get("multiscales") or [{}])[0]
    datasets = ms.get("datasets") or []
    if not datasets:
        raise SystemExit(f"{image.name}: no multiscale levels")
    sy = sx = 1.0
    ty = tx = 0.0
    for t in datasets[0].get("coordinateTransformations") or []:
        if t.get("type") == "scale":
            sy, sx = float(t["scale"][-2]), float(t["scale"][-1])
        elif t.get("type") == "translation":
            ty, tx = float(t["translation"][-2]), float(t["translation"][-1])
    channels = [c.get("label") for c in root.get("attributes", {}).get("ome", {}).get("omero", {}).get("channels", [])]
    arr = da.from_zarr(str(image / str(datasets[0]["path"])))
    if arr.ndim == 2:
        arr = arr[None]
    levels = len(datasets)
    log(f"  images/{name}: {arr.shape} from {image.name}, {levels} level(s), {sx:g} units/px at ({tx:g}, {ty:g}) ...")
    element = Image2DModel.parse(
        arr,
        dims=("c", "y", "x"),
        c_coords=channels if channels and len(channels) == arr.shape[0] and all(channels) else None,
        transformations={"global": Sequence([Scale([sy, sx], axes=("y", "x")), Translation([ty, tx], axes=("y", "x"))])},
        scale_factors=[2] * (levels - 1) if levels > 1 else None,
        chunks=(1, CHUNK, CHUNK),
    )
    drop_consolidated(store / "zarr.json")
    sdata = spatialdata.read_zarr(store, selection=("tables",))
    sdata.images[name] = element
    sdata.write_element(name)
    return name


def region_transforms(store: Path, dims: int = 2) -> list[dict] | None:
    """The transforms the cells' region (the shapes element the default
    table annotates) declares into the store's named systems. Points are
    drawn in the cells' frame, so a points element made from molecules in
    the cells' units must declare the same (store-layout.md §8); with no
    region, identity into `global`."""
    table = read_json(store / "tables" / default_table_name(store) / "zarr.json") or {}
    region = table.get("attributes", {}).get("region")
    if not isinstance(region, str):
        return None
    el = read_json(store / "shapes" / region / "zarr.json") or {}
    out = []
    for t in el.get("attributes", {}).get("coordinateTransformations") or []:
        kept = {k: v for k, v in t.items() if k != "input"}
        for key in ("scale", "translation"):
            if key in kept and len(kept[key]) != dims:
                kept[key] = list(kept[key])[:dims]
        out.append(kept)
    return out or None


def wrap_molecules(store: Path, molecules: Path, name: str, workers: int, log, columns: str | None = None) -> str:
    """The molecules become `points/<name>`, a tiled points element the
    viewer streams by view (docs/points-zarr-spec.md), flat on the cells'
    plane and in the cells' frame: straight from a transcripts table, a
    store's own parquet points element, or a per-gene folder."""
    transforms = region_transforms(store)
    log(f"  points/{name}: from {molecules.name}" + (", transforms of the cells' region" if transforms else "") + " ...")
    if molecules.is_dir() and (molecules / "genes").is_dir():
        import legacy_to_points

        legacy_to_points.convert(molecules, store, name=name, drop_z=True, transformations=transforms, log=lambda *_: None)
    else:
        import transcripts_to_points

        transcripts_to_points.convert(molecules, store, name=name, columns=columns, transformations=transforms, log=log)
    return name


def expand_molecules(folder: Path, molecules: Path, columns: str | None, workers: int, log) -> Path:
    """The expanded layout keeps molecules as a per-gene folder beside the
    store: a transcripts table is written as one (process_single_molecule.py)."""
    if molecules.is_dir():
        return molecules
    import process_single_molecule as psm
    import transcripts_to_points

    gene, x, y, z = transcripts_to_points.resolve_columns(transcripts_to_points.column_names(molecules), columns)
    out = folder / "molecules"
    log(f"  {out.name}/: one file per gene from {molecules.name} ...")
    psm.process_single_molecule_data(
        input_file=str(molecules), output_folder=str(out), dataset_type="custom",
        gene_col=gene, x_col=x, y_col=y, z_col=None, num_workers=workers,
    )
    return out


# ── the table ─────────────────────────────────────────────────────────────

def drop_consolidated(zarr_json: Path) -> None:
    """Stale consolidated metadata is removed before a group is edited:
    anndata refuses to write into a consolidated group, and the root's is
    written again, last, once everything is in place."""
    meta = read_json(zarr_json)
    if meta and "consolidated_metadata" in meta:
        del meta["consolidated_metadata"]
        zarr_json.write_text(json.dumps(meta, indent=2))


def default_table_name(store: Path) -> str:
    """The table the viewer reads (store-layout.md): the curated
    `metadata_general["Default Table"]`, else `table`, else the first by name."""
    root = read_json(store / "zarr.json") or {}
    tables = sorted(p.name for p in (store / "tables").iterdir() if (p / "zarr.json").is_file()) if (store / "tables").is_dir() else []
    if not tables:
        raise SystemExit(f"{store.name} has no tables/")
    return validate_store.default_table(root.get("attributes", {}), tables) or tables[0]


def table_group(store: Path, name: str = "table"):
    import zarr

    drop_consolidated(store / "zarr.json")
    drop_consolidated(store / "tables" / name / "zarr.json")
    root = zarr.open_group(store, mode="r+", use_consolidated=False)
    return root[f"tables/{name}"]


def ensure_spatial(table, log) -> None:
    """`obsm/spatial` from `obsm/X_spatial`, or from coordinate columns of
    `obs` (center_x / center_y (/ center_z), x / y...) — where the legacy
    pipeline and the browser also look — when it is not there."""
    import anndata as ad
    import numpy as np
    from meyes_common import obs_coordinate_columns

    obsm = table.require_group("obsm")
    if "spatial" in obsm:
        return
    if "X_spatial" in obsm:
        ad.io.write_elem(obsm, "spatial", ad.io.read_elem(obsm["X_spatial"]))
        log("  obsm/spatial written from obsm/X_spatial")
        return
    obs = ad.io.read_elem(table["obs"])
    columns = obs_coordinate_columns(obs.columns)
    if columns is None:
        raise SystemExit(
            "the table has no obsm/spatial, obsm/X_spatial nor coordinate columns in obs "
            "(center_x / center_y, x / y...): where are the cells?"
        )
    ad.io.write_elem(obsm, "spatial", obs[columns].to_numpy(dtype=np.float64))
    log(f"  obsm/spatial written from obs columns {', '.join(columns)}")


def ensure_csc(table, log) -> bool:
    """A CSC copy of X in `layers/X_csc` unless X is CSC already, or a CSC
    layer of X's shape exists. Returns True when a layer was written."""
    import anndata as ad
    import scipy.sparse as sp

    if "X" not in table:
        log("  no X: no gene expression")
        return False
    x = table["X"]
    enc = x.attrs.get("encoding-type", "array")
    if str(enc).startswith("csc"):
        log("  X is CSC already")
        return False
    shape = list(x.attrs.get("shape") or x.shape)
    layers = table.require_group("layers")
    for name in layers:
        layer = layers[name]
        if str(layer.attrs.get("encoding-type", "")).startswith("csc") and list(layer.attrs.get("shape", [])) == shape:
            log(f"  layers/{name} is a CSC copy already")
            return False
    log(f"  X is {enc}: writing layers/X_csc ...")
    csc = sp.csc_matrix(ad.io.read_elem(x))
    ad.io.write_elem(layers, "X_csc", csc)
    return True


def filter_genes(adata, log):
    """`adata` without control probes, blanks and unassigned codewords: the
    same rule as every other reader (meyes_common.should_filter_gene)."""
    keep = [not should_filter_gene(g) for g in adata.var_names]
    dropped = len(keep) - sum(keep)
    if dropped:
        log(f"  {dropped} control / blank features of {len(keep)} left out")
        adata = adata[:, keep].copy()
    return adata


def write_table(tables, adata, log) -> None:
    """`tables/table` from `adata`. An `uns` entry the Zarr writer refuses (a
    string too long for an array: the BIL files' `uns['cell_polygons']`) is
    left out with a warning; everything else is written as it is."""
    import anndata as ad

    uns, adata.uns = dict(adata.uns), {}
    ad.io.write_elem(tables, "table", adata)
    group = tables["table"].require_group("uns")
    for key, value in uns.items():
        try:
            ad.io.write_elem(group, key, value)
        except Exception as e:  # noqa: BLE001 — whatever the writer refuses
            if key in group:
                del group[key]
            log(f"WARN uns['{key}'] left out: the Zarr writer refused it ({type(e).__name__}: {str(e)[:120]})")


def write_store_from_table(source: Path, kind: str, store: Path, log) -> None:
    """A lone table becomes `store/tables/table`: moved when it is Zarr v3,
    read and written again when it is Zarr v2 or an .h5ad."""
    import anndata as ad
    import zarr

    if store.exists():
        raise SystemExit(f"{store} exists already; remove it or pass --name")
    root = zarr.open_group(store, mode="w", zarr_format=3)
    root.attrs["spatialdata_attrs"] = SPATIALDATA_ATTRS
    tables = root.create_group("tables")
    adata = None
    if kind == "table3":
        log(f"  moving {source.name} into {store.name}/tables/table")
        shutil.move(str(source), str(store / "tables" / "table"))
        # Rewritten only when it holds control probes to drop
        var = ad.io.read_elem(zarr.open_group(store / "tables" / "table" / "var", mode="r"))
        if any(should_filter_gene(g) for g in var.index):
            adata = ad.io.read_elem(zarr.open_group(store / "tables" / "table", mode="r"))
            shutil.rmtree(store / "tables" / "table")
    else:
        log(f"  reading {source.name} and writing it as Zarr v3 ...")
        adata = ad.read_h5ad(source) if kind == "h5ad" else ad.read_zarr(source)
    if adata is not None:
        write_table(tables, filter_genes(adata, log), log)
    # What the spatialdata library writes on a table, so it reads this one back
    drop_consolidated(store / "tables" / "table" / "zarr.json")
    table = zarr.open_group(store / "tables" / "table", mode="r+", use_consolidated=False)
    table.attrs.update(
        {"spatialdata-encoding-type": "ngff:regions_table", "region": None, "region_key": None,
         "instance_key": None, "version": SPATIALDATA_ATTRS["version"]}
    )


def merge_labels(table, csv: Path, log) -> None:
    """Per-cell labels as `obs` columns: a MapMyCells result, or any CSV with
    one row per cell (`#` lines are comments). Rows are matched by their
    `cell_id` when those are the table's cell ids (its index or instance
    key), else taken in the table's order when the counts agree. Every other
    column that is not blank is written, text as categories."""
    import anndata as ad
    import pandas as pd

    labels = pd.read_csv(csv, comment="#")
    obs = ad.io.read_elem(table["obs"])
    keys = [obs.index.astype(str)]
    instance_key = table.attrs.get("instance_key")
    if instance_key and instance_key in obs.columns:
        keys.append(obs[instance_key].astype(str))
    ids = labels["cell_id"].astype(str) if "cell_id" in labels.columns else None
    matched = next((k for k in keys if ids is not None and set(ids) <= set(k)), None)
    if matched is not None:
        labels = labels.set_index(ids).reindex(matched.to_numpy())
        log(f"  {csv.name}: matched by cell_id ({labels.index.notna().sum():,} cells)")
    elif len(labels) == len(obs):
        log(f"  {csv.name}: no cell_id the table knows; rows taken in the table's order")
    else:
        raise SystemExit(f"{csv.name} has {len(labels):,} rows and the table {len(obs):,} cells, and no cell_id the table knows")

    group = table["obs"]
    order = list(group.attrs.get("column-order", []))
    for column in labels.columns:
        if column == "cell_id":
            continue
        values = labels[column]
        if values.isna().all() or (values.astype(str).str.strip() == "").all():
            continue
        values = values.to_numpy() if pd.api.types.is_numeric_dtype(values) else pd.Categorical(values.astype("string"))
        if column in group:
            log(f"  obs/{column} replaced by {csv.name}'s")
            del group[column]
        ad.io.write_elem(group, column, values)
        if column not in order:
            order.append(column)
    group.attrs["column-order"] = order


def has_categorical_column(table) -> bool:
    """A column the viewer would offer as cell types: categorical, with more
    than one category, and not the join keys (`region`, the instance key)."""
    obs = table["obs"]
    skip = {obs.attrs.get("_index", "_index"), table.attrs.get("region_key"), table.attrs.get("instance_key")}
    for k in obs:
        if k in skip or not hasattr(obs[k], "attrs") or obs[k].attrs.get("encoding-type") != "categorical":
            continue
        categories = obs[k].get("categories")
        if categories is None or categories.shape[0] > 1:
            return True
    return False


# ── the folder ────────────────────────────────────────────────────────────

def relative(part: Path | None, folder: Path, base: str | None) -> str | None:
    if part is None:
        return None
    name = part.name if part != folder else "."
    return f"{base.rstrip('/')}/{name}" if base else name


def prepare(
    folder: Path,
    *,
    name: str | None = None,
    url: str | None = None,
    cells: Path | None = None,
    molecules: Path | None = None,
    image: Path | None = None,
    cluster: bool = False,
    resolution: float = 0.5,
    molecules_from_points: bool = False,
    workers: int = 1,
    expanded: bool = False,
    molecule_columns: str | None = None,
    images: bool = True,
    labels: Path | None = None,
    log=print,
) -> validate_store.Report:
    folder = folder.resolve()
    parts = find_parts(folder, molecule_columns)

    def pick(kind: str, given: Path | None) -> Path | None:
        if given is not None:
            if kind_of(given, molecule_columns) is None:
                raise SystemExit(f"{given} is not a {kind} part this script knows")
            return given.resolve()
        if len(parts[kind]) > 1:
            names = ", ".join(p.name for p in parts[kind])
            raise SystemExit(f"several {kind} candidates ({names}); pass --{kind}")
        return parts[kind][0] if parts[kind] else None

    cells, molecules = pick("cells", cells), pick("molecules", molecules)
    # An image is used when there is one, unless images are left out
    image = pick("image", image) if images else None
    log(f"{folder}")
    log(f"  cells:     {cells.name if cells else '-'}")
    log(f"  molecules: {molecules.name if molecules else '-'}")
    log(f"  image:     {image.name if image else '-'}")
    if cells is None and molecules is None and image is None:
        raise SystemExit("nothing found: no table, store, molecule folder or image")

    store: Path | None = None
    cells_box: Box | None = None
    if cells is not None:
        kind = kind_of(cells)
        if kind == "store":
            store = cells
        else:
            store = folder / f"{name or folder.name}.zarr"
            write_store_from_table(cells, kind, store, log)
        table_name = default_table_name(store)
        log(f"  table: tables/{table_name}")
        table = table_group(store, table_name)
        ensure_spatial(table, log)
        ensure_csc(table, log)
        if labels is not None:
            merge_labels(table, labels, log)
        cells_box = table_box(table)
        if cluster and not has_categorical_column(table):
            import spatialdata_cluster

            log("  no categorical column: Leiden clustering ...")
            spatialdata_cluster.cluster(store, folder / "annotations", table=table_name, resolution=resolution, log=log)
        repair_shapes(store, log)
        repair_image_chunks(store, log)
        if molecules is None and molecules_from_points:
            parquet_points = [
                el for el in element_dirs(store, "points")
                if (el / "points.parquet").exists() and not (read_json(el / "zarr.json") or {}).get("attributes", {}).get("points:version")
            ]
            if parquet_points and expanded:
                import spatialdata_points_to_sm

                molecules = folder / "molecules"
                log(f"  points/{parquet_points[0].name} → {molecules.name} (flat) ...")
                spatialdata_points_to_sm.convert(store, molecules, points=parquet_points[0].name, workers=workers, drop_z=True)
            elif parquet_points and not (store / "points" / "molecules").exists():
                wrap_molecules(store, parquet_points[0], "molecules", workers, log)
        if expanded and molecules is not None:
            molecules = expand_molecules(folder, molecules, molecule_columns, workers, log)
        if not expanded and (image is not None or molecules is not None):
            # One SpatialData store holding everything; the parts it was made
            # from go to sources/ (where a drop does not look) when they sit in
            # this folder; a file given from elsewhere stays where it is
            sources = folder / "sources"

            def put_aside(part: Path) -> None:
                if part.parent == folder:
                    sources.mkdir(exist_ok=True)
                    shutil.move(str(part), str(sources / part.name))
                    log(f"  {part.name} moved to {sources.name}/")

            if image is not None:
                wrap_image(store, image, log)
                put_aside(image)
                image = None
            if molecules is not None:
                wrap_molecules(store, molecules, "molecules", workers, log, columns=molecule_columns)
                put_aside(molecules)
                molecules = None
        import zarr

        zarr.consolidate_metadata(str(store))
        log(f"  consolidated metadata written: {store.name}")

    mapping: dict = {}
    if store is not None:
        mapping["cells"] = relative(store, folder, url)
    if molecules is not None:
        mapping["linkColumn"] = "__all__"
        mapping["links"] = {"__all__": relative(molecules, folder, url)}
    if image is not None:
        mapping["image"] = relative(image, folder, url)
    (folder / "mapping.json").write_text(json.dumps(mapping, indent=2) + "\n")
    log(f"  mapping.json: {json.dumps(mapping)}")

    r = validate_store.Report()
    if store is not None:
        r = validate_store.validate(store, folder / "mapping.json")
    else:
        validate_store.check_mapping(folder / "mapping.json", r)
    if molecules is not None:
        r.ok(f"molecules/{molecules.name}", "manifest and genes/ present")
        check_against_cells(f"molecules/{molecules.name}", molecules_box(molecules), cells_box, r)
    if image is not None:
        check_image_folder(image, r)
        check_against_cells(f"images/{image.name}", image_box(image), cells_box, r)
    for level, where, msg in r.lines:
        log(f"{level:4} {where}: {msg}")
    return r


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("--name", help="name of the store written from a lone table (default: the folder's)")
    ap.add_argument("--url", help="the folder's public URL: mapping.json then holds absolute links")
    ap.add_argument("--cells", type=Path, help="the cells part, when the folder holds several")
    ap.add_argument("--molecules", type=Path)
    ap.add_argument("--image", type=Path)
    ap.add_argument("--cluster", action="store_true", help="Leiden clusters + DEGs when the table has no categorical column")
    ap.add_argument("--resolution", type=float, default=0.5)
    ap.add_argument("--molecules-from-points", action="store_true",
                    help="no molecule folder but the store has a points element: write one from it (flat)")
    ap.add_argument("--workers", type=int, default=4, help="parallel gene-file writers for --molecules-from-points")
    ap.add_argument("--expanded", action="store_true",
                    help="leave the image and the molecule folder beside the store, named in mapping.json, instead of inside it")
    ap.add_argument("--molecule-columns", help="a transcripts table's columns: xenium, merscope, or GENE,X,Y[,Z] (default: recognised)")
    ap.add_argument("--no-images", action="store_true", help="leave an image in the folder out of the store")
    ap.add_argument("--labels", type=Path, help="per-cell labels (a MapMyCells CSV) to add to the table's obs")
    args = ap.parse_args(argv)
    r = prepare(
        args.folder,
        name=args.name,
        url=args.url,
        cells=args.cells,
        molecules=args.molecules,
        image=args.image,
        cluster=args.cluster,
        resolution=args.resolution,
        molecules_from_points=args.molecules_from_points,
        workers=args.workers,
        expanded=args.expanded,
        molecule_columns=args.molecule_columns,
        images=not args.no_images,
        labels=args.labels,
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
