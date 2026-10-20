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

  <name>.zarr/     a SpatialData store (Zarr v3) with the table in tables/table:
                   `layers/X_csc` added when X is not already CSC (the viewer
                   reads one gene as one slice; a CSR X cannot be read per gene),
                   `obsm/spatial` present, consolidated metadata written last
  mapping.json     names the store, the molecules and the image, relative to the
                   folder (--url BASE writes absolute URLs for a bucket)
  annotations/     with --cluster, Leiden clusters and per-cluster DEG stats when
                   the table has no categorical column (scripts/spatialdata_cluster.py)

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

SPATIALDATA_ATTRS = {"version": "0.2"}


# ── what is in the folder ────────────────────────────────────────────────

def read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def kind_of(path: Path) -> str | None:
    """`store`, `table3`, `table2`, `h5ad`, `molecules`, `image`, or None."""
    if path.is_file():
        return "h5ad" if path.suffix == ".h5ad" else None
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


def find_parts(folder: Path) -> dict[str, list[Path]]:
    """Every part in the folder, or the folder itself, grouped by what it is."""
    found: dict[str, list[Path]] = {"cells": [], "molecules": [], "image": []}
    candidates = [folder] + sorted(p for p in folder.iterdir() if not p.name.startswith("."))
    for p in candidates:
        kind = kind_of(p)
        if kind in ("store", "table3", "table2", "h5ad"):
            found["cells"].append(p)
        elif kind == "molecules":
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
    """`obsm/spatial` from `obsm/X_spatial` when only the latter is there."""
    import anndata as ad

    obsm = table.require_group("obsm")
    if "spatial" in obsm:
        return
    if "X_spatial" not in obsm:
        raise SystemExit("the table has no obsm/spatial (nor obsm/X_spatial): where are the cells?")
    ad.io.write_elem(obsm, "spatial", ad.io.read_elem(obsm["X_spatial"]))
    log("  obsm/spatial written from obsm/X_spatial")


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
    if kind == "table3":
        log(f"  moving {source.name} into {store.name}/tables/table")
        shutil.move(str(source), str(store / "tables" / "table"))
        return
    log(f"  reading {source.name} and writing it as Zarr v3 ...")
    adata = ad.read_h5ad(source) if kind == "h5ad" else ad.read_zarr(source)
    ad.io.write_elem(tables, "table", adata)


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
    log=print,
) -> validate_store.Report:
    folder = folder.resolve()
    parts = find_parts(folder)

    def pick(kind: str, given: Path | None) -> Path | None:
        if given is not None:
            if kind_of(given) is None:
                raise SystemExit(f"{given} is not a {kind} part this script knows")
            return given.resolve()
        if len(parts[kind]) > 1:
            names = ", ".join(p.name for p in parts[kind])
            raise SystemExit(f"several {kind} candidates ({names}); pass --{kind}")
        return parts[kind][0] if parts[kind] else None

    cells, molecules, image = pick("cells", cells), pick("molecules", molecules), pick("image", image)
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
        cells_box = table_box(table)
        if cluster and not has_categorical_column(table):
            import spatialdata_cluster

            log("  no categorical column: Leiden clustering ...")
            spatialdata_cluster.cluster(store, folder / "annotations", table=table_name, resolution=resolution, log=log)
        repair_shapes(store, log)
        repair_image_chunks(store, log)
        if molecules is None and molecules_from_points and element_dirs(store, "points"):
            import spatialdata_points_to_sm

            molecules = folder / "molecules"
            log(f"  points → {molecules.name} (flat) ...")
            spatialdata_points_to_sm.convert(store, molecules, workers=workers, drop_z=True)
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
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
