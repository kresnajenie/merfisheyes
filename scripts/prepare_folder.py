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


# ── the table ─────────────────────────────────────────────────────────────

def drop_consolidated(zarr_json: Path) -> None:
    """Stale consolidated metadata is removed before a group is edited:
    anndata refuses to write into a consolidated group, and the root's is
    written again, last, once everything is in place."""
    meta = read_json(zarr_json)
    if meta and "consolidated_metadata" in meta:
        del meta["consolidated_metadata"]
        zarr_json.write_text(json.dumps(meta, indent=2))


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
    obs = table["obs"]
    index = obs.attrs.get("_index", "_index")
    return any(
        obs[k].attrs.get("encoding-type") == "categorical"
        for k in obs
        if k != index and hasattr(obs[k], "attrs")
    )


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
    if cells is not None:
        kind = kind_of(cells)
        if kind == "store":
            store = cells
        else:
            store = folder / f"{name or folder.name}.zarr"
            write_store_from_table(cells, kind, store, log)
        table = table_group(store)
        ensure_spatial(table, log)
        ensure_csc(table, log)
        if cluster and not has_categorical_column(table):
            import spatialdata_cluster

            log("  no categorical column: Leiden clustering ...")
            spatialdata_cluster.cluster(store, folder / "annotations", resolution=resolution, log=log)
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
        r.ok("molecules", f"{molecules.name}: manifest and genes/ present")
    if image is not None:
        root = read_json(image / "zarr.json") or {}
        levels = root.get("attributes", {}).get("ome", {}).get("multiscales", [{}])[0].get("datasets", [])
        r.ok("image", f"{image.name}: {len(levels)} level(s)")
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
    )
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
