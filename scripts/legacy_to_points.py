#!/usr/bin/env python3
"""Convert a legacy molecule folder (manifest.json.gz + genes/<gene>.bin.gz,
as scripts/process_single_molecule.py writes it) into a multiscale points
element (docs/points-zarr-spec.md), in a new store or inside an existing
SpatialData store, so the viewer reads it a view at a time and the folder
can be retired.

  python scripts/legacy_to_points.py <molecule folder> <store.zarr>
                                     [--name transcripts] [--levels 8]
                                     [--tile-size 100] [--drop-z]

A folder written with normalized coordinates (no `coordinate_range:
raw_rounded_2dp`) is scaled back to microns with its scaling factor. The
unassigned files (`<gene>_uuuuuuuuuu.bin.gz`) are kept as a boolean
column `unassigned` of the element. Needs numpy and zarr >= 3.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import zarr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import points_zarr as pz  # noqa: E402

UNASSIGNED_SUFFIX = "_uuuuuuuuuu"


def sanitize_gene_name(gene: str) -> str:
    """A gene's file name, as the processor wrote it."""
    import re

    name = re.sub(r"[^a-zA-Z0-9]", "_", gene)
    name = re.sub(r"_+", "_", name)
    return name.strip("_")


def read_manifest(folder: Path) -> dict:
    with gzip.open(folder / "manifest.json.gz", "rt") as f:
        return json.load(f)


def read_gene(folder: Path, gene: str, unassigned=False) -> np.ndarray | None:
    """x, y, z of a gene's molecules as an (n, 3) float32 array; None when
    the folder has no file for it."""
    f = folder / "genes" / f"{sanitize_gene_name(gene)}{UNASSIGNED_SUFFIX if unassigned else ''}.bin.gz"
    if not f.is_file():
        return None
    with gzip.open(f, "rb") as g:
        raw = np.frombuffer(g.read(), dtype="<f4")
    return raw.reshape(-1, 3)


def convert(
    folder: Path,
    store: Path,
    *,
    name="transcripts",
    levels=8,
    tile_size=100.0,
    drop_z=False,
    transformations: list[dict] | None = None,
    log=print,
) -> Path:
    """`transformations`: the element's transforms into the store's named
    systems, as the cells' region declares them (points are drawn in the
    cells' frame); default identity into `global`."""
    manifest = read_manifest(folder)
    genes: list[str] = manifest["genes"]["unique_gene_names"]
    dims = int(manifest["statistics"]["spatial_dimensions"])
    raw = manifest.get("processing", {}).get("coordinate_range") == "raw_rounded_2dp"
    scale = 1.0 if raw else float(manifest.get("processing", {}).get("scaling_factor", 1))
    has_unassigned = bool(manifest.get("has_unassigned", False))
    if drop_z:
        dims = 2
    # Pass 1: the extent, from every file
    lo = np.array([np.inf, np.inf])
    hi = np.array([-np.inf, -np.inf])
    for gene in genes:
        for u in (False, True) if has_unassigned else (False,):
            xyz = read_gene(folder, gene, u)
            if xyz is None or len(xyz) == 0:
                continue
            lo = np.minimum(lo, xyz[:, :2].min(axis=0) * scale)
            hi = np.maximum(hi, xyz[:, :2].max(axis=0) * scale)
    if not np.all(np.isfinite(lo)):
        raise SystemExit("the folder holds no molecules")
    pad = tile_size * 0.01
    origin = (float(lo[0]) - pad, float(lo[1]) - pad)
    extent = (float(hi[0] - lo[0]) + 2 * pad, float(hi[1] - lo[1]) + 2 * pad)
    log(f"{len(genes)} genes, {dims}D, extent {extent[0]:.0f} × {extent[1]:.0f} from {origin}")

    store = Path(store)
    out = store / "points" / name
    if out.exists():
        raise SystemExit(f"{out} exists")
    root = zarr.open_group(store, mode="a")
    root.attrs.setdefault("spatialdata_attrs", {"version": "0.2"})
    root.require_group("points")
    columns = [("unassigned", "u1", None)] if has_unassigned else []
    w = pz.PointsWriter(
        out,
        store.parent / f".{name}.tmp",
        features=genes,
        dims=dims,
        columns=columns,
        copied=[c[0] for c in columns],
        origin=origin,
        extent=extent,
        tile_size=tile_size,
        n_levels=levels,
        transformations=transformations,
    )
    # Pass 2: every gene's molecules
    total = 0
    for i, gene in enumerate(genes):
        for u in (False, True) if has_unassigned else (False,):
            xyz = read_gene(folder, gene, u)
            if xyz is None or len(xyz) == 0:
                continue
            cols = {
                "x": xyz[:, 0] * scale,
                "y": xyz[:, 1] * scale,
                "f": np.full(len(xyz), i, dtype=w.fdtype),
            }
            if dims == 3:
                cols["z"] = xyz[:, 2] * scale
            if has_unassigned:
                cols["unassigned"] = np.full(len(xyz), int(u), dtype=np.uint8)
            w.spill(cols, seed=i * 2 + int(u), tag=f"{i}-{int(u)}")
            total += len(xyz)
    log(f"{total:,} molecules spilled; building levels")
    w.finish()
    zarr.consolidate_metadata(str(store))
    log(f"wrote {out}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("folder", type=Path)
    ap.add_argument("store", type=Path)
    ap.add_argument("--name", default="transcripts")
    ap.add_argument("--levels", type=int, default=8)
    ap.add_argument("--tile-size", type=float, default=100.0)
    ap.add_argument("--drop-z", action="store_true", help="write 2D molecules")
    args = ap.parse_args(argv)
    convert(args.folder, args.store, name=args.name, levels=args.levels, tile_size=args.tile_size, drop_z=args.drop_z)


if __name__ == "__main__":
    main()
