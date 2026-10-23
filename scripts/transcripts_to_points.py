#!/usr/bin/env python3
"""Write a transcripts table straight into the tiled points element the viewer streams.

  python scripts/transcripts_to_points.py <transcripts> <store.zarr> [--name molecules]
                                          [--columns xenium | merscope | GENE,X,Y[,Z]]
                                          [--keep-z] [--levels 8] [--tile-size 100]

<transcripts> is one of
  - a .parquet or .csv file, one row per molecule, with a gene column and x, y
    (, z) columns: Xenium's (feature_name, x_location, y_location, z_location)
    and MERSCOPE's (gene, global_x, global_y, global_z) are recognised; any
    other names are given with --columns;
  - a SpatialData points element (`points/<name>` of a store, or its
    `points.parquet` directory as spatialdata-io writes it).

The element written is docs/points-zarr-spec.md: level 0 holds every
molecule, each further level a uniform quarter of the one below; each level
in two row orders (by tile, by gene) with an index, so the viewer reads one
region or one gene by HTTP range requests. Control probes and blanks are
left out (the same filter as every other reader here). Coordinates are
written as they are in the table; the element's transforms into the store's
coordinate systems are the caller's (`convert(transformations=…)`), identity
into `global` by default. Molecules are flat on the plane unless --keep-z.

Two passes over the table, in batches: the genes and the extent, then the
molecules. 70 M rows take about a minute and 1 GB.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import points_zarr as pz  # noqa: E402
from meyes_common import match_molecule_columns, should_filter_gene  # noqa: E402

KNOWN = {
    "xenium": ("feature_name", "x_location", "y_location", "z_location"),
    "merscope": ("gene", "global_x", "global_y", "global_z"),
}


def parquet_source(source: Path) -> tuple[list[Path], str | None]:
    """The parquet files of `source` and, for a points element, its feature key."""
    if source.is_file():
        return [source], None
    if (source / "points.parquet").is_dir():
        attrs = {}
        try:
            attrs = json.loads((source / "zarr.json").read_text()).get("attributes", {})
        except (OSError, ValueError):
            pass
        source = source / "points.parquet"
        return sorted(source.glob("*.parquet")), (attrs.get("spatialdata_attrs") or {}).get("feature_key")
    if source.is_dir():
        files = sorted(source.glob("*.parquet"))
        if files:
            feature = None
            try:
                feature = json.loads((source.parent / "zarr.json").read_text())["attributes"]["spatialdata_attrs"]["feature_key"]
            except (OSError, ValueError, KeyError):
                pass
            return files, feature
    raise SystemExit(f"{source}: not a parquet / csv file nor a points element")


def column_names(source: Path) -> list[str]:
    if source.suffix.lower() == ".csv":
        import pandas as pd

        return list(pd.read_csv(source, nrows=0).columns)
    import pyarrow.parquet as pq

    files, _ = parquet_source(source)
    return pq.read_schema(files[0]).names


def resolve_columns(names: list[str], columns: str | None, feature_key: str | None = None) -> tuple[str, str, str, str | None]:
    """(gene, x, y, z | None) for the table: a known platform's, a points
    element's, or the names given as GENE,X,Y[,Z]."""
    if columns in KNOWN:
        g, x, y, z = KNOWN[columns]
        return g, x, y, z if z in names else None
    if columns:
        parts = [c.strip() for c in columns.split(",")]
        if len(parts) not in (3, 4):
            raise SystemExit("--columns is a platform name or GENE,X,Y[,Z]")
        missing = [c for c in parts if c not in names]
        if missing:
            raise SystemExit(f"columns {missing} are not in the table; it has {names}")
        return parts[0], parts[1], parts[2], parts[3] if len(parts) == 4 else None
    for g, x, y, z in KNOWN.values():
        if {g, x, y} <= set(names):
            return g, x, y, z if z in names else None
    if feature_key and {feature_key, "x", "y"} <= set(names):
        return feature_key, "x", "y", "z" if "z" in names else None
    # A lab's own table: the gene and coordinates by name, in any case
    found = match_molecule_columns(names)
    if found:
        return found["gene"], found["x"], found["y"], found.get("z")
    raise SystemExit(
        f"cannot tell the gene and coordinate columns of {names}: give --columns GENE,X,Y[,Z]"
    )


def batches(source: Path, cols: list[str], batch_size: int):
    """Dicts of numpy arrays, `batch_size` rows at a time."""
    if source.suffix.lower() == ".csv":
        import pandas as pd

        for chunk in pd.read_csv(source, usecols=cols, chunksize=batch_size):
            yield {c: chunk[c].to_numpy() for c in cols}
        return
    import pyarrow as pa
    import pyarrow.parquet as pq

    files, _ = parquet_source(source)
    for f in files:
        for batch in pq.ParquetFile(f).iter_batches(batch_size=batch_size, columns=cols):
            out = {}
            for c in cols:
                col = batch.column(c)
                # Xenium stores gene names as bytes; dictionary-encoded names are decoded
                if pa.types.is_dictionary(col.type):
                    col = col.dictionary_decode()
                if pa.types.is_binary(col.type) or pa.types.is_large_binary(col.type):
                    col = col.cast(pa.string())
                out[c] = col.to_numpy(zero_copy_only=False)
            yield out


def convert(
    source,
    store,
    *,
    name: str = "molecules",
    columns: str | None = None,
    keep_z: bool = False,
    levels: int = 8,
    tile_size: float = 100.0,
    transformations: list[dict] | None = None,
    batch_size: int = 2_000_000,
    log=print,
) -> Path:
    import pandas as pd
    import zarr

    source, store = Path(source), Path(store)
    feature_key = parquet_source(source)[1] if source.suffix.lower() != ".csv" else None
    gene_col, x_col, y_col, z_col = resolve_columns(column_names(source), columns, feature_key)
    dims = 3 if (z_col and keep_z) else 2
    cols = [gene_col, x_col, y_col] + ([z_col] if dims == 3 else [])
    if dims == 3 and transformations:
        # The cells' transforms are 2-D: z, in the same unit as x, is scaled
        # as x is and not moved
        transformations = [
            {**t, **{k: [*t[k], t[k][0] if k == "scale" else 0.0] for k in ("scale", "translation") if len(t.get(k) or []) == 2}}
            for t in transformations
        ]

    # Pass 1: the genes and the extent
    genes: set[str] = set()
    lo = np.array([np.inf, np.inf])
    hi = np.array([-np.inf, -np.inf])
    rows = 0
    for b in batches(source, [gene_col, x_col, y_col], batch_size):
        genes.update(str(g) for g in pd.unique(pd.Series(b[gene_col]).astype(str)))
        x, y = b[x_col].astype(np.float64), b[y_col].astype(np.float64)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.any():
            lo = np.minimum(lo, [x[ok].min(), y[ok].min()])
            hi = np.maximum(hi, [x[ok].max(), y[ok].max()])
        rows += len(x)
    if rows == 0 or not np.all(np.isfinite(lo)):
        raise SystemExit(f"{source.name}: no molecules with coordinates")
    kept = sorted(g for g in genes if not should_filter_gene(g))
    dropped = len(genes) - len(kept)
    code = {g: i for i, g in enumerate(kept)}
    pad = tile_size * 0.01
    origin = (float(lo[0]) - pad, float(lo[1]) - pad)
    extent = (float(hi[0] - lo[0]) + 2 * pad, float(hi[1] - lo[1]) + 2 * pad)
    log(f"{rows:,} molecules, {len(kept)} genes ({dropped} control/blank dropped), {dims}D, extent {extent[0]:.0f} × {extent[1]:.0f}")

    out = store / "points" / name
    if out.exists():
        raise SystemExit(f"{out} exists")
    root = zarr.open_group(store, mode="a")
    root.attrs.setdefault("spatialdata_attrs", {"version": "0.2"})
    root.require_group("points")
    w = pz.PointsWriter(
        out,
        store.parent / f".{name}.tmp",
        features=kept,
        dims=dims,
        origin=origin,
        extent=extent,
        tile_size=tile_size,
        n_levels=levels,
        transformations=transformations,
    )
    # Pass 2: the molecules, in batches, each gene by its code
    total = 0
    for i, b in enumerate(batches(source, cols, batch_size)):
        g = pd.Series(b[gene_col]).astype(str).map(code)
        keep = g.notna().to_numpy() & np.isfinite(b[x_col].astype(np.float64)) & np.isfinite(b[y_col].astype(np.float64))
        if not keep.any():
            continue
        spill = {
            "x": b[x_col][keep].astype(np.float64),
            "y": b[y_col][keep].astype(np.float64),
            "f": g.to_numpy()[keep].astype(w.fdtype),
        }
        if dims == 3:
            spill["z"] = b[z_col][keep].astype(np.float64)
        w.spill(spill, seed=i, tag=f"b{i}")
        total += int(keep.sum())
    log(f"{total:,} molecules spilled; building levels")
    w.finish()
    zarr.consolidate_metadata(str(store))
    log(f"wrote {out}")
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("transcripts", type=Path)
    ap.add_argument("store", type=Path)
    ap.add_argument("--name", default="molecules", help="the element's name under points/ (default molecules)")
    ap.add_argument("--columns", help="xenium, merscope, or GENE,X,Y[,Z] (default: recognised from the header)")
    ap.add_argument("--keep-z", action="store_true", help="keep z (for a 3D view of the molecules on their own)")
    ap.add_argument("--levels", type=int, default=8)
    ap.add_argument("--tile-size", type=float, default=100.0)
    args = ap.parse_args(argv)
    convert(args.transcripts, args.store, name=args.name, columns=args.columns, keep_z=args.keep_z,
            levels=args.levels, tile_size=args.tile_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
