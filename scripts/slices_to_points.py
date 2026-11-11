#!/usr/bin/env python3
"""Write every slice's molecules into one tiled points element, in the combined cells' frame.

  python scripts/slices_to_points.py <combined_output> <store.zarr>
                                     [--name molecules] [--keep-z]
                                     [--columns ...] [--levels 8] [--tile-size 100]

The cells of a multi-slice dataset are combined into one table by
scripts/combine_slices_v3.py: it grid-shifts each slice so the slices do not
overlap, tags every cell with `_sample_id`, and records each slice's exact
shift — and the path to that slice's molecules — in
`<combined_output>/sample_offsets.json`.

This reads that file and does the same thing to the molecules: each slice's
transcripts table is moved by its slice's (dx, dy) and written into one
points element, so molecules and cells share one coordinate system and one
store. Every molecule keeps its slice's id in a `sample` column, matching the
cells' `_sample_id`.

The molecules are never concatenated into one table on disk: they are
streamed slice by slice into the writer, which spills to a temporary
directory and builds the levels once. A slice with no transcripts file is
skipped (cells only).

Exit codes: 0 written; 3 nothing to write (no slice has molecules).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import transcripts_to_points as t2p  # noqa: E402

OFFSETS = "sample_offsets.json"


def read_offsets(combined: Path) -> list[dict]:
    """The slices of `combined`, each with its shift and its molecules."""
    path = combined / OFFSETS if combined.is_dir() else combined
    if not path.exists():
        raise SystemExit(
            f"{path}: no {OFFSETS}. Run scripts/combine_slices_v3.py first — it "
            "writes the per-slice shifts this stage places the molecules by."
        )
    data = json.loads(path.read_text())
    samples = data.get("samples") or []
    if not samples:
        raise SystemExit(f"{path}: no slices")
    return samples


def sources(samples: list[dict], log=print) -> list[dict]:
    """The slices that have molecules, as convert_many takes them."""
    out = []
    for s in samples:
        tx = s.get("transcripts")
        if not tx:
            log(f"  {s['sample_id']}: no molecules, skipped")
            continue
        path = Path(tx)
        if not path.exists():
            raise SystemExit(
                f"{path}: the molecules {OFFSETS} names are not there. "
                "The combine ran somewhere else, or the slice has moved."
            )
        out.append({
            "path": path,
            "dx": float(s["dx"]),
            "dy": float(s["dy"]),
            "sample_id": str(s["sample_id"]),
        })
    return out


def refresh_stats(store: Path) -> None:
    """Bring `stats.json` back in line with the store.

    scripts/prepare.py writes that sidecar when it writes the cells, before
    this stage adds the molecules — so left alone it says `"molecules": 0`
    about a store holding hundreds of millions. It is what a catalog reads,
    so it is rewritten here from the store itself.
    """
    import prepare

    path = store.parent / "stats.json"
    if not path.exists():
        return
    stats = json.loads(path.read_text())
    stats.update(prepare.store_stats(store))
    path.write_text(json.dumps(stats, indent=2) + "\n")
    print(f"stats.json: molecules {stats['molecules']:,}, genes {stats['molecule_genes']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("combined", type=Path, help=f"the combined_output folder (holding {OFFSETS})")
    ap.add_argument("store", type=Path, help="the SpatialData store to write points/<name> into")
    ap.add_argument("--name", default="molecules", help="the element's name under points/ (default molecules)")
    ap.add_argument("--columns", help="xenium, merscope, GENE,X,Y[,Z], or gene=G,x=X,y=Y[,z=Z][,cell=C]")
    ap.add_argument("--keep-z", action="store_true", help="keep z (every slice must have it)")
    ap.add_argument("--levels", type=int, default=8)
    ap.add_argument("--tile-size", type=float, default=100.0)
    args = ap.parse_args(argv)

    samples = read_offsets(args.combined)
    print(f"{len(samples)} slice(s) in {args.combined}")
    srcs = sources(samples)
    if not srcs:
        print("no slice has molecules: nothing to write", file=sys.stderr)
        return 3

    t2p.convert_many(
        srcs,
        args.store,
        name=args.name,
        columns=args.columns,
        keep_z=args.keep_z,
        levels=args.levels,
        tile_size=args.tile_size,
    )
    refresh_stats(args.store)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
