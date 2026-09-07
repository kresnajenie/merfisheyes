#!/usr/bin/env python3
"""Write a tiny renderable preview beside a labelled-molecule dataset.

    python export_preview.py /path/to/MER6-2_E3_1_lm [--points 5000] [--column cell]

Enough points to recognise the embryo and explore it in a thumbnail, small
enough that a rail of them costs less than one real dataset.

    preview/points.bin.gz   u32 version | u32 nPoints | u32 dims
                            f32 min[dims] | f32 scale[dims]
                            u16 coords[nPoints * dims]   (same q16 scheme as the
                                                          full coordinate file)
                            u8  colorIndex[nPoints]
    preview/index.json      the palette, in colour-index order

Colours are stored as an index into a palette, not RGB: one byte per point
instead of three, and it keeps the preview consistent with the viewer's own
colours if the palette is later edited.

Sampling is uniform, so the preview reproduces the embryo's density. Small cells
are consequently sparse — a polar body holding 0.02% of the molecules gets ~1
point at n=5000 — which is the right trade for a silhouette but means a preview
should never be read as showing every cell.
"""

import argparse
import gzip
import json
import struct
from pathlib import Path

import numpy as np


def read_obs(path: Path):
    raw = gzip.open(path, "rb").read()
    _, n, _, width, _, dict_len = struct.unpack_from("<IIBBHI", raw, 0)
    labels = json.loads(raw[16 : 16 + dict_len].decode())
    off = 16 + dict_len + ((-dict_len) % 4)
    codes = np.frombuffer(
        raw, {1: np.uint8, 2: np.uint16, 4: np.uint32}[width], n, off
    )

    return labels, codes


def read_coords(path: Path):
    raw = gzip.open(path, "rb").read()
    n, dims = struct.unpack_from("<II", raw, 0)
    off = 8
    mn = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    sc = np.frombuffer(raw, np.float32, dims, off)
    off += 4 * dims
    codes = np.frombuffer(raw, np.uint16, n * dims, off).reshape(n, dims)

    return mn + codes.astype(np.float32) * sc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dataset_dir")
    ap.add_argument("--points", type=int, default=5000)
    ap.add_argument("--column", default="cell")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    d = Path(args.dataset_dir)
    xyz = read_coords(d / "coords/spatial.q16.bin.gz")
    labels, codes = read_obs(d / f"obs/{args.column}.bin.gz")
    palette = json.loads((d / f"palettes/{args.column}.json").read_text())

    n = min(args.points, len(xyz))
    rng = np.random.default_rng(args.seed)  # fixed: the preview is deterministic
    pick = rng.choice(len(xyz), size=n, replace=False)
    pick.sort()

    P = xyz[pick]
    C = codes[pick].astype(np.uint8)

    if len(labels) > 255:
        raise SystemExit(
            f"{args.column} has {len(labels)} values; the colour index is one byte"
        )

    lo, hi = P.min(axis=0).astype(np.float32), P.max(axis=0).astype(np.float32)
    scale = ((hi - lo) / 65535.0).astype(np.float32)
    scale[scale == 0] = 1.0
    q = np.rint((P - lo) / scale).clip(0, 65535).astype(np.uint16)

    (d / "preview").mkdir(parents=True, exist_ok=True)
    with gzip.open(d / "preview/points.bin.gz", "wb") as f:
        f.write(struct.pack("<III", 1, n, P.shape[1]))
        f.write(lo.tobytes())
        f.write(scale.tobytes())
        f.write(np.ascontiguousarray(q).tobytes())
        f.write(C.tobytes())

    (d / "preview/index.json").write_text(json.dumps({
        "version": 1,
        "points": n,
        "column": args.column,
        "sampling": "uniform",
        # Index order matches the obs dictionary, so colourIndex maps straight in.
        "palette": [palette.get(l, "#808080") for l in labels],
        "labels": labels,
    }, indent=2))

    manifest_path = d / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["has_preview"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2))

    size = (d / "preview/points.bin.gz").stat().st_size
    idx = (d / "preview/index.json").stat().st_size
    print(f"  {d.name}: {n:,} points, {len(labels)} {args.column} values — "
          f"{size / 1024:.0f} KB + {idx / 1024:.0f} KB index")


if __name__ == "__main__":
    main()
