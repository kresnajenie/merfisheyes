#!/usr/bin/env python3
"""Write a TIFF (or any image tifffile reads) as the OME-Zarr image MERFISHEYES draws.

  python scripts/image_to_ome_zarr.py dapi.tif dapi.zarr --um-per-pixel 0.108
                                      [--origin X0 Y0] [--channels DAPI,PolyT]
                                      [--chunk 1024]

The image must be uint8 or uint16, one plane per channel: a 2-D array, a
(c, y, x) stack, or an RGB/RGBA (y, x, c) image. It is written as OME-Zarr
0.5 in Zarr v3 with axes c, y, x, chunks of 1024² (one tile is at most
2 MB), and a pyramid halved until the coarsest level fits in one tile.

Each level carries a `scale` (micrometres per pixel) and a `translation`
(where the centre of pixel (0, 0) sits among the cells, `--origin`, in the
same units as the cells' `obsm/spatial`), which is how the viewer places
it under the cells (docs/specs/image-profile.md). The whole image is held
in memory while it is written.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def as_cyx(image: np.ndarray) -> np.ndarray:
    """One plane per channel, channels first."""
    if image.ndim == 2:
        return image[None]
    if image.ndim == 3 and image.shape[0] <= 8:
        return image
    if image.ndim == 3 and image.shape[-1] <= 8:
        return np.moveaxis(image, -1, 0)
    raise SystemExit(f"cannot tell the channels of an image of shape {image.shape}: give (y, x), (c, y, x) or (y, x, c)")


def halve(level: np.ndarray) -> np.ndarray:
    """Every channel averaged 2×2, the odd row and column dropped."""
    c, h, w = level.shape
    h, w = h // 2 * 2, w // 2 * 2
    return level[:, :h, :w].reshape(c, h // 2, 2, w // 2, 2).mean(axis=(2, 4)).astype(level.dtype)


def write_ome_zarr(
    image: np.ndarray,
    out: Path,
    *,
    um_per_pixel: float,
    origin: tuple[float, float] = (0.0, 0.0),
    channels: list[str] | None = None,
    chunk: int = 1024,
    log=print,
) -> Path:
    import zarr

    image = as_cyx(np.asarray(image))
    if image.dtype not in (np.uint8, np.uint16):
        raise SystemExit(f"the image is {image.dtype}; the viewer reads uint8 or uint16 — convert it first")
    x0, y0 = origin
    names = channels or [f"channel {i}" for i in range(image.shape[0])]
    if len(names) != image.shape[0]:
        raise SystemExit(f"{len(names)} channel names for {image.shape[0]} channels")

    root = zarr.open_group(out, mode="w", zarr_format=3)
    datasets, level = [], image
    for i in range(32):
        root.create_array(str(i), data=level, chunks=(1, chunk, chunk))
        size = um_per_pixel * 2**i
        shift = (size - um_per_pixel) / 2  # a coarser pixel's centre moves out
        datasets.append(
            {
                "path": str(i),
                "coordinateTransformations": [
                    {"type": "scale", "scale": [1.0, size, size]},
                    {"type": "translation", "translation": [0.0, y0 + shift, x0 + shift]},
                ],
            }
        )
        log(f"  level {i}: {level.shape[1]} × {level.shape[2]} px, {size:.4g} µm/px")
        if max(level.shape[1:]) <= chunk:
            break
        level = halve(level)
    root.attrs["ome"] = {
        "version": "0.5",
        "multiscales": [
            {
                "name": out.stem,
                "axes": [
                    {"name": "c", "type": "channel"},
                    {"name": "y", "type": "space", "unit": "micrometer"},
                    {"name": "x", "type": "space", "unit": "micrometer"},
                ],
                "datasets": datasets,
            }
        ],
        "omero": {"channels": [{"label": n} for n in names]},
    }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("image", type=Path, help="a TIFF, or any file tifffile reads")
    ap.add_argument("out", type=Path, help="the OME-Zarr folder to write (ending in .zarr)")
    ap.add_argument("--um-per-pixel", type=float, required=True, help="pixel size, in the cells' units")
    ap.add_argument("--origin", type=float, nargs=2, default=(0.0, 0.0), metavar=("X0", "Y0"),
                    help="where the centre of pixel (0, 0) sits among the cells (default 0 0)")
    ap.add_argument("--channels", help="channel names, comma separated (default: channel 0, 1, …)")
    ap.add_argument("--chunk", type=int, default=1024)
    args = ap.parse_args(argv)
    try:
        import tifffile
    except ImportError:
        raise SystemExit("pip install tifffile")
    image = tifffile.imread(args.image)
    print(f"{args.image}: {image.shape} {image.dtype}")
    write_ome_zarr(
        image,
        args.out,
        um_per_pixel=args.um_per_pixel,
        origin=tuple(args.origin),
        channels=args.channels.split(",") if args.channels else None,
        chunk=args.chunk,
    )
    print(f"Wrote {args.out}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
