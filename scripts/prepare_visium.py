#!/usr/bin/env python3
"""Turn a Visium or Visium HD Space Ranger output into a folder MERFISHEYES opens as it is.

  python scripts/prepare_visium.py <spaceranger_outs> <folder> [--name NAME]
                                   [--fullres IMAGE] [--bins 8 16] [--cluster] [--url BASE]

Reads the output with spatialdata-io (`spatialdata_io.visium`, or
`spatialdata_io.visium_hd` when the folder holds `binned_outputs/`) and
writes

  <folder>/<name>.zarr/   the SpatialData store (Zarr v3): the spots (Visium) or
                          the bins (Visium HD, --bins, 8 and 16 µm by default; 2 µm
                          is too many for a browser) as shapes, each table with a
                          CSC copy of X, the hires and lowres images and, with
                          --fullres, the full-resolution microscope image in 1024²
                          chunks with a pyramid
  <folder>/mapping.json   names the store

then runs scripts/prepare_folder.py's checks and scripts/validate_store.py.
A Visium store has no shared `global` system: the spots and each image
declare their own, and the viewer picks the spots' (docs/specs/store-layout.md).
The spots and bins are drawn at their real size as the cells; there are no
molecules. For Visium HD the 8 µm table is the one the viewer opens
(`metadata_general["Default Table"]`, as SpatialData-DB writes it).

Needs: pip install -r scripts/requirements-prepare.txt spatialdata-io
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_folder  # noqa: E402

CHUNK = 1024
IMAGE_KWARGS = {"chunks": (1, CHUNK, CHUNK), "scale_factors": [2, 2, 2, 2]}


def read_visium(path: Path, name: str, *, fullres: Path | None, bins: list[int], log=print):
    try:
        import spatialdata_io
    except ImportError:
        raise SystemExit("pip install spatialdata-io   (Python 3.12 or newer)")

    hd = (path / "binned_outputs").is_dir()
    log(f"reading {path} as {'Visium HD' if hd else 'Visium'} with spatialdata-io ...")
    if hd:
        sdata = spatialdata_io.visium_hd(
            path,
            dataset_id=name,
            bin_size=bins,
            bins_as_squares=True,
            fullres_image_file=fullres,
            image_models_kwargs=IMAGE_KWARGS,
        )
    else:
        sdata = spatialdata_io.visium(
            path,
            dataset_id=name,
            fullres_image_file=fullres,
            image_models_kwargs=IMAGE_KWARGS,
        )
    return sdata, hd


def set_default_table(store: Path, table: str) -> None:
    """`metadata_general["Default Table"]` in the root attributes (the
    viewer's rule for which table it opens; the same key SpatialData-DB uses)."""
    import zarr

    prepare_folder.drop_consolidated(store / "zarr.json")
    root = zarr.open_group(store, mode="r+", use_consolidated=False)
    general = dict(root.attrs.get("metadata_general") or {})
    general["Default Table"] = table
    root.attrs["metadata_general"] = general


def prepare_visium(
    path: Path,
    folder: Path,
    *,
    name: str | None = None,
    fullres: Path | None = None,
    bins: list[int] | None = None,
    cluster: bool = False,
    url: str | None = None,
    images: bool = True,
    label_csv: Path | list[Path] | None = None,
    log=print,
):
    path, folder = Path(path), Path(folder)
    name = name or path.name.removesuffix("_outs")
    bins = bins or [8, 16]
    folder.mkdir(parents=True, exist_ok=True)
    store = folder / f"{name}.zarr"
    if store.exists():
        raise SystemExit(f"{store} exists already; remove it or pass --name")
    if fullres is not None and not Path(fullres).is_absolute() and not (path / fullres).exists():
        raise SystemExit(f"--fullres {fullres}: not found in {path}")

    sdata, hd = read_visium(path, name, fullres=fullres, bins=bins, log=log)
    if not images:
        for key in list(sdata.images):
            del sdata.images[key]
    log(f"  {sdata}")
    log(f"writing {store} ...")
    sdata.write(store)
    if hd:
        default = f"square_{min(bins):03d}um"
        set_default_table(store, default)
        log(f"  default table: {default}")

    log("checking the folder ...")
    return prepare_folder.prepare(folder, cells=store, cluster=cluster, url=url, labels=label_csv, log=log)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("outs", type=Path, help="the Space Ranger output folder (holds spatial/ and the counts h5)")
    ap.add_argument("folder", type=Path, help="the folder to write; it is what you drop, save or sync")
    ap.add_argument("--name", help="the store's name (default: the output folder's, without _outs)")
    ap.add_argument("--fullres", type=Path, help="the full-resolution microscope image (a file name inside the output folder, or a path)")
    ap.add_argument("--bins", type=int, nargs="+", default=[8, 16], help="Visium HD bin sizes to read, in µm (default 8 16)")
    ap.add_argument("--cluster", action="store_true", help="Leiden clusters + DEGs when the table has no cell types")
    ap.add_argument("--url", help="the folder's public URL, for absolute links in mapping.json")
    args = ap.parse_args(argv)
    r = prepare_visium(
        args.outs, args.folder, name=args.name, fullres=args.fullres, bins=args.bins, cluster=args.cluster, url=args.url
    )
    print(json.dumps({"failed": r.failed}))
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
