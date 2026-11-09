#!/usr/bin/env python3
"""Prepare spatial data for MERFISHEYES: whatever you have becomes one SpatialData store.

  python scripts/prepare.py INPUT OUT [--name NAME] [--no-images] [--keep-z]
                            [--cluster] [--labels CSV]...
                            [--molecule-columns GENE,X,Y[,Z] | gene=G,x=X,y=Y[,z=Z][,cell=C]] [--workers N]

Run it on your own machine to prepare data before you drop or upload it. The
server's worker is moving to this same script (docs/plans/ONE-PIPELINE.md,
step 2); once it has, "Upload & process on server" and a local run give the
same store. Install what it needs once:

  python3.12 -m venv .venv && source .venv/bin/activate
  pip install -r scripts/requirements.txt

INPUT is recognised by what is in it, never by its name:

  an AnnData table     a .h5ad file, or an AnnData Zarr folder (v2 or v3)
  a Xenium export      the folder holding experiment.xenium
  a MERSCOPE region    the folder holding cell_by_gene.csv and cell_metadata.csv
  a Visium output      Space Ranger's outs/ (spatial/ and a feature matrix)
  molecules only       a .parquet / .csv with a gene column and x, y (, z)
                       columns: Xenium's and MERSCOPE's names, or your own
                       (Gene, X, Y, Z... in any case; else --molecule-columns)
  a folder of parts    a table, a molecule table and an OME-Zarr image side
                       by side (what scripts/prepare_folder.py takes)

OUT is written and nothing in INPUT is changed. OUT holds exactly:

  NAME.zarr/      the SpatialData store: tables/table (X, with a CSC copy in
                  layers/X_csc unless X is CSC; obsm/spatial; control probes
                  left out), points/molecules (tiled, when there are molecules;
                  flat on the cells' plane unless --keep-z keeps their z),
                  shapes/cell_boundaries (when the input has polygons),
                  images/... (when the input has an image; --no-images skips it);
                  the DEG panel's stats for every categorical column in the
                  table's varm / uns (scripts/de_stats.py)
  mapping.json    {"cells": "NAME.zarr"}
  stats.json      what was written: cells, genes, molecules, the input kind, and
                  the pipeline's version and packages

Drag OUT onto the homepage to look at it, then Upload & Save to keep it.

Exit codes: 0 written and validated; 2 INPUT not recognised; 3 INPUT
incomplete or unusable (the message names what is missing); 4 the store was
written but scripts/validate_store.py failed it.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import prepare_folder  # noqa: E402

PACKAGES = ["anndata", "numpy", "pandas", "pyarrow", "scipy", "spatialdata", "spatialdata-io", "zarr"]


class Unrecognised(Exception):
    """INPUT is none of the kinds this script knows."""


# ── what INPUT is ────────────────────────────────────────────────────────

def recognise(path: Path, molecule_columns: str | None = None) -> str:
    """`h5ad`, `table` (AnnData Zarr), `xenium`, `merscope`, `visium`,
    `transcripts` or `folder`."""
    if path.is_file():
        kind = prepare_folder.kind_of(path, molecule_columns)
        if kind == "h5ad":
            return "h5ad"
        if kind == "transcripts":
            return "transcripts"
        raise Unrecognised(
            f"{path.name}: neither an .h5ad nor a molecule table (a .parquet / .csv "
            "with a gene column and x, y columns; name them with --molecule-columns)"
        )
    if not path.is_dir():
        raise Unrecognised(f"{path}: no such file or folder")
    if (path / "experiment.xenium").exists():
        return "xenium"
    if (path / "cell_by_gene.csv").exists() and (path / "cell_metadata.csv").exists():
        return "merscope"
    if (path / "spatial").is_dir() and (
        any(path.glob("*feature_bc_matrix*")) or (path / "binned_outputs").is_dir()
    ):
        return "visium"
    kind = prepare_folder.kind_of(path, molecule_columns)
    if kind in ("table3", "table2"):
        return "table"
    parts = prepare_folder.find_parts(path, molecule_columns)
    if kind == "store" or any(parts.values()):
        return "folder"
    raise Unrecognised(
        f"{path}: not an AnnData table, a Xenium export, a MERSCOPE region, a "
        "Visium output, a molecule table nor a folder holding any of them"
    )


# ── writing OUT ──────────────────────────────────────────────────────────

def stage(source: Path, folder: Path) -> Path:
    """`source` inside `folder` for prepare_folder, which works in a folder
    of parts: read-only parts are linked; a Zarr table or store, which it
    edits or moves, is copied, so INPUT is never changed."""
    target = folder / source.name
    edited = source.is_dir() and prepare_folder.kind_of(source) in ("store", "table3")
    if edited:
        shutil.copytree(source, target)
    else:
        target.symlink_to(source.resolve(), target_is_directory=source.is_dir())
    return target


def write_transcripts(source: Path, out: Path, name: str, columns: str | None, log, keep_z: bool = False) -> None:
    """Molecules only: a store holding the tiled points element."""
    import transcripts_to_points

    store = out / f"{name}.zarr"
    transcripts_to_points.convert(source, store, columns=columns, keep_z=keep_z, log=log)
    (out / "mapping.json").write_text(json.dumps({"cells": store.name}, indent=2) + "\n")


def mapped_store(out: Path) -> Path | None:
    """The store mapping.json names, when it does."""
    mapping = prepare_folder.read_json(out / "mapping.json") or {}
    cells = mapping.get("cells")
    return out / cells if cells else None


def pipeline_version() -> str:
    """The pipeline's version: MERFISHEYES_PIPELINE_VERSION (the worker image
    sets it), else the git commit of these scripts."""
    if os.environ.get("MERFISHEYES_PIPELINE_VERSION"):
        return os.environ["MERFISHEYES_PIPELINE_VERSION"]
    try:
        return subprocess.run(
            ["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001 — not a checkout
        return "unknown"


def package_versions() -> dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    found = {}
    for name in PACKAGES:
        try:
            found[name] = version(name)
        except PackageNotFoundError:
            pass
    return found


def store_stats(store: Path) -> dict:
    """What the store holds: cells and genes of its table, molecules and
    genes of its points element, its images."""
    import zarr

    root = zarr.open_group(store, mode="r")
    stats: dict = {"cells": 0, "genes": 0, "molecules": 0, "molecule_genes": 0, "dims": 2, "images": []}
    if "tables" in root:
        table = root["tables"][prepare_folder.default_table_name(store)]
        x = table["X"]
        shape = x.attrs.get("shape") or list(x.shape)
        stats["cells"], stats["genes"] = int(shape[0]), int(shape[1])
        spatial = table["obsm"]["spatial"]
        stats["dims"] = 3 if (spatial.shape[1] if hasattr(spatial, "shape") else 2) >= 3 else 2
        uns = table["uns"] if "uns" in table else {}
        stats["de_stats"] = sorted(uns["merfisheyes_de"]) if "merfisheyes_de" in uns else []
    if "points" in root:
        import points_zarr

        for name in root["points"]:
            attrs = root["points"][name].attrs
            if not attrs.get("points:version"):
                continue
            stats["molecules"] += int(attrs.get("points:n_points") or 0)
            # The genes: the categories of the element's feature column
            frame = root["points"][name]["0"]["feature_major"]
            feature = points_zarr._col(frame, attrs["spatialdata_attrs"]["feature_key"])
            stats["molecule_genes"] = max(stats["molecule_genes"], int(feature["enum_values"].shape[0]))
            if not stats["cells"]:
                stats["dims"] = len(attrs.get("axes") or ["x", "y"])
    if "images" in root:
        stats["images"] = sorted(root["images"])
    return stats


def prepare(
    source: Path,
    out: Path,
    *,
    name: str | None = None,
    images: bool = True,
    cluster: bool = False,
    labels: Path | list[Path] | None = None,
    molecule_columns: str | None = None,
    workers: int = 4,
    keep_z: bool = False,
    log=print,
) -> int:
    source = Path(source).resolve()
    out = Path(out).resolve()
    try:
        kind = recognise(source, molecule_columns)
    except Unrecognised as e:
        log(f"error: {e}")
        return 2
    if out.exists() and any(out.iterdir()):
        log(f"error: {out} is not empty; give a new folder")
        return 2
    out.mkdir(parents=True, exist_ok=True)
    name = name or source.stem.removesuffix("_outs")
    log(f"{source.name}: {kind} → {out}/{name}.zarr")

    label_csv = [Path(c).resolve() for c in ([labels] if isinstance(labels, (str, Path)) else labels or [])]
    try:
        if kind == "xenium":
            import prepare_xenium

            report = prepare_xenium.prepare_xenium(
                source, out, name=name, images=images, workers=workers, cluster=cluster,
                label_csv=label_csv, keep_z=keep_z, log=log,
            )
        elif kind == "merscope":
            import prepare_merscope

            report = prepare_merscope.prepare_merscope(
                source, out, name=name, workers=workers, cluster=cluster, label_csv=label_csv,
                keep_z=keep_z, log=log,
            )
        elif kind == "visium":
            import prepare_visium

            report = prepare_visium.prepare_visium(
                source, out, name=name, cluster=cluster, images=images, label_csv=label_csv, log=log,
            )
        elif kind == "transcripts":
            write_transcripts(source, out, name, molecule_columns, log, keep_z=keep_z)
            import validate_store

            report = validate_store.validate(out / f"{name}.zarr", out / "mapping.json")
            for level, where, msg in report.lines:
                log(f"{level:4} {where}: {msg}")
        else:
            # A table or a folder of parts, staged in OUT and prepared there
            parts = [source] if kind in ("h5ad", "table") else [
                p for p in sorted(source.iterdir()) if not p.name.startswith(".")
            ]
            for part in parts:
                stage(part, out)
            report = prepare_folder.prepare(
                out, name=name, cluster=cluster, molecule_columns=molecule_columns,
                workers=workers, images=images, labels=label_csv, keep_z=keep_z, log=log,
            )
            # What OUT was made from stays in INPUT; a store that was there
            # already is OUT's store (the one mapping.json names)
            kept = mapped_store(out)
            for staged in [out / p.name for p in parts] + [out / "sources"]:
                if staged == kept:
                    continue
                if staged.is_symlink() or staged.is_file():
                    staged.unlink()
                elif staged.is_dir():
                    shutil.rmtree(staged)
    except SystemExit as e:
        log(f"error: {e}")
        return 3

    store = mapped_store(out) or out / f"{name}.zarr"
    if not store.exists():
        # Nothing became a store: a merfisheyes molecule folder or an image
        # with no table to go into is not read on its own yet. Exiting 0 here
        # left an OUT whose mapping.json named nothing, which the server took
        # for a finished dataset.
        log(
            f"error: no store was written ({store.name} is missing): a merfisheyes "
            "molecule folder or an image is only read beside a table for now"
        )
        return 3
    stats = {
        "input": kind,
        "store": store.name,
        **(store_stats(store) if store.exists() else {}),
        "pipeline": pipeline_version(),
        "packages": package_versions(),
    }
    (out / "stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    log(f"stats.json: {json.dumps({k: v for k, v in stats.items() if k != 'packages'})}")
    return 4 if report.failed else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="Exit codes: 0 done; 2 INPUT not recognised; 3 INPUT incomplete; 4 the store failed validation.",
    )
    ap.add_argument("input", type=Path, help="what you have: a file or a folder (see above)")
    ap.add_argument("out", type=Path, help="a new folder to write; it is what you drop or upload")
    ap.add_argument("--name", help="the store's name (default: INPUT's)")
    ap.add_argument("--no-images", action="store_true", help="leave an image that INPUT has out of the store")
    ap.add_argument("--keep-z", action="store_true", help="keep the molecules' z (3-D); by default they are flat")
    ap.add_argument("--cluster", action="store_true", help="Leiden clusters + DEGs when the table has no cell types")
    ap.add_argument("--labels", type=Path, action="append",
                    help="per-cell labels to add to the table (a MapMyCells CSV, an artifact mask); "
                         "repeat for several, merged in the order given")
    ap.add_argument("--molecule-columns", help="a molecule table's columns as GENE,X,Y[,Z] when they are not recognised, or gene=G,x=X,y=Y[,z=Z][,cell=C] to also split assigned and unassigned molecules by C")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)
    return prepare(
        args.input,
        args.out,
        name=args.name,
        images=not args.no_images,
        cluster=args.cluster,
        labels=args.labels,
        molecule_columns=args.molecule_columns,
        workers=args.workers,
        keep_z=args.keep_z,
    )


if __name__ == "__main__":
    sys.exit(main())
