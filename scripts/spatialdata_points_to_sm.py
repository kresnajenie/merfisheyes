#!/usr/bin/env python3
"""
SpatialData points → single molecule format

Reads the transcripts of a SpatialData zarr store (points/<element>) and
writes the per-gene format the single molecule viewer loads
(manifest.json.gz + genes/{gene}.bin.gz). The store itself is only read.

Usage:
    python spatialdata_points_to_sm.py store.zarr output_folder/
    python spatialdata_points_to_sm.py store.zarr output_folder/ --points transcripts --workers 16
    python spatialdata_points_to_sm.py store.zarr output_folder/ --drop-z

Coordinates are written as stored in the points element (the same units as
the table's obsm/spatial), so molecules overlay the store's cells directly.
"""

import argparse
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

import process_single_molecule as psm


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("store", help="SpatialData zarr store")
    parser.add_argument("output_folder")
    parser.add_argument("--points", help="points element (default: the only one)")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--drop-z",
        action="store_true",
        help="write 2D molecules, flat on the plane of the cells and the image",
    )
    args = parser.parse_args()

    points_dir = Path(args.store) / "points"
    elements = sorted(p.name for p in points_dir.iterdir() if p.is_dir())
    element = args.points or (elements[0] if len(elements) == 1 else None)
    if element not in elements:
        sys.exit(f"Pick a points element with --points; found: {elements}")

    attrs = json.loads((points_dir / element / "zarr.json").read_text())["attributes"]
    feature_key = attrs["spatialdata_attrs"]["feature_key"]
    axes = attrs["axes"]

    # The element's parquet is a directory of part files; read_table handles
    # that, read_schema doesn't.
    psm.pq.read_schema = lambda path: pq.ParquetDataset(path).schema

    psm.process_single_molecule_data(
        input_file=str(points_dir / element / "points.parquet"),
        output_folder=args.output_folder,
        dataset_type="xenium",
        gene_col=feature_key,
        x_col="x",
        y_col="y",
        z_col="z" if "z" in axes and not args.drop_z else None,
        num_workers=args.workers,
    )


if __name__ == "__main__":
    main()
