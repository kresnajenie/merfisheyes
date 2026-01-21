#!/usr/bin/env python3
"""
Single Molecule Data to S3-Compatible Format Processor
Converts .parquet or .csv single molecule data to the binary format expected by SingleMoleculeDataset.ts

Usage:
    python process_single_molecule.py input.parquet output_folder/
    python process_single_molecule.py input.csv output_folder/ --dataset-type xenium

Options:
    --dataset-type    Dataset type for column mappings: xenium, merscope, custom (default: merscope)
    --gene-col        Custom gene column name (overrides dataset-type)
    --x-col           Custom x coordinate column name (overrides dataset-type)
    --y-col           Custom y coordinate column name (overrides dataset-type)
    --z-col           Custom z coordinate column name (overrides dataset-type, optional for 2D)
"""

import argparse
import gzip
import json
import os
import re
import struct
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# Column mappings from lib/config/moleculeColumnMappings.ts
COLUMN_MAPPINGS = {
    "xenium": {
        "gene": "feature_name",
        "x": "x_location",
        "y": "y_location",
        "z": "z_location",
    },
    "merscope": {
        "gene": "gene",
        "x": "global_x",
        "y": "global_y",
        "z": "global_z",
    },
    "custom": {
        "gene": "feature_name",
        "x": "x_location",
        "y": "y_location",
        "z": "z_location",
    },
}


def should_filter_gene(gene: str) -> bool:
    """
    Filter control probes and unassigned genes
    Mirrors lib/utils/gene-filters.ts shouldFilterGene()
    """
    g = gene.lower()
    patterns = [
        r"negative[\s_-]*control",
        r"neg[\s_-]*ctrl",
        r"unassigned",
        r"deprecated",
        r"codeword",
        r"blank",
        r"negcontrol",
    ]

    for pattern in patterns:
        if re.search(pattern, g):
            return True
    return False


def sanitize_gene_name(gene_name: str) -> str:
    """
    Sanitize gene name for filename
    Mirrors lib/utils/SingleMoleculeProcessor.ts sanitizeGeneName()
    """
    # Replace non-alphanumeric with underscore
    sanitized = re.sub(r'[^a-zA-Z0-9]', '_', gene_name)
    # Replace multiple underscores with single
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized


def normalize_coordinates(coordinates: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray]:
    """
    Normalize coordinates to [-1, 1] range
    Mirrors lib/utils/coordinates.ts normalizeCoordinates()

    Returns:
        normalized: Normalized coordinates
        scaling_factor: The max absolute value used for scaling
        center: The center point (mean of each dimension)
    """
    if len(coordinates) == 0:
        return coordinates, 1.0, np.array([0.0, 0.0, 0.0])

    # Calculate center (mean of each dimension)
    center = np.mean(coordinates, axis=0)

    # Center the coordinates
    centered = coordinates - center

    # Find max absolute value
    max_abs = np.max(np.abs(centered))

    if max_abs == 0:
        return centered, 1.0, center

    # Scale to [-1, 1]
    normalized = centered / max_abs

    return normalized, max_abs, center


def read_parquet_file(
    file_path: str,
    gene_col: str,
    x_col: str,
    y_col: str,
    z_col: Optional[str],
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Read parquet file and extract gene names and coordinates

    Returns:
        genes: Array of gene names
        coords: Array of coordinates (N x dimensions)
        dimensions: 2 or 3
    """
    print(f"Reading parquet file: {file_path}")

    # Read parquet file
    table = pq.read_table(file_path)
    df = table.to_pandas()

    # Verify required columns exist
    available_cols = set(df.columns)
    required_cols = {gene_col, x_col, y_col}
    if z_col:
        required_cols.add(z_col)

    missing_cols = required_cols - available_cols
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {sorted(available_cols)}\n"
            f"Expected columns: gene='{gene_col}', x='{x_col}', y='{y_col}'"
            + (f", z='{z_col}'" if z_col else " (2D dataset)")
        )

    # Extract data
    genes = df[gene_col].values
    x = df[x_col].values
    y = df[y_col].values

    # Check if 2D or 3D
    if z_col and z_col in df.columns:
        z = df[z_col].values
        coords = np.column_stack([x, y, z])
        dimensions = 3
    else:
        z = np.zeros_like(x)
        coords = np.column_stack([x, y, z])
        dimensions = 2

    print(f"  Total molecules: {len(genes):,}")
    print(f"  Dimensions: {dimensions}D")

    return genes, coords, dimensions


def read_csv_file(
    file_path: str,
    gene_col: str,
    x_col: str,
    y_col: str,
    z_col: Optional[str],
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Read CSV file and extract gene names and coordinates

    Returns:
        genes: Array of gene names
        coords: Array of coordinates (N x dimensions)
        dimensions: 2 or 3
    """
    print(f"Reading CSV file: {file_path}")

    # Read CSV file
    df = pd.read_csv(file_path)

    # Verify required columns exist
    available_cols = set(df.columns)
    required_cols = {gene_col, x_col, y_col}
    if z_col:
        required_cols.add(z_col)

    missing_cols = required_cols - available_cols
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {sorted(available_cols)}\n"
            f"Expected columns: gene='{gene_col}', x='{x_col}', y='{y_col}'"
            + (f", z='{z_col}'" if z_col else " (2D dataset)")
        )

    # Extract data
    genes = df[gene_col].values
    x = df[x_col].values
    y = df[y_col].values

    # Check if 2D or 3D
    if z_col and z_col in df.columns:
        z = df[z_col].values
        coords = np.column_stack([x, y, z])
        dimensions = 3
    else:
        z = np.zeros_like(x)
        coords = np.column_stack([x, y, z])
        dimensions = 2

    print(f"  Total molecules: {len(genes):,}")
    print(f"  Dimensions: {dimensions}D")

    return genes, coords, dimensions


def process_single_molecule_data(
    input_file: str,
    output_folder: str,
    dataset_type: str = "merscope",
    gene_col: Optional[str] = None,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    z_col: Optional[str] = None,
):
    """
    Process single molecule data and create S3-compatible folder structure
    """
    start_time = time.time()

    # Determine column names
    if dataset_type not in COLUMN_MAPPINGS:
        raise ValueError(
            f"Invalid dataset type: {dataset_type}\n"
            f"Valid types: {', '.join(COLUMN_MAPPINGS.keys())}"
        )

    mapping = COLUMN_MAPPINGS[dataset_type]
    gene_col = gene_col or mapping["gene"]
    x_col = x_col or mapping["x"]
    y_col = y_col or mapping["y"]
    z_col = z_col or mapping["z"]

    print(f"Processing single molecule data...")
    print(f"  Dataset type: {dataset_type}")
    print(f"  Gene column: {gene_col}")
    print(f"  X column: {x_col}")
    print(f"  Y column: {y_col}")
    print(f"  Z column: {z_col}")
    print()

    # Step 1: Read input file (10%)
    progress = 10
    print(f"[{progress:3d}%] Reading input file...")

    file_ext = Path(input_file).suffix.lower()
    if file_ext == ".parquet":
        genes, coords, dimensions = read_parquet_file(
            input_file, gene_col, x_col, y_col, z_col
        )
    elif file_ext == ".csv":
        genes, coords, dimensions = read_csv_file(
            input_file, gene_col, x_col, y_col, z_col
        )
    else:
        raise ValueError(f"Unsupported file type: {file_ext} (expected .parquet or .csv)")

    total_molecules = len(genes)

    # Step 2: Normalize coordinates (30%)
    progress = 30
    print(f"[{progress:3d}%] Normalizing coordinates...")

    normalized_coords, scaling_factor, center = normalize_coordinates(coords)

    print(f"  Scaling factor: {scaling_factor:.2f}")
    print(f"  Center: {center}")

    # Step 3: Build gene index (50-70%)
    progress = 50
    print(f"[{progress:3d}%] Building gene index...")

    gene_index: Dict[str, List[int]] = {}

    progress_interval = max(1, total_molecules // 20)  # Report every 5%

    for i, gene in enumerate(genes):
        if gene not in gene_index:
            gene_index[gene] = []
        gene_index[gene].append(i)

        # Progress update every 5%
        if i > 0 and i % progress_interval == 0:
            progress = 50 + int((i / total_molecules) * 20)
            elapsed = time.time() - start_time
            print(f"[{progress:3d}%] Indexing molecules... ({i:,}/{total_molecules:,}) [{elapsed:.1f}s]")

    progress = 70
    elapsed = time.time() - start_time
    print(f"[{progress:3d}%] Gene index built. Unique genes: {len(gene_index):,} [{elapsed:.1f}s]")

    # Step 4: Filter genes (85%)
    progress = 85
    print(f"[{progress:3d}%] Filtering control probes and unassigned genes...")

    filtered_genes = {gene for gene in gene_index.keys() if not should_filter_gene(gene)}
    filtered_count = len(gene_index) - len(filtered_genes)

    print(f"  Filtered out {filtered_count} genes")
    print(f"  Remaining genes: {len(filtered_genes):,}")

    # Step 5: Create output folder structure (90%)
    progress = 90
    print(f"[{progress:3d}%] Creating output folder structure...")

    output_path = Path(output_folder)
    genes_folder = output_path / "genes"
    genes_folder.mkdir(parents=True, exist_ok=True)

    # Step 6: Write gene files (90-98%)
    progress = 90
    print(f"[{progress:3d}%] Writing gene files...")

    unique_genes = sorted(filtered_genes)
    total_genes = len(unique_genes)

    for idx, gene in enumerate(unique_genes):
        indices = gene_index[gene]

        # Extract coordinates for this gene (already normalized)
        gene_coords = normalized_coords[indices].flatten()  # Flat array: [x1,y1,z1, x2,y2,z2, ...]

        # Convert to Float32
        gene_coords_float32 = gene_coords.astype(np.float32)

        # Sanitize gene name for filename
        sanitized_name = sanitize_gene_name(gene)

        # Write gzipped binary file
        gene_file = genes_folder / f"{sanitized_name}.bin.gz"
        with gzip.open(gene_file, 'wb') as f:
            f.write(gene_coords_float32.tobytes())

        # Progress update
        if idx > 0 and idx % max(1, total_genes // 20) == 0:
            progress = 90 + int((idx / total_genes) * 8)
            elapsed = time.time() - start_time
            print(f"[{progress:3d}%] Writing gene files... ({idx:,}/{total_genes:,}) [{elapsed:.1f}s]")

    # Step 7: Write manifest (98%)
    progress = 98
    elapsed = time.time() - start_time
    print(f"[{progress:3d}%] Writing manifest... [{elapsed:.1f}s]")

    manifest = {
        "uniqueGenes": unique_genes,
        "dimensions": dimensions,
        "scalingFactor": float(scaling_factor),
        "totalMolecules": total_molecules,
        "geneCount": len(unique_genes),
    }

    manifest_json = json.dumps(manifest, indent=2)
    manifest_file = output_path / "manifest.json.gz"

    with gzip.open(manifest_file, 'wt', encoding='utf-8') as f:
        f.write(manifest_json)

    # Done!
    progress = 100
    elapsed = time.time() - start_time
    elapsed_str = f"{elapsed:.2f}s" if elapsed < 60 else f"{int(elapsed // 60)}m {elapsed % 60:.1f}s"

    print(f"[{progress:3d}%] Processing complete! [{elapsed_str}]")
    print()
    print("Output structure:")
    print(f"  {output_folder}/")
    print(f"    manifest.json.gz")
    print(f"    genes/")
    print(f"      {sanitize_gene_name(unique_genes[0])}.bin.gz")
    print(f"      ...")
    print(f"      ({len(unique_genes)} total gene files)")
    print()
    print("Statistics:")
    print(f"  Total molecules: {total_molecules:,}")
    print(f"  Unique genes: {len(unique_genes):,}")
    print(f"  Dimensions: {dimensions}D")
    print(f"  Scaling factor: {scaling_factor:.2f}")
    print()
    print("Ready for S3 upload! Upload the entire output folder to:")
    print(f"  s3://your-bucket/datasets/{{datasetId}}/")


def main():
    parser = argparse.ArgumentParser(
        description="Process single molecule data to S3-compatible format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input_file",
        help="Input parquet or CSV file",
    )
    parser.add_argument(
        "output_folder",
        help="Output folder for processed files",
    )
    parser.add_argument(
        "--dataset-type",
        choices=["xenium", "merscope", "custom"],
        default="merscope",
        help="Dataset type for column mappings (default: merscope)",
    )
    parser.add_argument(
        "--gene-col",
        help="Custom gene column name (overrides dataset-type)",
    )
    parser.add_argument(
        "--x-col",
        help="Custom x coordinate column name (overrides dataset-type)",
    )
    parser.add_argument(
        "--y-col",
        help="Custom y coordinate column name (overrides dataset-type)",
    )
    parser.add_argument(
        "--z-col",
        help="Custom z coordinate column name (overrides dataset-type, optional for 2D)",
    )

    args = parser.parse_args()

    # Verify input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    try:
        process_single_molecule_data(
            input_file=args.input_file,
            output_folder=args.output_folder,
            dataset_type=args.dataset_type,
            gene_col=args.gene_col,
            x_col=args.x_col,
            y_col=args.y_col,
            z_col=args.z_col,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
