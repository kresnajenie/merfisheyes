#!/usr/bin/env python3
"""
Single Molecule Data to S3-Compatible Format Processor
Converts .parquet or .csv single molecule data to the binary format expected by SingleMoleculeDataset.ts

Usage:
    # Single file mode
    python process_single_molecule.py input.parquet output_folder/
    python process_single_molecule.py input.csv output_folder/ --dataset-type xenium
    python process_single_molecule.py input.parquet output_folder/ --manifest-only

    # Directory mode (processes all detected_transcripts.csv in subdirectories)
    python process_single_molecule.py /path/to/dataset_folder/ output_folder/ --s3-prefix https://bucket.s3.region.amazonaws.com/prefix

Options:
    --dataset-type    Dataset type for column mappings: xenium, merscope, custom (default: merscope)
    --gene-col        Custom gene column name (overrides dataset-type)
    --x-col           Custom x coordinate column name (overrides dataset-type)
    --y-col           Custom y coordinate column name (overrides dataset-type)
    --z-col           Custom z coordinate column name (overrides dataset-type, optional for 2D)
    --manifest-only   Only generate manifest.json.gz without creating gene files (faster)
    --cell-id-col     Column name for cell assignment (auto-detected for MERSCOPE)
                      Molecules with value -1 are treated as unassigned
    --s3-prefix       S3 base URL prefix for mapping file (directory mode only)
    --link-column     Cluster column name to read for linking (default: _sample_id)
"""

import argparse
import datetime
import gzip
import json
import os
import re
import struct
import sys
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# ─────────────────────────────────────────────
# FUZZY FILE MATCHING (mirrors combine_slices_v3.py)
# ─────────────────────────────────────────────

TRANSCRIPT_PATTERN = {
    'keywords': ['detected', 'transcript'],
    'alternative_keywords': [['transcript']],
    'exclude': ['metadata', 'cell_by_gene'],
    'extensions': ['.csv'],
    'description': 'Detected transcripts',
}


def normalize_filename(filename: str) -> str:
    """Lowercase, strip extension, replace separators with spaces."""
    name = Path(filename).stem.lower()
    return re.sub(r'[-_.]', ' ', name)


def match_transcript_file(filename: str) -> Tuple[int, str]:
    """
    Returns (score, reason). score > 0 means a match for detected_transcripts.
    """
    normalized = normalize_filename(filename)
    ext = Path(filename).suffix.lower()

    if ext not in TRANSCRIPT_PATTERN['extensions']:
        return 0, f"Extension {ext} not in {TRANSCRIPT_PATTERN['extensions']}"

    for excl in TRANSCRIPT_PATTERN.get('exclude', []):
        if excl in normalized:
            return 0, f"Contains excluded keyword '{excl}'"

    # Primary keywords — all must be present
    keywords = TRANSCRIPT_PATTERN['keywords']
    hits = [kw for kw in keywords
            if kw in normalized or kw + 's' in normalized
            or (kw.endswith('s') and kw[:-1] in normalized)]

    if len(hits) == len(keywords):
        return 100, f"Primary match: all keywords {keywords} found"

    # Alternative keyword combos
    for alt_combo in TRANSCRIPT_PATTERN.get('alternative_keywords', []):
        alt_hits = [kw for kw in alt_combo
                    if kw in normalized or kw + 's' in normalized
                    or (kw.endswith('s') and kw[:-1] in normalized)]
        if len(alt_hits) == len(alt_combo):
            return 80, f"Alternative match: keywords {alt_combo} found"

    return 0, f"No match"


def find_transcript_in_dir(directory: str) -> Optional[Path]:
    """
    Scan a directory for a file matching detected_transcripts using fuzzy matching.
    Returns the path to the best match, or None.
    """
    matches = []
    try:
        for entry in os.scandir(directory):
            if not entry.is_file(follow_symlinks=False):
                continue
            score, reason = match_transcript_file(entry.name)
            if score > 0:
                matches.append((entry.path, score, reason))
    except PermissionError:
        return None

    if not matches:
        return None

    if len(matches) > 1:
        print(f"  WARNING: Multiple transcript files in {directory}:")
        for path, score, reason in matches:
            print(f"    {Path(path).name} (score: {score}, {reason})")
        print(f"  Using highest score match.")

    best = max(matches, key=lambda x: x[1])
    return Path(best[0])


def scan_dir_once(directory: str) -> Tuple[List[str], List[str]]:
    """
    Single os.scandir pass. Returns (csv_names, subdirs).
    """
    csv_names = []
    subdirs = []
    try:
        for entry in os.scandir(directory):
            if entry.is_dir(follow_symlinks=False):
                subdirs.append(entry.path)
            elif entry.is_file(follow_symlinks=False) and entry.name.lower().endswith('.csv'):
                csv_names.append(entry.name)
    except PermissionError:
        pass
    return csv_names, subdirs


def has_transcript_file(csv_names: List[str]) -> bool:
    """Check if a list of CSV filenames contains a detected_transcripts file."""
    for name in csv_names:
        score, _ = match_transcript_file(name)
        if score > 0:
            return True
    return False


def discover_transcript_dirs_bfs(target: Path) -> List[Tuple[str, Path]]:
    """
    BFS from each child of target to find directories containing detected_transcripts.csv.
    Returns list of (sample_id, directory_path) tuples.
    sample_id is the top-level directory name relative to target (matches _sample_id in combine_slices_v3.py).
    """
    results = []

    # Get top-level children
    top_children = []
    for entry in os.scandir(target):
        if entry.is_dir(follow_symlinks=False):
            top_children.append((entry.name, entry.path))
    top_children.sort()

    print(f"Scanning {len(top_children)} top-level directories...")

    for idx, (child_name, child_path) in enumerate(top_children, 1):
        print(f"  [{idx}/{len(top_children)}] Searching '{child_name}'...")

        # BFS within this child branch
        queue = deque([child_path])
        dirs_scanned = 0

        while queue:
            current = queue.popleft()
            csv_names, subdirs = scan_dir_once(current)
            dirs_scanned += 1

            if has_transcript_file(csv_names):
                results.append((child_name, Path(current)))
                print(f"    FOUND at {current} (scanned {dirs_scanned} dirs)")
                break

            subdirs.sort()
            queue.extend(subdirs)
        else:
            print(f"    not found (scanned {dirs_scanned} dirs)")

    return results


# Column mappings from lib/config/moleculeColumnMappings.ts
COLUMN_MAPPINGS = {
    "xenium": {
        "gene": "feature_name",
        "x": "x_location",
        "y": "y_location",
        "z": "z_location",
        "cell_id": None,
    },
    "merscope": {
        "gene": "gene",
        "x": "global_x",
        "y": "global_y",
        "z": "global_z",
        "cell_id": "cell_id",
    },
    "custom": {
        "gene": "feature_name",
        "x": "x_location",
        "y": "y_location",
        "z": "z_location",
        "cell_id": None,
    },
}


UNASSIGNED_SUFFIX = "_uuuuuuuuuu"
UNASSIGNED_CELL_ID = -1


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


def round_coordinates(coordinates: np.ndarray) -> np.ndarray:
    """
    Round coordinates to 2 decimal places.
    Matches browser behavior: Math.round(x * 100) / 100

    Returns:
        rounded: Coordinates rounded to 2dp
    """
    if len(coordinates) == 0:
        return coordinates

    return np.round(coordinates, 2)


def _write_gene_file_worker(args):
    """Worker function for parallel gene file writing. Must be top-level for pickling."""
    gene, sanitized_name, genes_folder, assigned_coords_bytes, unassigned_coords_bytes = args

    # Write assigned molecules file
    gene_file = os.path.join(genes_folder, f"{sanitized_name}.bin.gz")
    with gzip.open(gene_file, 'wb') as f:
        f.write(assigned_coords_bytes)

    # Write unassigned molecules file if present
    if unassigned_coords_bytes is not None:
        unassigned_file = os.path.join(genes_folder, f"{sanitized_name}{UNASSIGNED_SUFFIX}.bin.gz")
        with gzip.open(unassigned_file, 'wb') as f:
            f.write(unassigned_coords_bytes)

    return gene


def _process_sample_worker(args):
    """Worker function for parallel sample processing. Must be top-level for pickling."""
    (sample_id, transcript_file, sample_output,
     dataset_type, gene_col, x_col, y_col, z_col,
     cell_id_col, manifest_only, num_workers) = args

    try:
        process_single_molecule_data(
            input_file=str(transcript_file),
            output_folder=sample_output,
            dataset_type=dataset_type,
            gene_col=gene_col,
            x_col=x_col,
            y_col=y_col,
            z_col=z_col,
            cell_id_col=cell_id_col,
            manifest_only=manifest_only,
            num_workers=num_workers,
        )
        return sample_id, True, None
    except Exception as e:
        return sample_id, False, str(e)


def read_parquet_file(
    file_path: str,
    gene_col: str,
    x_col: str,
    y_col: str,
    z_col: Optional[str],
    cell_id_col: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, int, Optional[np.ndarray]]:
    """
    Read parquet file and extract gene names and coordinates

    Returns:
        genes: Array of gene names
        coords: Array of coordinates (N x dimensions)
        dimensions: 2 or 3
        cell_ids: Array of cell IDs (or None if cell_id_col not provided/found)
    """
    print(f"Reading parquet file: {file_path}")

    # Read parquet file
    table = pq.read_table(file_path)
    df = table.to_pandas()

    # Verify required columns exist (only gene, x, y are required)
    available_cols = set(df.columns)
    required_cols = {gene_col, x_col, y_col}

    missing_cols = required_cols - available_cols
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {sorted(available_cols)}\n"
            f"Expected columns: gene='{gene_col}', x='{x_col}', y='{y_col}'"
            + (f", z='{z_col}' (optional for 3D)" if z_col else "")
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
        print(f"  Detected 3D dataset (z column '{z_col}' found)")
    else:
        z = np.zeros_like(x)
        coords = np.column_stack([x, y, z])
        dimensions = 2
        if z_col:
            print(f"  Detected 2D dataset (z column '{z_col}' not found in file)")
        else:
            print(f"  Detected 2D dataset (no z column specified)")

    # Read cell_id column if provided
    cell_ids = None
    if cell_id_col and cell_id_col in df.columns:
        cell_ids = df[cell_id_col].values
        n_unassigned = np.sum(cell_ids == UNASSIGNED_CELL_ID)
        print(f"  Cell ID column '{cell_id_col}' found: {n_unassigned:,} unassigned ({n_unassigned/len(genes)*100:.1f}%)")
    elif cell_id_col:
        print(f"  Cell ID column '{cell_id_col}' not found in file, skipping assignment split")

    print(f"  Total molecules: {len(genes):,}")
    print(f"  Dimensions: {dimensions}D")

    return genes, coords, dimensions, cell_ids


def read_csv_file(
    file_path: str,
    gene_col: str,
    x_col: str,
    y_col: str,
    z_col: Optional[str],
    cell_id_col: Optional[str] = None,
) -> Tuple[np.ndarray, np.ndarray, int, Optional[np.ndarray]]:
    """
    Read CSV file and extract gene names and coordinates

    Returns:
        genes: Array of gene names
        coords: Array of coordinates (N x dimensions)
        dimensions: 2 or 3
        cell_ids: Array of cell IDs (or None if cell_id_col not provided/found)
    """
    print(f"Reading CSV file: {file_path}")

    # Read CSV file
    df = pd.read_csv(file_path)

    # Verify required columns exist (only gene, x, y are required)
    available_cols = set(df.columns)
    required_cols = {gene_col, x_col, y_col}

    missing_cols = required_cols - available_cols
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}\n"
            f"Available columns: {sorted(available_cols)}\n"
            f"Expected columns: gene='{gene_col}', x='{x_col}', y='{y_col}'"
            + (f", z='{z_col}' (optional for 3D)" if z_col else "")
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
        print(f"  Detected 3D dataset (z column '{z_col}' found)")
    else:
        z = np.zeros_like(x)
        coords = np.column_stack([x, y, z])
        dimensions = 2
        if z_col:
            print(f"  Detected 2D dataset (z column '{z_col}' not found in file)")
        else:
            print(f"  Detected 2D dataset (no z column specified)")

    # Read cell_id column if provided
    cell_ids = None
    if cell_id_col and cell_id_col in df.columns:
        cell_ids = df[cell_id_col].values
        n_unassigned = np.sum(cell_ids == UNASSIGNED_CELL_ID)
        print(f"  Cell ID column '{cell_id_col}' found: {n_unassigned:,} unassigned ({n_unassigned/len(genes)*100:.1f}%)")
    elif cell_id_col:
        print(f"  Cell ID column '{cell_id_col}' not found in file, skipping assignment split")

    print(f"  Total molecules: {len(genes):,}")
    print(f"  Dimensions: {dimensions}D")

    return genes, coords, dimensions, cell_ids


def process_single_molecule_data(
    input_file: str,
    output_folder: str,
    dataset_type: str = "merscope",
    gene_col: Optional[str] = None,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    z_col: Optional[str] = None,
    cell_id_col: Optional[str] = None,
    manifest_only: bool = False,
    num_workers: int = 1,
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
    cell_id_col = cell_id_col or mapping.get("cell_id")

    print(f"Processing single molecule data...")
    print(f"  Dataset type: {dataset_type}")
    print(f"  Gene column: {gene_col}")
    print(f"  X column: {x_col}")
    print(f"  Y column: {y_col}")
    print(f"  Z column: {z_col}")
    if cell_id_col:
        print(f"  Cell ID column: {cell_id_col}")
    print(f"  Workers: {num_workers}")
    print()

    # Step 1: Read input file (10%)
    progress = 10
    print(f"[{progress:3d}%] Reading input file...")

    file_ext = Path(input_file).suffix.lower()
    if file_ext == ".parquet":
        genes, coords, dimensions, cell_ids = read_parquet_file(
            input_file, gene_col, x_col, y_col, z_col, cell_id_col
        )
    elif file_ext == ".csv":
        genes, coords, dimensions, cell_ids = read_csv_file(
            input_file, gene_col, x_col, y_col, z_col, cell_id_col
        )
    else:
        raise ValueError(f"Unsupported file type: {file_ext} (expected .parquet or .csv)")

    total_molecules = len(genes)

    # Step 2: Round coordinates to 2dp (30%)
    progress = 30
    print(f"[{progress:3d}%] Rounding coordinates to 2 decimal places...")

    rounded_coords = round_coordinates(coords)

    print(f"  Coordinate range: x=[{rounded_coords[:,0].min():.2f}, {rounded_coords[:,0].max():.2f}], y=[{rounded_coords[:,1].min():.2f}, {rounded_coords[:,1].max():.2f}]")

    # Step 3: Build gene index using vectorized pandas groupby (50-70%)
    progress = 50
    print(f"[{progress:3d}%] Building gene index (vectorized)...")

    has_unassigned = cell_ids is not None

    if has_unassigned:
        # Split into assigned vs unassigned using boolean mask
        unassigned_mask = cell_ids == UNASSIGNED_CELL_ID
        assigned_mask = ~unassigned_mask

        # Assigned gene index via groupby
        assigned_indices = np.where(assigned_mask)[0]
        assigned_genes = genes[assigned_indices]
        df_assigned = pd.DataFrame({'gene': assigned_genes, 'idx': assigned_indices})
        gene_index = {gene: group['idx'].values for gene, group in df_assigned.groupby('gene')}
        del df_assigned

        # Unassigned gene index via groupby
        unassigned_indices_arr = np.where(unassigned_mask)[0]
        unassigned_genes = genes[unassigned_indices_arr]
        df_unassigned = pd.DataFrame({'gene': unassigned_genes, 'idx': unassigned_indices_arr})
        unassigned_index = {gene: group['idx'].values for gene, group in df_unassigned.groupby('gene')}
        del df_unassigned
    else:
        # All molecules are assigned
        df_all = pd.DataFrame({'gene': genes, 'idx': np.arange(len(genes))})
        gene_index = {gene: group['idx'].values for gene, group in df_all.groupby('gene')}
        del df_all
        unassigned_index = {}

    progress = 70
    elapsed = time.time() - start_time
    print(f"[{progress:3d}%] Gene index built. Unique genes: {len(gene_index):,} [{elapsed:.1f}s]")
    if has_unassigned:
        total_assigned = sum(len(v) for v in gene_index.values())
        total_unassigned = sum(len(v) for v in unassigned_index.values())
        print(f"  Assigned molecules: {total_assigned:,}")
        print(f"  Unassigned molecules: {total_unassigned:,}")
        print(f"  Genes with unassigned molecules: {len(unassigned_index):,}")

    # Step 4: Filter genes (85%)
    progress = 85
    print(f"[{progress:3d}%] Filtering control probes and unassigned genes...")

    # Collect all gene names from both indices
    all_gene_names = set(gene_index.keys()) | set(unassigned_index.keys())
    filtered_genes = {gene for gene in all_gene_names if not should_filter_gene(gene)}
    filtered_count = len(all_gene_names) - len(filtered_genes)

    print(f"  Filtered out {filtered_count} genes")
    print(f"  Remaining genes: {len(filtered_genes):,}")

    # Step 5: Create output folder structure (90%)
    progress = 90
    print(f"[{progress:3d}%] Creating output folder structure...")

    output_path = Path(output_folder)

    unique_genes = sorted(filtered_genes)
    total_genes = len(unique_genes)

    if not manifest_only:
        genes_folder = output_path / "genes"
        genes_folder.mkdir(parents=True, exist_ok=True)

        # Step 6: Write gene files (90-98%)
        # Prepare worker args: pre-serialize coordinates to bytes in main process
        # (numpy slicing is fast, avoids sending huge arrays to workers)
        print(f"[{progress:3d}%] Preparing gene file data...")

        worker_args = []
        for gene in unique_genes:
            sanitized_name = sanitize_gene_name(gene)

            # Assigned coordinates
            if gene in gene_index:
                indices = gene_index[gene]
                assigned_bytes = rounded_coords[indices].flatten().astype(np.float32).tobytes()
            else:
                assigned_bytes = b''

            # Unassigned coordinates
            if has_unassigned and gene in unassigned_index:
                u_indices = unassigned_index[gene]
                unassigned_bytes = rounded_coords[u_indices].flatten().astype(np.float32).tobytes()
            else:
                unassigned_bytes = None

            worker_args.append((gene, sanitized_name, str(genes_folder), assigned_bytes, unassigned_bytes))

        total_files = len(worker_args)

        elapsed = time.time() - start_time
        print(f"[{progress:3d}%] Writing {total_files} gene files (workers={num_workers})... [{elapsed:.1f}s]")

        if num_workers > 1 and total_files > 1:
            effective_workers = min(num_workers, total_files)
            with ProcessPoolExecutor(max_workers=effective_workers) as pool:
                futures = {pool.submit(_write_gene_file_worker, a): a[0] for a in worker_args}
                done_count = 0
                for future in as_completed(futures):
                    future.result()
                    done_count += 1
                    if done_count % max(1, total_files // 10) == 0 or done_count == total_files:
                        progress = 90 + int((done_count / total_files) * 8)
                        elapsed = time.time() - start_time
                        print(f"[{progress:3d}%] Writing gene files... ({done_count:,}/{total_files:,}) [{elapsed:.1f}s]")
        else:
            # Serial fallback
            for idx, args in enumerate(worker_args):
                _write_gene_file_worker(args)
                if (idx + 1) % max(1, total_files // 10) == 0 or (idx + 1) == total_files:
                    progress = 90 + int(((idx + 1) / total_files) * 8)
                    elapsed = time.time() - start_time
                    print(f"[{progress:3d}%] Writing gene files... ({idx + 1:,}/{total_files:,}) [{elapsed:.1f}s]")
    else:
        # Skip gene file writing
        output_path.mkdir(parents=True, exist_ok=True)
        print(f"[{90:3d}%] Skipping gene files (--manifest-only mode)...")

    # Step 7: Write manifest (98%)
    progress = 98
    elapsed = time.time() - start_time
    print(f"[{progress:3d}%] Writing manifest... [{elapsed:.1f}s]")

    # Get source filename without extension
    source_name = Path(input_file).stem

    # Build per-gene molecule counts
    gene_molecule_counts = {}
    for gene in unique_genes:
        assigned = len(gene_index.get(gene, []))
        entry = {"assigned": assigned}
        if has_unassigned:
            unassigned = len(unassigned_index.get(gene, []))
            entry["unassigned"] = unassigned
        gene_molecule_counts[gene] = entry

    # Create manifest in exact browser format
    manifest = {
        "version": "1.0",
        "created_at": datetime.datetime.utcnow().isoformat() + "Z",
        "dataset_id": f"temp_{int(time.time() * 1000)}",
        "name": source_name,
        "type": dataset_type,
        "statistics": {
            "total_molecules": total_molecules,
            "unique_genes": len(unique_genes),
            "spatial_dimensions": dimensions,
        },
        "genes": {
            "unique_gene_names": unique_genes,
            "molecule_counts": gene_molecule_counts,
        },
        "has_unassigned": has_unassigned,
        "processing": {
            "compression": "gzip",
            "coordinate_format": "float32_flat_array",
            "coordinate_range": "raw_rounded_2dp",
            "scaling_factor": 1.0,
            "created_by": "MERFISH Eyes - Single Molecule Viewer",
            "source_file": source_name,
        },
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
    if not manifest_only:
        print(f"    genes/")
        print(f"      {sanitize_gene_name(unique_genes[0])}.bin.gz")
        if has_unassigned:
            print(f"      {sanitize_gene_name(unique_genes[0])}{UNASSIGNED_SUFFIX}.bin.gz")
        print(f"      ...")
        print(f"      ({total_files} total gene files)")
    else:
        print(f"    (genes/ folder not created - use without --manifest-only to generate)")
    print()
    print("Statistics:")
    print(f"  Total molecules: {total_molecules:,}")
    print(f"  Unique genes: {len(unique_genes):,}")
    print(f"  Dimensions: {dimensions}D")
    print(f"  Coordinates: raw, rounded to 2dp")
    if has_unassigned:
        total_assigned = sum(len(v) for v in gene_index.values())
        total_unassigned = sum(len(v) for v in unassigned_index.values())
        print(f"  Assigned molecules: {total_assigned:,}")
        print(f"  Unassigned molecules: {total_unassigned:,}")
    print()
    print("Ready for S3 upload! Upload the entire output folder to:")
    print(f"  s3://your-bucket/datasets/{{datasetId}}/")


def process_directory(
    target_dir: str,
    output_folder: str,
    s3_prefix: str,
    link_column: str = "_sample_id",
    dataset_type: str = "merscope",
    gene_col: Optional[str] = None,
    x_col: Optional[str] = None,
    y_col: Optional[str] = None,
    z_col: Optional[str] = None,
    cell_id_col: Optional[str] = None,
    manifest_only: bool = False,
    num_workers: int = 1,
    sample_workers: int = 1,
):
    """
    Process all detected_transcripts.csv files found in a directory tree.
    Generates per-sample output folders and a mapping.json linking _sample_id to S3 URLs.

    num_workers: workers per sample (gene file writing parallelism)
    sample_workers: number of samples to process concurrently
    """
    total_start = time.time()

    target = Path(target_dir).resolve()
    if not target.is_dir():
        print(f"Error: '{target}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    # Discover sample directories
    print(f"Target directory: {target}")
    print(f"Sample workers: {sample_workers}")
    print(f"Write workers per sample: {num_workers}")
    print()
    sample_dirs = discover_transcript_dirs_bfs(target)

    if not sample_dirs:
        print("Error: No detected_transcripts.csv files found in any subdirectory.", file=sys.stderr)
        sys.exit(1)

    print(f"\nFound {len(sample_dirs)} sample(s) with transcript files:")
    for sample_id, dir_path in sample_dirs:
        print(f"  {sample_id} → {dir_path}")
    print()

    # Create output directory
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    # Strip trailing slash from s3_prefix
    s3_prefix = s3_prefix.rstrip("/")

    # Resolve transcript files for each sample before processing
    sample_tasks = []
    for sample_id, dir_path in sample_dirs:
        transcript_file = find_transcript_in_dir(str(dir_path))
        if transcript_file is None:
            print(f"  WARNING: Could not find transcript file in {dir_path}, skipping {sample_id}")
            continue
        sample_output = str(output_path / sample_id)
        sample_tasks.append((sample_id, transcript_file, sample_output))

    # Process samples
    mapping_links = {}
    failed = []

    if sample_workers > 1 and len(sample_tasks) > 1:
        # Parallel sample processing
        effective_sample_workers = min(sample_workers, len(sample_tasks))
        print(f"Processing {len(sample_tasks)} samples with {effective_sample_workers} sample workers...")

        worker_args = [
            (sample_id, str(transcript_file), sample_output,
             dataset_type, gene_col, x_col, y_col, z_col,
             cell_id_col, manifest_only, num_workers)
            for sample_id, transcript_file, sample_output in sample_tasks
        ]

        with ProcessPoolExecutor(max_workers=effective_sample_workers) as pool:
            futures = {pool.submit(_process_sample_worker, a): a[0] for a in worker_args}
            done_count = 0
            for future in as_completed(futures):
                sample_id, success, error = future.result()
                done_count += 1
                if success:
                    mapping_links[sample_id] = f"{s3_prefix}/{sample_id}"
                    print(f"[{done_count}/{len(sample_tasks)}] {sample_id}: OK")
                else:
                    failed.append(sample_id)
                    print(f"[{done_count}/{len(sample_tasks)}] {sample_id}: FAILED - {error}")
    else:
        # Serial sample processing
        for idx, (sample_id, transcript_file, sample_output) in enumerate(sample_tasks, 1):
            print(f"\n{'='*60}")
            print(f"[{idx}/{len(sample_tasks)}] Processing sample: {sample_id}")
            print(f"{'='*60}")
            print(f"  Transcript file: {transcript_file}")

            try:
                process_single_molecule_data(
                    input_file=str(transcript_file),
                    output_folder=sample_output,
                    dataset_type=dataset_type,
                    gene_col=gene_col,
                    x_col=x_col,
                    y_col=y_col,
                    z_col=z_col,
                    cell_id_col=cell_id_col,
                    manifest_only=manifest_only,
                    num_workers=num_workers,
                )
                mapping_links[sample_id] = f"{s3_prefix}/{sample_id}"
            except Exception as e:
                print(f"  ERROR processing {sample_id}: {e}", file=sys.stderr)
                failed.append(sample_id)

    # Write mapping.json
    mapping = {
        "linkColumn": link_column,
        "links": mapping_links,
    }

    mapping_file = output_path / "mapping.json"
    with open(mapping_file, "w") as f:
        json.dump(mapping, f, indent=2)

    # Summary
    elapsed = time.time() - total_start
    elapsed_str = f"{elapsed:.2f}s" if elapsed < 60 else f"{int(elapsed // 60)}m {elapsed % 60:.1f}s"

    print(f"\n{'='*60}")
    print(f"Directory processing complete! [{elapsed_str}]")
    print(f"{'='*60}")
    print(f"  Samples processed: {len(mapping_links)}/{len(sample_tasks)}")
    if failed:
        print(f"  Failed: {failed}")
    print(f"  Mapping file: {mapping_file}")
    print(f"\nMapping contents:")
    print(json.dumps(mapping, indent=2))


def main():
    parser = argparse.ArgumentParser(
        description="Process single molecule data to S3-compatible format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "input",
        help="Input parquet/CSV file OR directory containing sample sub-folders",
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
    parser.add_argument(
        "--cell-id-col",
        help="Column name for cell assignment (auto-detected as 'cell_id' for MERSCOPE). "
             "Molecules with value -1 are treated as unassigned and written to separate files. "
             "Use this flag to override the auto-detected column name",
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Only generate manifest.json.gz without creating gene files (faster)",
    )
    parser.add_argument(
        "--s3-prefix",
        help="S3 base URL prefix for mapping file (required for directory mode). "
             "Example: https://bucket.s3.region.amazonaws.com/prefix",
    )
    parser.add_argument(
        "--link-column",
        default="_sample_id",
        help="Cluster column name that the viewer reads for linking (default: _sample_id)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel workers for gene file writing within each sample (default: 1)",
    )
    parser.add_argument(
        "--sample-workers",
        type=int,
        default=1,
        help="Number of samples to process concurrently in directory mode (default: 1)",
    )

    args = parser.parse_args()

    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: Input path not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    if input_path.is_dir():
        # Directory mode
        if not args.s3_prefix:
            print("Error: --s3-prefix is required when input is a directory.", file=sys.stderr)
            sys.exit(1)

        process_directory(
            target_dir=args.input,
            output_folder=args.output_folder,
            s3_prefix=args.s3_prefix,
            link_column=args.link_column,
            dataset_type=args.dataset_type,
            gene_col=args.gene_col,
            x_col=args.x_col,
            y_col=args.y_col,
            z_col=args.z_col,
            cell_id_col=args.cell_id_col,
            manifest_only=args.manifest_only,
            num_workers=args.workers,
            sample_workers=args.sample_workers,
        )
    else:
        # Single file mode
        try:
            process_single_molecule_data(
                input_file=args.input,
                output_folder=args.output_folder,
                dataset_type=args.dataset_type,
                gene_col=args.gene_col,
                x_col=args.x_col,
                y_col=args.y_col,
                z_col=args.z_col,
                cell_id_col=args.cell_id_col,
                manifest_only=args.manifest_only,
                num_workers=args.workers,
            )
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
