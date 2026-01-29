#!/usr/bin/env python3
"""
Spatial Transcriptomics Data Processor
Converts H5AD, Xenium, and MERSCOPE datasets to chunked binary format for MERFISH Eyes

Usage:
    python process_spatial_data.py input.h5ad output_folder/ [--chunk-size GENES_PER_CHUNK]
    python process_spatial_data.py xenium_folder/ output_folder/ [--chunk-size GENES_PER_CHUNK]
    python process_spatial_data.py merscope_folder/ output_folder/ [--chunk-size GENES_PER_CHUNK]

Options:
    --chunk-size    Number of genes per chunk (default: auto, or 1 for single gene chunks)
"""

import argparse
import gzip
import json
import os
import shutil
import struct
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


# Color palette from lib/utils/color-palette.ts
DEFAULT_COLOR_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#00d9ff",
    "#e377c2", "#ffeb3b", "#bcbd22", "#17becf", "#ff006e", "#00ff00",
    "#ff1744", "#00e5ff", "#ffc107", "#e91e63", "#4caf50", "#ff9800",
    "#9c27b0", "#00bcd4", "#ffeb3b", "#f44336", "#3f51b5", "#8bc34a",
    "#ff5722", "#673ab7", "#03a9f4", "#cddc39", "#ff9100", "#7c4dff",
    "#00e676", "#ff3d00", "#651fff", "#1de9b6", "#ff6e40", "#d500f9",
    "#00b0ff", "#76ff03", "#ff1744", "#00e5ff",
]


def get_color_from_palette(index: int) -> str:
    """Get a color from the palette by index (cycles through)"""
    return DEFAULT_COLOR_PALETTE[index % len(DEFAULT_COLOR_PALETTE)]


def generate_color_palette(values: List[str]) -> Dict[str, str]:
    """Generate a color palette for unique values"""
    palette = {}
    for i, value in enumerate(values):
        palette[str(value)] = get_color_from_palette(i)
    return palette


def determine_chunk_size(num_genes: int, custom_chunk_size: Optional[int] = None) -> int:
    """Determine chunk size based on number of genes (matches TS logic)"""
    if custom_chunk_size is not None:
        return custom_chunk_size

    if num_genes < 100:
        return 50
    elif num_genes < 500:
        return 100
    elif num_genes < 2000:
        return 200
    elif num_genes < 10000:
        return 500
    else:
        return 1000


def is_categorical(values: np.ndarray, column_name: Optional[str] = None) -> bool:
    """
    Determine if data is categorical or numerical

    Algorithm:
    1. Always treat columns named "leiden" or "louvain" as categorical
    2. Check if values look like floats → numerical
    3. For integer-like values: if ≥80% are unique → numerical, otherwise categorical
    4. For non-numeric strings → categorical

    Args:
        values: Array of column values
        column_name: Optional column name (for special cases like "leiden")

    Returns:
        True if categorical, False if numerical
    """
    if not isinstance(values, np.ndarray):
        values = np.asarray(values)

    series = pd.Series(values)

    # Special case: leiden/louvain columns are always categorical
    if column_name:
        lower_name = column_name.lower()
        if "leiden" in lower_name or "louvain" in lower_name:
            return True

    # Filter out null/NaN values for analysis
    valid_series = series.dropna()
    if len(valid_series) == 0:
        return False

    # Check if values look like floats
    if pd.api.types.is_float_dtype(valid_series):
        # Check if any values have fractional parts
        float_count = 0
        sample_size = min(1000, len(valid_series))
        sample = valid_series.sample(n=sample_size, random_state=42) if len(valid_series) > sample_size else valid_series

        for val in sample:
            if isinstance(val, (float, np.floating)) and not float(val).is_integer():
                float_count += 1

        float_ratio = float_count / len(sample)

        # If majority (>50%) are floats with fractional parts → numerical
        if float_ratio > 0.5:
            return False

    # Check if values are numeric strings or numbers
    if series.dtype == "object":
        # Try to convert to numeric
        numeric_series = pd.to_numeric(valid_series, errors='coerce')
        numeric_ratio = numeric_series.notna().sum() / len(valid_series)

        # If not numbers at all (>50% non-numeric) → categorical (string labels)
        if numeric_ratio < 0.5:
            return True

        # Check for float-like strings (contain decimal point or scientific notation)
        float_count = 0
        sample_size = min(1000, len(valid_series))
        sample = valid_series.sample(n=sample_size, random_state=42) if len(valid_series) > sample_size else valid_series

        for val in sample:
            val_str = str(val).strip()
            if '.' in val_str or 'e' in val_str.lower():
                float_count += 1

        float_ratio = float_count / len(sample)
        if float_ratio > 0.5:
            return False

    # At this point, values are integer-like numbers or integer strings
    # Check uniqueness ratio: if ≥80% unique → numerical, otherwise categorical
    unique_count = valid_series.nunique()
    unique_ratio = unique_count / len(valid_series)

    return unique_ratio < 0.8


def normalize_coordinates(coords: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Normalize coordinates to [-1, 1] range
    Returns: (normalized_coords, scaling_factor)
    """
    # Flatten all coordinates to find global min/max
    all_coords = coords.flatten()
    min_val = np.min(all_coords)
    max_val = np.max(all_coords)
    range_val = max_val - min_val

    if range_val == 0:
        return coords, 1.0

    # Normalize to [-1, 1]
    normalized = ((coords - min_val) / range_val) * 2 - 1

    return normalized.astype(np.float32), float(range_val / 2)


def write_coordinate_binary(coords: np.ndarray, output_path: Path):
    """
    Write coordinates to binary format with gzip compression
    Format: [num_points: uint32, dimensions: uint32, [coordinates: float32]]
    """
    num_points, dimensions = coords.shape

    # Create binary data
    data = bytearray()

    # Header
    data.extend(struct.pack('<I', num_points))  # uint32 little-endian
    data.extend(struct.pack('<I', dimensions))  # uint32 little-endian

    # Coordinates (flatten row-major)
    for point in coords:
        for coord in point:
            data.extend(struct.pack('<f', coord))  # float32 little-endian

    # Compress and write
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, 'wb') as f:
        f.write(data)

    print(f"  ✓ Wrote {output_path.name}: {num_points} points, {dimensions}D")


def write_sparse_gene_data(indices: np.ndarray, values: np.ndarray) -> bytes:
    """
    Write sparse gene data to binary format
    Format: [num_cells: uint32, num_non_zero: uint32, [indices: uint32[]], [values: float32[]]]
    """
    num_cells = len(indices)
    num_non_zero = len(indices)

    data = bytearray()

    # Sparse data header
    data.extend(struct.pack('<I', num_cells))  # uint32
    data.extend(struct.pack('<I', num_non_zero))  # uint32

    # Indices (cell indices with non-zero values)
    for idx in indices:
        data.extend(struct.pack('<I', int(idx)))  # uint32

    # Values (expression values)
    for val in values:
        data.extend(struct.pack('<f', float(val)))  # float32

    return bytes(data)


def write_expression_chunk(
    chunk_id: int,
    gene_indices: List[int],
    gene_data_list: List[Tuple[np.ndarray, np.ndarray]],
    total_cells: int,
    output_path: Path
):
    """
    Write expression chunk to binary format with gzip compression
    """
    num_genes = len(gene_indices)

    # Build gene data sections first to know offsets
    gene_data_sections = []
    for indices, values in gene_data_list:
        gene_data = write_sparse_gene_data(indices, values)
        gene_data_sections.append(gene_data)

    # Build the binary chunk
    data = bytearray()

    # Header (16 bytes)
    data.extend(struct.pack('<I', 1))  # version
    data.extend(struct.pack('<I', num_genes))
    data.extend(struct.pack('<I', chunk_id))
    data.extend(struct.pack('<I', total_cells))

    # Calculate gene table offset (after header)
    gene_table_start = 16
    gene_table_size = num_genes * 24  # 24 bytes per gene
    first_gene_data_offset = gene_table_start + gene_table_size

    # Write gene table (24 bytes per gene)
    current_offset = first_gene_data_offset
    for i, (gene_idx, gene_data) in enumerate(zip(gene_indices, gene_data_sections)):
        indices, values = gene_data_list[i]
        num_non_zero = len(indices)
        data_size = len(gene_data)

        data.extend(struct.pack('<I', gene_idx))  # gene_index: uint32
        data.extend(struct.pack('<I', current_offset))  # data_offset: uint32
        data.extend(struct.pack('<I', data_size))  # data_size: uint32
        data.extend(struct.pack('<I', data_size))  # uncompressed_size: uint32
        data.extend(struct.pack('<I', num_non_zero))  # num_non_zero: uint32
        data.extend(struct.pack('<I', 0))  # reserved: uint32

        current_offset += data_size

    # Write gene data sections
    for gene_data in gene_data_sections:
        data.extend(gene_data)

    # Compress and write
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, 'wb') as f:
        f.write(data)

    print(f"  ✓ Wrote chunk_{chunk_id:05d}.bin.gz: {num_genes} genes, {len(data)} bytes (uncompressed)")


def detect_input_format(input_path: Path) -> str:
    """Detect if input is H5AD, Xenium, or MERSCOPE"""
    # Debug output
    print(f"🔍 Detecting format for: {input_path}")
    print(f"   - Exists: {input_path.exists()}")
    print(f"   - Is file: {input_path.is_file()}")
    print(f"   - Is dir: {input_path.is_dir()}")

    if input_path.is_file():
        if input_path.suffix == '.h5ad':
            return 'h5ad'
        else:
            raise ValueError(f"Unknown file format: {input_path}")

    elif input_path.is_dir():
        # List files in directory for debugging
        files_in_dir = [f.name for f in input_path.iterdir()]
        print(f"   - Files in directory: {files_in_dir}")

        # Check for Xenium files
        xenium_files = ['cells.csv', 'cells.csv.gz', 'cell_feature_matrix']
        has_xenium = any((input_path / f).exists() for f in xenium_files)
        print(f"   - Has Xenium files: {has_xenium}")

        # Check for MERSCOPE files
        merscope_files = ['cell_metadata.csv', 'cell_by_gene.csv']
        has_merscope = any((input_path / f).exists() for f in merscope_files)
        print(f"   - Has MERSCOPE files: {has_merscope}")

        # Debug: Check each file individually
        for f in merscope_files:
            file_path = input_path / f
            print(f"   - Checking {f}: exists={file_path.exists()}, is_file={file_path.is_file()}")

        if has_xenium:
            return 'xenium'
        elif has_merscope:
            return 'merscope'
        else:
            raise ValueError("Could not detect format. Expected H5AD file, Xenium folder, or MERSCOPE folder.")

    else:
        raise ValueError(f"Input path does not exist: {input_path}")


def load_h5ad_data(input_path: Path):
    """Load H5AD file using anndata and extract spatial coordinates"""
    try:
        import anndata as ad
    except ImportError:
        raise ImportError("anndata is required for H5AD files. Install with: pip install anndata")

    print(f"📂 Loading H5AD file: {input_path.name}")
    adata = ad.read_h5ad(input_path)
    print(f"  ✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Extract spatial coordinates
    spatial_coords = None

    # Try obsm first
    for key in ['X_spatial', 'spatial']:
        if key in adata.obsm:
            spatial_coords = adata.obsm[key]
            print(f"  ✓ Found spatial coordinates in obsm['{key}']")
            break

    # Fallback to obs columns (matching TypeScript H5adAdapter logic)
    if spatial_coords is None:
        print("  ⚠ No spatial coordinates in obsm, checking obs columns...")
        x_candidates = [
            'center_x', 'centerX', 'centerx', 'x', 'X',
            'x_centroid', 'centroid_x', 'x_location',
            'center_X', 'spatial_x', 'spatial_X'
        ]
        y_candidates = [
            'center_y', 'centerY', 'centery', 'y', 'Y',
            'y_centroid', 'centroid_y', 'y_location',
            'center_Y', 'spatial_y', 'spatial_Y'
        ]
        z_candidates = [
            'center_z', 'centerZ', 'centerz', 'z', 'Z',
            'z_centroid', 'centroid_z', 'z_location',
            'center_Z', 'spatial_z', 'spatial_Z'
        ]

        x_col = next((col for col in x_candidates if col in adata.obs.columns), None)
        y_col = next((col for col in y_candidates if col in adata.obs.columns), None)

        if x_col and y_col:
            print(f"  ✓ Found coordinates in obs: {x_col}, {y_col}")
            z_col = next((col for col in z_candidates if col in adata.obs.columns), None)
            if z_col:
                print(f"  ✓ Found z coordinate: {z_col}")
                spatial_coords = np.column_stack([
                    adata.obs[x_col].values,
                    adata.obs[y_col].values,
                    adata.obs[z_col].values
                ])
            else:
                spatial_coords = np.column_stack([
                    adata.obs[x_col].values,
                    adata.obs[y_col].values
                ])

    if spatial_coords is None:
        raise ValueError(
            "No spatial coordinates found! Checked: obsm['X_spatial'], obsm['spatial'], "
            "and obs columns (center_x/y/z, x/y/z, etc.)"
        )

    # Extract expression matrix
    expr_matrix = adata.X
    gene_names = list(adata.var_names)

    # Extract observation columns (clusters, metadata)
    obs_columns = {}
    # Exclude coordinate columns from obs
    coord_cols = ['center_x', 'centerX', 'center_y', 'centerY', 'center_z', 'centerZ',
                  'x', 'y', 'z', 'X', 'Y', 'Z',
                  'x_centroid', 'y_centroid', 'z_centroid',
                  'centroid_x', 'centroid_y', 'centroid_z']
    for col in adata.obs.columns:
        if col not in coord_cols:
            obs_columns[col] = adata.obs[col].values

    return spatial_coords, expr_matrix, gene_names, obs_columns, adata


def load_xenium_data(input_path: Path):
    """Load Xenium folder structure"""
    print(f"📂 Loading Xenium folder: {input_path.name}")

    # Load cells.csv or cells.csv.gz
    cells_file = None
    for filename in ['cells.csv', 'cells.csv.gz']:
        candidate = input_path / filename
        if candidate.exists():
            cells_file = candidate
            break

    if cells_file is None:
        raise FileNotFoundError("cells.csv or cells.csv.gz not found in Xenium folder")

    print(f"  Loading {cells_file.name}...")
    if cells_file.suffix == '.gz':
        cells_df = pd.read_csv(cells_file, compression='gzip')
    else:
        cells_df = pd.read_csv(cells_file)

    print(f"  ✓ Loaded {len(cells_df)} cells")

    # Try to load expression matrix (cell_feature_matrix or cell_by_gene)
    expr_matrix = None
    gene_names = None
    features_file = None

    # Option 1: Try wide format in cells.csv (genes as columns)
    # Skip coordinate and metadata columns
    coordinate_cols = ['cell_id', 'x_centroid', 'y_centroid', 'z_centroid',
                       'nucleus_x', 'nucleus_y', 'cell_centroid_x', 'cell_centroid_y']
    gene_cols = [col for col in cells_df.columns if col not in coordinate_cols]

    # Check if we have numeric gene columns
    numeric_cols = []
    for col in gene_cols:
        if pd.api.types.is_numeric_dtype(cells_df[col]):
            numeric_cols.append(col)

    if len(numeric_cols) > 10:  # Likely wide format with genes
        print(f"  Detected wide format: {len(numeric_cols)} gene columns")
        expr_matrix = cells_df[numeric_cols].values
        gene_names = numeric_cols

    # Option 2: Try loading features.tsv
    if gene_names is None:
        for path in [input_path / 'features.tsv', input_path / 'features.tsv.gz',
                     input_path / 'cell_feature_matrix' / 'features.tsv',
                     input_path / 'cell_feature_matrix' / 'features.tsv.gz']:
            if path.exists():
                features_file = path
                break

        if features_file:
            print(f"  Loading {features_file.name}...")
            if features_file.suffix == '.gz':
                features_df = pd.read_csv(features_file, sep='\t', header=None, compression='gzip')
            else:
                features_df = pd.read_csv(features_file, sep='\t', header=None)
            # Per Xenium spec, column 1 (0-indexed) holds the gene symbols; fallback to column 0 if missing
            gene_names = features_df[1].tolist() if features_df.shape[1] > 1 else features_df[0].tolist()
            print(f"  ✓ Loaded {len(gene_names)} genes from features file")

    # Option 3: Load matrix.mtx or cell_feature_matrix.h5
    if expr_matrix is None:
        matrix_file = None
        temp_extract_dir: Optional[Path] = None

        for path in [
            input_path / "matrix.mtx",
            input_path / "matrix.mtx.gz",
            input_path / "cell_feature_matrix" / "matrix.mtx",
            input_path / "cell_feature_matrix" / "matrix.mtx.gz",
        ]:
            if path.exists():
                matrix_file = path
                break

        if matrix_file is None:
            for archive in [
                input_path / "cell_feature_matrix.zip",
                input_path / "cell_feature_matrix.tar",
                input_path / "cell_feature_matrix.tar.gz",
            ]:
                if archive.exists():
                    print(f"  Extracting {archive.name}...")
                    import tempfile
                    import tarfile
                    import zipfile

                    temp_extract_dir = Path(tempfile.mkdtemp(prefix="xenium_cfm_"))
                    try:
                        if archive.suffix == ".zip":
                            with zipfile.ZipFile(archive, "r") as zf:
                                zf.extractall(temp_extract_dir)
                        else:
                            mode = "r:gz" if archive.name.endswith(".tar.gz") else "r"
                            with tarfile.open(archive, mode) as tf:
                                tf.extractall(temp_extract_dir)
                        for candidate in temp_extract_dir.rglob("matrix.mtx*"):
                            matrix_file = candidate
                            break
                        if not features_file:
                            for candidate in temp_extract_dir.rglob("features.tsv*"):
                                if candidate.exists():
                                    features_file = candidate
                                    break
                        if matrix_file:
                            print(f"  ✓ Found matrix in archive: {matrix_file}")
                        else:
                            print("  ⚠ matrix.mtx not found in archive")
                    except Exception as e:
                        print(f"  ⚠ Failed to extract {archive.name}: {e}")
                        matrix_file = None
                        if temp_extract_dir:
                            shutil.rmtree(temp_extract_dir, ignore_errors=True)
                            temp_extract_dir = None
                    if matrix_file:
                        break

        if features_file and gene_names is None and features_file.exists():
            print(f"  Loading {features_file.name}...")
            if features_file.suffix.endswith('.gz'):
                features_df = pd.read_csv(features_file, sep='\t', header=None, compression='gzip')
            else:
                features_df = pd.read_csv(features_file, sep='\t', header=None)
            # Per Xenium spec, column 1 (0-indexed) holds the gene symbols; fallback to column 0 if missing
            gene_names = features_df[1].tolist() if features_df.shape[1] > 1 else features_df[0].tolist()
            print(f"  ✓ Loaded {len(gene_names)} genes from features file")

        if matrix_file:
            print(f"  Loading expression matrix from {matrix_file.name}...")
            matrix = None
            if matrix_file.suffix == '.gz' or matrix_file.name.endswith('.mtx.gz'):
                with gzip.open(matrix_file, 'rb') as f:
                    matrix = mmread(f)
            elif matrix_file.suffix == '.mtx':
                matrix = mmread(str(matrix_file))
            elif matrix_file.suffix == '.h5':
                try:
                    import anndata as ad

                    adata_matrix = ad.read_10x_h5(str(matrix_file))
                    matrix = adata_matrix.X
                    if gene_names is None:
                        gene_names = adata_matrix.var_names.tolist()
                except Exception as e:
                    print(f"  ⚠ Failed to read {matrix_file.name}: {e}")

            if matrix is not None:
                expr_matrix = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
                print(f"  ✓ Loaded sparse matrix: {expr_matrix.shape[0]} x {expr_matrix.shape[1]}")

                # Determine orientation (genes vs cells)
                if expr_matrix.shape[0] == len(cells_df):
                    pass  # Already cells x genes
                elif expr_matrix.shape[1] == len(cells_df):
                    expr_matrix = expr_matrix.transpose().tocsr()
                    print("  ↺ Transposed matrix to align with cells")
                else:
                    print("  ⚠ Matrix dimensions do not match cell count; continuing with placeholder")
                    expr_matrix = None

                # Ensure gene names align with matrix columns
                if expr_matrix is not None:
                    if gene_names is None:
                        if features_file and features_file.exists():
                            pass  # already loaded above if file present
                        else:
                            gene_names = [f"Gene{i+1}" for i in range(expr_matrix.shape[1])]
                            print(f"  ⚠ Features file missing; generated {len(gene_names)} placeholder gene names")
                    elif len(gene_names) != expr_matrix.shape[1]:
                        print(
                            f"  ⚠ Gene name count ({len(gene_names)}) does not match matrix columns ({expr_matrix.shape[1]}). "
                            "Truncating to match."
                        )
                        gene_names = gene_names[:expr_matrix.shape[1]]
        if temp_extract_dir:
            shutil.rmtree(temp_extract_dir, ignore_errors=True)

    if expr_matrix is None:
        print("  ⚠ No expression matrix found, creating placeholder")
        gene_names = ['Gene1']  # Placeholder
        expr_matrix = np.zeros((len(cells_df), 1))

    return cells_df, expr_matrix, gene_names


def load_merscope_data(input_path: Path):
    """Load MERSCOPE folder structure"""
    print(f"📂 Loading MERSCOPE folder: {input_path.name}")

    # Load cell_metadata.csv
    metadata_file = input_path / 'cell_metadata.csv'
    if not metadata_file.exists():
        raise FileNotFoundError("cell_metadata.csv not found in MERSCOPE folder")

    print(f"  Loading {metadata_file.name}...")
    metadata_df = pd.read_csv(metadata_file)
    print(f"  ✓ Loaded {len(metadata_df)} cells")

    # Load cell_by_gene.csv (expression matrix)
    expr_file = input_path / 'cell_by_gene.csv'
    expr_matrix = None
    gene_names = None

    if expr_file.exists():
        print(f"  Loading {expr_file.name}...")
        expr_df = pd.read_csv(expr_file, index_col=0)
        gene_names = expr_df.columns.tolist()
        expr_matrix = expr_df.values
        print(f"  ✓ Loaded expression matrix: {len(gene_names)} genes")
    else:
        print("  ⚠ cell_by_gene.csv not found, creating placeholder")
        gene_names = ['Gene1']
        expr_matrix = np.zeros((len(metadata_df), 1))

    return metadata_df, expr_matrix, gene_names


def process_dataset(
    input_path: Path,
    output_dir: Path,
    custom_chunk_size: Optional[int] = None,
    data_format: Optional[str] = None
):
    """Main processing function that handles all formats"""

    print(f"\n{'='*60}")
    print(f"Processing Spatial Transcriptomics Data")
    print(f"Input: {input_path}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}\n")

    # Auto-detect format if not specified
    if data_format is None:
        data_format = detect_input_format(input_path)

    print(f"📋 Detected format: {data_format.upper()}\n")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data based on format
    if data_format == 'h5ad':
        spatial_coords, expr_matrix, gene_names, obs_columns, adata = load_h5ad_data(input_path)
        dataset_name = input_path.stem
        dataset_type = 'h5ad'

        # Process embeddings (UMAP, etc.)
        embeddings = {}
        MAX_EMBEDDING_DIMS = 3  # Maximum dimensions to save for embeddings
        for key in adata.obsm.keys():
            if key not in ['X_spatial', 'spatial'] and key.startswith('X_'):
                embedding_name = key[2:].lower()  # Remove 'X_' prefix
                embedding_coords = adata.obsm[key]

                # Limit to first 3 dimensions
                if embedding_coords.shape[1] > MAX_EMBEDDING_DIMS:
                    embedding_coords = embedding_coords[:, :MAX_EMBEDDING_DIMS]
                    print(f"  ℹ {embedding_name}: Limiting to first {MAX_EMBEDDING_DIMS} dimensions (from {adata.obsm[key].shape[1]})")

                embeddings[embedding_name] = embedding_coords

    elif data_format == 'xenium':
        cells_df, expr_matrix, gene_names = load_xenium_data(input_path)
        dataset_name = input_path.name
        dataset_type = 'xenium'
        embeddings = {}  # No embeddings for Xenium

        # Extract spatial coordinates
        coord_candidates_x = ['x_centroid', 'cell_centroid_x', 'center_x']
        coord_candidates_y = ['y_centroid', 'cell_centroid_y', 'center_y']
        coord_candidates_z = ['z_centroid', 'cell_centroid_z', 'center_z']

        cells_columns_lower = {col.lower(): col for col in cells_df.columns}

        def find_column(columns_lower: Dict[str, str], candidates: List[str]) -> Optional[str]:
            for cand in candidates:
                actual = columns_lower.get(cand.lower())
                if actual:
                    return actual
            return None

        x_col = find_column(cells_columns_lower, coord_candidates_x)
        y_col = find_column(cells_columns_lower, coord_candidates_y)
        z_col = find_column(cells_columns_lower, coord_candidates_z)

        if not x_col or not y_col:
            raise ValueError("Could not find x/y coordinate columns in cells.csv")

        if z_col and z_col in cells_df.columns:
            spatial_coords = np.column_stack([
                cells_df[x_col].values,
                cells_df[y_col].values,
                cells_df[z_col].values
            ])
        else:
            spatial_coords = np.column_stack([
                cells_df[x_col].values,
                cells_df[y_col].values
            ])

        # Extract observation columns (clusters, metadata)
        obs_columns = {}
        exclude_cols = {col for col in [x_col, y_col, z_col] if col}
        for cand in ['cell_id', 'id', 'barcode', 'barcodes', 'cell']:
            actual = find_column(cells_columns_lower, [cand])
            if actual:
                exclude_cols.add(actual)
        for col in cells_df.columns:
            if col not in exclude_cols and col not in gene_names:
                obs_columns[col] = cells_df[col].values

        analysis_dir = input_path / 'analysis' / 'clustering'
        analysis_extract_dir: Optional[Path] = None

        if not analysis_dir.exists():
            print("  ℹ analysis/clustering not found on disk; searching archives...")
            for archive in [
                input_path / 'analysis.tar.gz',
                input_path / 'analysis.tgz',
                input_path / 'analysis.tar',
                input_path / 'analysis.zip',
            ]:
                if not archive.exists():
                    continue
                print(f"  Extracting {archive.name} to inspect clustering results...")
                import tempfile
                import tarfile
                import zipfile

                analysis_extract_dir = Path(tempfile.mkdtemp(prefix="xenium_analysis_"))
                try:
                    if archive.suffix == '.zip':
                        with zipfile.ZipFile(archive, "r") as zf:
                            zf.extractall(analysis_extract_dir)
                    else:
                        mode = "r:gz" if archive.name.endswith((".tar.gz", ".tgz")) else "r"
                        with tarfile.open(archive, mode) as tf:
                            tf.extractall(analysis_extract_dir)

                    candidate = next(
                        (
                            p
                            for p in analysis_extract_dir.rglob("clustering")
                            if p.is_dir() and p.parent.name == "analysis"
                        ),
                        None,
                    )
                    if candidate:
                        analysis_dir = candidate
                        print(f"  ✓ Found clustering directory inside archive: {analysis_dir}")
                        break
                    else:
                        print("  ⚠ Extracted archive but could not find analysis/clustering directory")
                        shutil.rmtree(analysis_extract_dir, ignore_errors=True)
                        analysis_extract_dir = None
                except Exception as e:
                    print(f"  ⚠ Failed to extract {archive.name}: {e}")
                    shutil.rmtree(analysis_extract_dir, ignore_errors=True)
                    analysis_extract_dir = None

        if analysis_dir.exists():
            print("  Looking for clustering results in analysis/clustering ...")
            cluster_files = sorted(
                p
                for p in analysis_dir.rglob('*')
                if p.is_file()
                and not p.name.startswith("._")
                and any(
                    p.name.lower().endswith(ext)
                    for ext in ('.csv', '.csv.gz', '.tsv', '.tsv.gz')
                )
            )
            preferred_subdir = "analysis/clustering/gene_expression_graphclust"
            cluster_files.sort(
                key=lambda p: (
                    0
                    if preferred_subdir in str(p.as_posix()).lower()
                    else 1,
                    str(p),
                )
            )
            print(f"  ℹ Found {len(cluster_files)} candidate cluster files")
            id_candidates = ['cell_id', 'cell', 'id', 'barcode', 'barcodes', 'cellid', 'cellId']
            label_candidates = [
                'cluster',
                'clusters',
                'cell_type',
                'celltype',
                'annotation',
                'annotations',
                'label',
                'labels',
                'group',
                'subclass',
                'class',
            ]
            cell_id_key = (
                find_column(cells_columns_lower, ['cell_id', 'id', 'barcode', 'barcodes', 'cell'])
                or cells_df.columns[0]
            )
            for cluster_file in cluster_files:
                print(f"    ℹ Trying {cluster_file}")
                try:
                    suffixes = [s.lower() for s in cluster_file.suffixes]
                    is_tsv = '.tsv' in suffixes
                    cluster_df = pd.read_csv(
                        cluster_file,
                        sep='\t' if is_tsv else ',',
                        compression='infer'
                    )
                except Exception as e:
                    print(f"    ⚠ Failed to read {cluster_file.name}: {e}")
                    continue
                cluster_columns_lower = {col.lower(): col for col in cluster_df.columns}
                id_col = find_column(cluster_columns_lower, id_candidates)
                label_col = find_column(cluster_columns_lower, label_candidates)
                if not id_col or not label_col:
                    print(f"    ⚠ Missing id/label columns in {cluster_file.name}")
                    continue
                joined = cells_df.merge(
                    cluster_df[[id_col, label_col]],
                    left_on=cell_id_key,
                    right_on=id_col,
                    how='left'
                )
                non_null = joined[label_col].notna().sum()
                if non_null == 0:
                    print(f"    ⚠ No matching rows when joining {cluster_file.name}")
                    continue
                obs_columns[label_col] = joined[label_col].fillna('').values
                print(f"    ✓ Imported clusters from {cluster_file.name} using column '{label_col}' ({non_null} matches)")
                break
            else:
                print("  ⚠ Exhausted cluster files without importing labels")
        if analysis_extract_dir:
            shutil.rmtree(analysis_extract_dir, ignore_errors=True)

    elif data_format == 'merscope':
        metadata_df, expr_matrix, gene_names = load_merscope_data(input_path)
        dataset_name = input_path.name
        dataset_type = 'merscope'
        embeddings = {}  # No embeddings for MERSCOPE

        # Extract spatial coordinates
        coord_candidates_x = ['center_x', 'centroid_x', 'x']
        coord_candidates_y = ['center_y', 'centroid_y', 'y']
        coord_candidates_z = ['center_z', 'centroid_z', 'z']

        x_col = next((c for c in coord_candidates_x if c in metadata_df.columns), None)
        y_col = next((c for c in coord_candidates_y if c in metadata_df.columns), None)
        z_col = next((c for c in coord_candidates_z if c in metadata_df.columns), None)

        if not x_col or not y_col:
            raise ValueError("Could not find x/y coordinate columns in cell_metadata.csv")

        if z_col and z_col in metadata_df.columns:
            spatial_coords = np.column_stack([
                metadata_df[x_col].values,
                metadata_df[y_col].values,
                metadata_df[z_col].values
            ])
        else:
            spatial_coords = np.column_stack([
                metadata_df[x_col].values,
                metadata_df[y_col].values
            ])

        # Extract observation columns
        obs_columns = {}
        exclude_cols = [x_col, y_col, z_col, 'cell_id', 'id', 'EntityID']
        for col in metadata_df.columns:
            if col not in exclude_cols:
                obs_columns[col] = metadata_df[col].values

    else:
        raise ValueError(f"Unknown format: {data_format}")

    # Continue with common processing
    num_cells = len(spatial_coords)
    num_genes = len(gene_names)
    spatial_dims = spatial_coords.shape[1]

    print(f"\n📊 Dataset Summary:")
    print(f"  Cells: {num_cells:,}")
    print(f"  Genes: {num_genes:,}")
    print(f"  Spatial dimensions: {spatial_dims}D\n")

    # Process using the same logic as H5AD
    # 1. Normalize and write spatial coordinates
    print("📍 Processing spatial coordinates...")
    normalized_spatial, scaling_factor = normalize_coordinates(spatial_coords)
    coords_dir = output_dir / 'coords'
    write_coordinate_binary(normalized_spatial, coords_dir / 'spatial.bin.gz')

    # 2. Process embeddings (if any)
    available_embeddings = []
    if embeddings:
        print("\n🗺️  Processing embeddings...")
        for emb_name, emb_coords in embeddings.items():
            normalized_emb, _ = normalize_coordinates(emb_coords)
            write_coordinate_binary(normalized_emb, coords_dir / f'{emb_name}.bin.gz')
            available_embeddings.append(emb_name)
            print(f"  ✓ {emb_name}: {emb_coords.shape[1]}D")
    else:
        print("\n  ℹ No embeddings to process")

    # 3. Process observation columns (clusters)
    print("\n🏷️  Processing observation columns...")
    obs_dir = output_dir / 'obs'
    palettes_dir = output_dir / 'palettes'
    obs_metadata = {}
    cluster_count = 0

    for col_name, col_values in obs_columns.items():
        categorical = is_categorical(col_values, col_name)
        series = pd.Series(col_values)
        series_for_unique = (
            series.astype("string") if series.dtype == "object" else series
        )
        unique_count = int(series_for_unique.nunique(dropna=False))

        obs_metadata[col_name] = {
            'type': 'categorical' if categorical else 'numerical',
            'unique_values': unique_count
        }

        # Write column values
        if categorical:
            cat_series = series.astype("string")
            values_list = cat_series.fillna("").tolist()
        else:
            numeric_series = pd.to_numeric(series, errors="coerce")
            values_list = [
                None if pd.isna(v) else float(v) for v in numeric_series.tolist()
            ]
        obs_file = obs_dir / f'{col_name}.json.gz'
        obs_file.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(obs_file, 'wt', encoding='utf-8') as f:
            json.dump(values_list, f)

        print(f"  ✓ {col_name}: {obs_metadata[col_name]['type']} ({obs_metadata[col_name]['unique_values']} unique)")

        # Generate palette for categorical
        if categorical:
            cluster_count += 1
            unique_values = [
                str(v)
                for v in cat_series.dropna().unique().tolist()
            ]
            palette = generate_color_palette(unique_values)

            palette_file = palettes_dir / f'{col_name}.json'
            palette_file.parent.mkdir(parents=True, exist_ok=True)
            with open(palette_file, 'w') as f:
                json.dump(palette, f, indent=2)

            print(f"    → Palette generated with {len(unique_values)} colors")

    # Write obs metadata
    with open(obs_dir / 'metadata.json', 'w') as f:
        json.dump(obs_metadata, f, indent=2)

    # 3. Process expression matrix
    print("\n🧬 Processing expression matrix...")
    chunk_size = determine_chunk_size(num_genes, custom_chunk_size)
    num_chunks = (num_genes + chunk_size - 1) // chunk_size

    print(f"  Chunk size: {chunk_size} genes/chunk")
    print(f"  Total chunks: {num_chunks}")

    # Build expression index
    expr_index = {
        'total_genes': num_genes,
        'num_chunks': num_chunks,
        'chunk_size': chunk_size,
        'genes': []
    }

    expr_dir = output_dir / 'expr'

    # Process chunks
    for chunk_id in range(num_chunks):
        start_gene = chunk_id * chunk_size
        end_gene = min(start_gene + chunk_size, num_genes)

        gene_indices = []
        gene_data_list = []

        for gene_idx in range(start_gene, end_gene):
            gene_name = gene_names[gene_idx]

            # Extract gene expression column
            if sparse.issparse(expr_matrix):
                gene_col = expr_matrix[:, gene_idx].toarray().flatten()
            else:
                gene_col = expr_matrix[:, gene_idx]

            # Get non-zero indices and values
            non_zero_mask = gene_col != 0
            non_zero_indices = np.where(non_zero_mask)[0]
            non_zero_values = gene_col[non_zero_mask]

            gene_indices.append(gene_idx)
            gene_data_list.append((non_zero_indices, non_zero_values))

            # Add to expression index
            expr_index['genes'].append({
                'name': gene_name,
                'chunk_id': chunk_id,
                'position_in_chunk': gene_idx - start_gene
            })

        # Write chunk
        chunk_file = expr_dir / f'chunk_{chunk_id:05d}.bin.gz'
        write_expression_chunk(chunk_id, gene_indices, gene_data_list, num_cells, chunk_file)

    # Write expression index
    with open(expr_dir / 'index.json', 'w') as f:
        json.dump(expr_index, f, indent=2)
    print(f"\n  ✓ Expression index written")

    # 4. Create manifest
    print("\n📋 Creating manifest...")
    manifest = {
        'version': '1.0',
        'dataset_id': 'local_dataset',
        'name': dataset_name,
        'type': dataset_type,
        'statistics': {
            'total_cells': int(num_cells),
            'total_genes': int(num_genes),
            'spatial_dimensions': int(spatial_dims),
            'available_embeddings': available_embeddings,
            'cluster_count': cluster_count
        },
        'processing': {
            'spatial_scaling_factor': float(scaling_factor),
            'chunk_size': chunk_size,
            'num_chunks': num_chunks,
            'created_by': 'process_spatial_data.py'
        },
        'files': {
            'coordinates': ['spatial'],
            'expression_chunks': num_chunks,
            'observation_columns': list(obs_metadata.keys())
        }
    }

    with open(output_dir / 'manifest.json', 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"  ✓ Manifest written")

    # Summary
    print(f"\n{'='*60}")
    print(f"✅ Processing complete!")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Convert spatial transcriptomics data to chunked binary format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # H5AD file
  python process_spatial_data.py data.h5ad output/

  # Xenium folder
  python process_spatial_data.py xenium_output/ output/

  # MERSCOPE folder
  python process_spatial_data.py merscope_output/ output/

  # Custom chunk size
  python process_spatial_data.py data.h5ad output/ --chunk-size 100
        """
    )

    parser.add_argument('input', type=Path, help='Input H5AD file or Xenium/MERSCOPE folder')
    parser.add_argument('output', type=Path, help='Output directory for chunked data')
    parser.add_argument('--chunk-size', type=int, default=None,
                        help='Genes per chunk (default: auto-determined)')
    parser.add_argument('--format', type=str, choices=['h5ad', 'xenium', 'merscope'],
                        help='Force input format (auto-detected if not specified)')

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"❌ Error: Input path not found: {args.input}")
        sys.exit(1)

    # Process
    try:
        process_dataset(args.input, args.output, args.chunk_size, args.format)
    except Exception as e:
        print(f"\n❌ Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
