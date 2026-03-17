
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
import gc
import gzip
import json
import os
import shutil
import struct
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


# ─────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────
_t_start = None

def fmt_elapsed(seconds):
    """Format elapsed time as human-readable string."""
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(seconds, 60)
    return f"{int(m)}m {s:.1f}s"

def log(msg, t0=None):
    """Print a timestamped log message. If t0 given, also prints elapsed."""
    elapsed = ""
    if t0 is not None:
        elapsed = f" [{fmt_elapsed(time.perf_counter() - t0)}]"
    print(f"[{time.strftime('%H:%M:%S')}]{elapsed} {msg}", flush=True)


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
    (Synced with lib/utils/column-type-detection.ts)

    Algorithm:
    1. Always treat columns named "leiden" or "louvain" as categorical
    2. Known numerical columns (n_genes, total_counts, etc.) → numerical
    3. Check if values look like floats → numerical
    4. For non-numeric strings → categorical
    5. For integer-like values: if ≥80% are unique → numerical, otherwise categorical

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

        # Known numerical columns (counts, QC metrics) — synced with column-type-detection.ts
        numerical_patterns = [
            "n_genes",
            "total_counts",
            "n_counts",
            "pct_counts",
            "log1p_",
            "n_cells",
            "doublet_score",
        ]
        if any(p in lower_name for p in numerical_patterns):
            return False

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

    # Header
    header = struct.pack('<II', num_points, dimensions)

    # Coordinates as contiguous float32 bytes (row-major)
    coord_bytes = np.ascontiguousarray(coords, dtype=np.float32).tobytes()

    # Compress and write
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, 'wb') as f:
        f.write(header)
        f.write(coord_bytes)

    print(f"  ✓ Wrote {output_path.name}: {num_points} points, {dimensions}D")


def write_sparse_gene_data(indices: np.ndarray, values: np.ndarray) -> bytes:
    """
    Write sparse gene data to binary format
    Format: [num_cells: uint32, num_non_zero: uint32, [indices: uint32[]], [values: float32[]]]
    """
    num_non_zero = len(indices)

    # Header
    header = struct.pack('<II', num_non_zero, num_non_zero)

    # Indices and values as contiguous typed arrays
    idx_bytes = np.asarray(indices, dtype=np.uint32).tobytes()
    val_bytes = np.asarray(values, dtype=np.float32).tobytes()

    return header + idx_bytes + val_bytes


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

    # Header (16 bytes)
    header = struct.pack('<IIII', 1, num_genes, chunk_id, total_cells)

    # Calculate gene table offset (after header)
    gene_table_start = 16
    gene_table_size = num_genes * 24  # 24 bytes per gene
    first_gene_data_offset = gene_table_start + gene_table_size

    # Build gene table (24 bytes per gene) as a single array
    gene_table = bytearray()
    current_offset = first_gene_data_offset
    for i, (gene_idx, gene_data) in enumerate(zip(gene_indices, gene_data_sections)):
        num_non_zero = len(gene_data_list[i][0])
        data_size = len(gene_data)

        gene_table.extend(struct.pack('<IIIIII',
            gene_idx, current_offset, data_size, data_size, num_non_zero, 0))

        current_offset += data_size

    # Compress and write
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, 'wb') as f:
        f.write(header)
        f.write(gene_table)
        for gene_data in gene_data_sections:
            f.write(gene_data)


def _write_chunk_worker(args):
    """Worker function for parallel chunk writing. Must be top-level for pickling."""
    chunk_id, gene_indices, gene_data_list, total_cells, output_path = args
    write_expression_chunk(chunk_id, gene_indices, gene_data_list, total_cells, Path(output_path))
    return chunk_id


def _process_obs_column_worker(args):
    """Worker function for parallel obs column processing. Must be top-level for pickling."""
    col_name, col_values, obs_dir, palettes_dir = args
    obs_dir = Path(obs_dir)
    palettes_dir = Path(palettes_dir)

    categorical = is_categorical(col_values, col_name)
    series = pd.Series(col_values)
    series_for_unique = (
        series.astype(str) if series.dtype == "object" else series
    )
    unique_count = int(series_for_unique.nunique(dropna=False))

    metadata_entry = {
        'type': 'categorical' if categorical else 'numerical',
        'unique_values': unique_count
    }

    # Write column values
    if categorical:
        cat_series = series.astype(str)
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

    # Generate palette for categorical
    palette_info = None
    if categorical:
        unique_values = [
            str(v)
            for v in list(cat_series.dropna().unique())
        ]
        palette = generate_color_palette(unique_values)

        palette_file = palettes_dir / f'{col_name}.json'
        palette_file.parent.mkdir(parents=True, exist_ok=True)
        with open(palette_file, 'w') as f:
            json.dump(palette, f, indent=2)

        palette_info = len(unique_values)

    return col_name, metadata_entry, palette_info


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

###############################################################################################################
        # CHANGED FROM ORIGINAL (fuzzy name logic)

        # Check for MERSCOPE files (fuzzy — handles prefixed names like cellpose_metadata.csv)
        files = [f.name for f in input_path.iterdir() if f.is_file()]
        has_merscope = any(
            ('metadata' in f and 'gene' not in f and f.endswith('.csv')) or
            ('cell' in f and 'gene' in f and f.endswith('.csv'))
            for f in files
        )
        print(f"   - Has MERSCOPE files: {has_merscope}")
        print(f"   - Files checked: {files}")
###############################################################################################################

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

    log(f"Loading H5AD file: {input_path.name}", _t_start)
    t_load = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    log(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes (read in {fmt_elapsed(time.perf_counter() - t_load)})", _t_start)

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
    log(f"Loading Xenium folder: {input_path.name}", _t_start)

    # Load cells.csv or cells.csv.gz
    cells_file = None
    for filename in ['cells.csv', 'cells.csv.gz']:
        candidate = input_path / filename
        if candidate.exists():
            cells_file = candidate
            break

    if cells_file is None:
        raise FileNotFoundError("cells.csv or cells.csv.gz not found in Xenium folder")

    log(f"  Loading {cells_file.name}...", _t_start)
    t_load = time.perf_counter()
    if cells_file.suffix == '.gz':
        cells_df = pd.read_csv(cells_file, compression='gzip')
    else:
        cells_df = pd.read_csv(cells_file)

    log(f"  Loaded {len(cells_df):,} cells (read in {fmt_elapsed(time.perf_counter() - t_load)})", _t_start)

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
    log(f"Loading MERSCOPE folder: {input_path.name}", _t_start)

    # Load cell_metadata.csv

###############################################################################################################
    # CHANGED FROM ORIGINAL (fuzzy name logic)
    metadata_file = next(
        (f for f in input_path.iterdir()
        if f.is_file() and 'metadata' in f.name and 'gene' not in f.name and f.suffix == '.csv'),
        None
    )
    if metadata_file is None:
        raise FileNotFoundError("No metadata CSV found in MERSCOPE folder")
###############################################################################################################

    log(f"  Loading {metadata_file.name}...", _t_start)
    t_load = time.perf_counter()
    metadata_df = pd.read_csv(metadata_file)
    log(f"  Loaded {len(metadata_df):,} cells (read in {fmt_elapsed(time.perf_counter() - t_load)})", _t_start)

    # Load cell_by_gene.csv (expression matrix)
###############################################################################################################
    # CHANGED FROM ORIGINAL (fuzzy name logic)
    expr_file = next(
    (f for f in input_path.iterdir()
     if f.is_file() and 'cell' in f.name and 'gene' in f.name and f.suffix == '.csv'),
    None
    )
###############################################################################################################
    # CHANGED FROM ORIGINAL (combine data check)
    expr_matrix = None
    gene_names = None

    if expr_file is not None:
        log(f"  Reading header of {expr_file.name}...", _t_start)
        t_load = time.perf_counter()
        # Only read header + count rows — do NOT load entire matrix into memory
        all_cols = pd.read_csv(expr_file, nrows=0).columns.tolist()
        index_col_name = all_cols[0]  # first column is the cell ID / index
        gene_names = all_cols[1:]     # remaining columns are genes

        # Count rows to validate (read only the index column)
        expr_row_count = 0
        for chunk in pd.read_csv(expr_file, usecols=[index_col_name], chunksize=100_000):
            expr_row_count += len(chunk)

        log(f"  Header read: {len(gene_names)} genes, {expr_row_count:,} rows (in {fmt_elapsed(time.perf_counter() - t_load)})", _t_start)

        if expr_row_count != len(metadata_df):
            raise ValueError(
                f"Row count mismatch: {metadata_file.name} has {len(metadata_df):,} rows, "
                f"{expr_file.name} has {expr_row_count:,} rows. "
                f"Both files must have the same number of rows in the same order."
            )
        log(f"  Expression matrix: {len(gene_names)} genes ({expr_row_count:,} rows, will load per-chunk)", _t_start)
    else:
        log("  cell_by_gene.csv not found, creating placeholder", _t_start)
        gene_names = ['Gene1']
        expr_file = None

    return metadata_df, expr_file, gene_names, index_col_name if expr_file is not None else None
###############################################################################################

def process_dataset(
    input_path: Path,
    output_dir: Path,
    custom_chunk_size: Optional[int] = None,
    data_format: Optional[str] = None,
    num_workers: int = 1
):
    """Main processing function that handles all formats"""
    global _t_start
    _t_start = time.perf_counter()

    log(f"{'='*60}")
    log(f"Processing Spatial Transcriptomics Data")
    log(f"Input: {input_path}")
    log(f"Output: {output_dir}")
    log(f"{'='*60}")

    # Auto-detect format if not specified
    if data_format is None:
        data_format = detect_input_format(input_path)

    log(f"Detected format: {data_format.upper()}", _t_start)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # merscope_expr_file is set only for MERSCOPE (per-chunk CSV loading)
    merscope_expr_file = None

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
                label_col = find_column(cluster_columns_lower, label_candidates)
                if not label_col:
                    print(f"    ⚠ No label column found in {cluster_file.name}")
                    continue
                if len(cluster_df) != len(cells_df):
                    print(f"    ⚠ Row count mismatch: {cluster_file.name} has {len(cluster_df):,} rows, cells has {len(cells_df):,} rows. Skipping.")
                    continue
                non_null = cluster_df[label_col].notna().sum()
                if non_null == 0:
                    print(f"    ⚠ No non-null values in {cluster_file.name}")
                    continue
                obs_columns[label_col] = cluster_df[label_col].fillna('').values
                print(f"    ✓ Imported clusters from {cluster_file.name} using column '{label_col}' ({non_null:,} values, matched positionally)")
                break
            else:
                print("  ⚠ Exhausted cluster files without importing labels")
        if analysis_extract_dir:
            shutil.rmtree(analysis_extract_dir, ignore_errors=True)

    elif data_format == 'merscope':
        metadata_df, merscope_expr_file, gene_names, merscope_index_col = load_merscope_data(input_path)
        dataset_name = input_path.name
        dataset_type = 'merscope'
        embeddings = {}  # No embeddings for MERSCOPE

        if merscope_expr_file is not None:
            expr_matrix = None  # Will be loaded per-chunk from CSV
        else:
            # No expression file found — placeholder
            merscope_expr_file = None
            expr_matrix = np.zeros((len(metadata_df), 1))

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
        exclude_cols = [x_col, y_col, z_col, 'cell_id', 'EntityID']
        for col in metadata_df.columns:
            if col not in exclude_cols:
                obs_columns[col] = metadata_df[col].values

        # Free metadata_df — we've extracted everything we need
        del metadata_df
        gc.collect()
        log("  Freed metadata_df from memory", _t_start)

    else:
        raise ValueError(f"Unknown format: {data_format}")

    # Continue with common processing
    num_cells = len(spatial_coords)
    num_genes = len(gene_names)
    spatial_dims = spatial_coords.shape[1]

    log(f"=== Dataset Summary ===", _t_start)
    log(f"  Cells: {num_cells:,}", _t_start)
    log(f"  Genes: {num_genes:,}", _t_start)
    log(f"  Spatial dimensions: {spatial_dims}D", _t_start)

    # Step 1: Normalize and write spatial coordinates
    log(f"=== STEP 1: Processing spatial coordinates ===", _t_start)
    t_step = time.perf_counter()
    normalized_spatial, scaling_factor = normalize_coordinates(spatial_coords)
    coords_dir = output_dir / 'coords'
    write_coordinate_binary(normalized_spatial, coords_dir / 'spatial.bin.gz')
    log(f"  Spatial coordinates done ({fmt_elapsed(time.perf_counter() - t_step)})", _t_start)

    # Step 2: Process embeddings (if any)
    available_embeddings = []
    if embeddings:
        log(f"=== STEP 2: Processing {len(embeddings)} embedding(s) ===", _t_start)
        for emb_name, emb_coords in embeddings.items():
            t_emb = time.perf_counter()
            normalized_emb, _ = normalize_coordinates(emb_coords)
            write_coordinate_binary(normalized_emb, coords_dir / f'{emb_name}.bin.gz')
            available_embeddings.append(emb_name)
            log(f"  {emb_name}: {emb_coords.shape[1]}D ({fmt_elapsed(time.perf_counter() - t_emb)})", _t_start)
    else:
        log(f"=== STEP 2: No embeddings to process ===", _t_start)

    # Step 3: Process observation columns (clusters)
    log(f"=== STEP 3: Processing {len(obs_columns)} observation column(s) (workers={num_workers}) ===", _t_start)
    t_step = time.perf_counter()
    obs_dir = output_dir / 'obs'
    palettes_dir = output_dir / 'palettes'
    obs_dir.mkdir(parents=True, exist_ok=True)
    palettes_dir.mkdir(parents=True, exist_ok=True)
    obs_metadata = {}
    cluster_count = 0

    if num_workers > 1 and len(obs_columns) > 1:
        # Parallel obs column processing
        worker_args = [
            (col_name, col_values, str(obs_dir), str(palettes_dir))
            for col_name, col_values in obs_columns.items()
        ]
        effective_workers = min(num_workers, len(obs_columns))
        log(f"  Launching {effective_workers} workers for {len(obs_columns)} obs columns...", _t_start)

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            futures = {pool.submit(_process_obs_column_worker, a): a[0] for a in worker_args}
            done_count = 0
            for future in as_completed(futures):
                col_name, metadata_entry, palette_info = future.result()
                obs_metadata[col_name] = metadata_entry
                if palette_info is not None:
                    cluster_count += 1
                done_count += 1
                palette_msg = f" (palette: {palette_info} colors)" if palette_info else ""
                log(f"  [{done_count}/{len(obs_columns)}] {col_name}: {metadata_entry['type']} ({metadata_entry['unique_values']} unique){palette_msg}", _t_start)
    else:
        # Serial fallback
        for col_idx, (col_name, col_values) in enumerate(obs_columns.items(), 1):
            col_name, metadata_entry, palette_info = _process_obs_column_worker(
                (col_name, col_values, str(obs_dir), str(palettes_dir)))
            obs_metadata[col_name] = metadata_entry
            if palette_info is not None:
                cluster_count += 1
            palette_msg = f" (palette: {palette_info} colors)" if palette_info else ""
            log(f"  [{col_idx}/{len(obs_columns)}] {col_name}: {metadata_entry['type']} ({metadata_entry['unique_values']} unique){palette_msg}", _t_start)

    # Write obs metadata
    with open(obs_dir / 'metadata.json', 'w') as f:
        json.dump(obs_metadata, f, indent=2)
    log(f"  Observation columns done ({fmt_elapsed(time.perf_counter() - t_step)})", _t_start)

    # Step 4: Process expression matrix
    chunk_size = determine_chunk_size(num_genes, custom_chunk_size)
    num_chunks = (num_genes + chunk_size - 1) // chunk_size

    log(f"=== STEP 4: Processing expression matrix ({num_genes:,} genes, {num_chunks} chunks, {chunk_size} genes/chunk, workers={num_workers}) ===", _t_start)
    if merscope_expr_file is not None:
        log(f"  Mode: single-pass CSV loading for MERSCOPE", _t_start)
    elif expr_matrix is not None:
        log(f"  Mode: in-memory matrix ({'sparse' if sparse.issparse(expr_matrix) else 'dense'})", _t_start)
    t_step = time.perf_counter()

    # Build expression index
    expr_index = {
        'total_genes': num_genes,
        'num_chunks': num_chunks,
        'chunk_size': chunk_size,
        'genes': []
    }

    expr_dir = output_dir / 'expr'
    expr_dir.mkdir(parents=True, exist_ok=True)

    # --- Extract sparse data for ALL genes ---
    # all_gene_sparse[gene_idx] = (non_zero_indices, non_zero_values)
    all_gene_sparse = [None] * num_genes

    if merscope_expr_file is not None:
        # MERSCOPE: single-pass CSV read — read the file ONCE, extract all genes
        log(f"  Reading entire CSV in one pass...", _t_start)
        t_read = time.perf_counter()
        ROW_CHUNK_SIZE = 500_000
        reader = pd.read_csv(merscope_expr_file, index_col=merscope_index_col, chunksize=ROW_CHUNK_SIZE)

        # Accumulate sparse data per gene across row-batches
        # Each gene gets a list of (offset, indices, values) from each batch
        gene_parts = [[] for _ in range(num_genes)]
        row_offset = 0

        for batch_idx, batch_df in enumerate(reader):
            batch_size = len(batch_df)
            log(f"    Batch {batch_idx + 1}: rows {row_offset:,}-{row_offset + batch_size - 1:,} ({fmt_elapsed(time.perf_counter() - t_read)})", _t_start)

            for gene_idx, gene_name in enumerate(gene_names):
                gene_col = batch_df[gene_name].values
                non_zero_mask = gene_col != 0
                if non_zero_mask.any():
                    local_indices = np.where(non_zero_mask)[0]
                    gene_parts[gene_idx].append((
                        local_indices + row_offset,
                        gene_col[non_zero_mask].astype(np.float32)
                    ))

            row_offset += batch_size
            del batch_df
            gc.collect()

        log(f"  CSV read complete ({row_offset:,} rows, {fmt_elapsed(time.perf_counter() - t_read)})", _t_start)

        # Concatenate parts for each gene
        log(f"  Concatenating sparse data for {num_genes} genes...", _t_start)
        t_concat = time.perf_counter()
        for gene_idx in range(num_genes):
            parts = gene_parts[gene_idx]
            if not parts:
                all_gene_sparse[gene_idx] = (np.array([], dtype=np.uint32), np.array([], dtype=np.float32))
            elif len(parts) == 1:
                all_gene_sparse[gene_idx] = (parts[0][0].astype(np.uint32), parts[0][1])
            else:
                all_gene_sparse[gene_idx] = (
                    np.concatenate([p[0] for p in parts]).astype(np.uint32),
                    np.concatenate([p[1] for p in parts])
                )
        del gene_parts
        gc.collect()
        log(f"  Concatenation done ({fmt_elapsed(time.perf_counter() - t_concat)})", _t_start)

    else:
        # H5AD / Xenium: expression matrix already in memory
        for gene_idx in range(num_genes):
            if sparse.issparse(expr_matrix):
                gene_col = expr_matrix[:, gene_idx].toarray().flatten()
            else:
                gene_col = expr_matrix[:, gene_idx]

            non_zero_mask = gene_col != 0
            non_zero_indices = np.where(non_zero_mask)[0]
            non_zero_values = gene_col[non_zero_mask].astype(np.float32)
            all_gene_sparse[gene_idx] = (non_zero_indices, non_zero_values)

    # --- Build expression index (must be ordered) ---
    for chunk_id in range(num_chunks):
        start_gene = chunk_id * chunk_size
        end_gene = min(start_gene + chunk_size, num_genes)
        for local_idx, gene_idx in enumerate(range(start_gene, end_gene)):
            expr_index['genes'].append({
                'name': gene_names[gene_idx],
                'chunk_id': chunk_id,
                'position_in_chunk': local_idx
            })

    # --- Write chunks (parallel or serial) ---
    log(f"  Writing {num_chunks} chunk files...", _t_start)
    t_write = time.perf_counter()

    if num_workers > 1 and num_chunks > 1:
        # Parallel chunk writing
        effective_workers = min(num_workers, num_chunks)
        log(f"  Launching {effective_workers} workers for {num_chunks} chunks...", _t_start)

        worker_args = []
        for chunk_id in range(num_chunks):
            start_gene = chunk_id * chunk_size
            end_gene = min(start_gene + chunk_size, num_genes)
            gene_indices = list(range(start_gene, end_gene))
            gene_data_list = [all_gene_sparse[gi] for gi in gene_indices]
            chunk_file = str(expr_dir / f'chunk_{chunk_id:05d}.bin.gz')
            worker_args.append((chunk_id, gene_indices, gene_data_list, num_cells, chunk_file))

        with ProcessPoolExecutor(max_workers=effective_workers) as pool:
            futures = {pool.submit(_write_chunk_worker, a): a[0] for a in worker_args}
            done_count = 0
            for future in as_completed(futures):
                cid = future.result()
                done_count += 1
                if done_count % max(1, num_chunks // 10) == 0 or done_count == num_chunks:
                    log(f"    {done_count}/{num_chunks} chunks written ({fmt_elapsed(time.perf_counter() - t_write)})", _t_start)
    else:
        # Serial chunk writing
        for chunk_id in range(num_chunks):
            t_chunk = time.perf_counter()
            start_gene = chunk_id * chunk_size
            end_gene = min(start_gene + chunk_size, num_genes)
            gene_indices = list(range(start_gene, end_gene))
            gene_data_list = [all_gene_sparse[gi] for gi in gene_indices]
            chunk_file = expr_dir / f'chunk_{chunk_id:05d}.bin.gz'
            write_expression_chunk(chunk_id, gene_indices, gene_data_list, num_cells, chunk_file)

            pct = (chunk_id + 1) / num_chunks * 100
            log(f"  [{chunk_id + 1}/{num_chunks}] chunk_{chunk_id:05d} ({end_gene - start_gene} genes, {pct:.0f}%, {fmt_elapsed(time.perf_counter() - t_chunk)})", _t_start)

    log(f"  Chunk writing done ({fmt_elapsed(time.perf_counter() - t_write)})", _t_start)

    del all_gene_sparse
    gc.collect()

    # Write expression index
    with open(expr_dir / 'index.json', 'w') as f:
        json.dump(expr_index, f, indent=2)
    log(f"  Expression matrix done ({fmt_elapsed(time.perf_counter() - t_step)})", _t_start)

    # Step 5: Create manifest
    log(f"=== STEP 5: Creating manifest ===", _t_start)
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
            'coordinates': ['spatial'] + available_embeddings,
            'expression_chunks': num_chunks,
            'observation_columns': list(obs_metadata.keys())
        }
    }

    with open(output_dir / 'manifest.json', 'w') as f:
        json.dump(manifest, f, indent=2)

    log(f"  Manifest written", _t_start)

    # Summary
    log(f"{'='*60}", _t_start)
    log(f"=== COMPLETE ===", _t_start)
    log(f"  Format: {data_format.upper()}", _t_start)
    log(f"  Cells: {num_cells:,}", _t_start)
    log(f"  Genes: {num_genes:,}", _t_start)
    log(f"  Chunks: {num_chunks}", _t_start)
    log(f"  Obs columns: {len(obs_metadata)}", _t_start)
    log(f"  Output: {output_dir}", _t_start)
    log(f"  Total time: {fmt_elapsed(time.perf_counter() - _t_start)}", _t_start)
    log(f"{'='*60}", _t_start)


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
    parser.add_argument('--workers', type=int, default=1,
                        help='Number of parallel workers for chunk writing and obs processing (default: 1)')

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"❌ Error: Input path not found: {args.input}")
        sys.exit(1)

    # Process
    try:
        process_dataset(args.input, args.output, args.chunk_size, args.format, args.workers)
    except Exception as e:
        print(f"\n❌ Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
