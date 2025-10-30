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
import struct
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional, Union

import numpy as np
import pandas as pd
from scipy import sparse


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


def is_categorical(values: np.ndarray) -> bool:
    """Determine if data is categorical (≤100 unique values)"""
    unique_count = len(np.unique(values))
    return unique_count <= 100


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
    if input_path.is_file():
        if input_path.suffix == '.h5ad':
            return 'h5ad'
        else:
            raise ValueError(f"Unknown file format: {input_path}")

    elif input_path.is_dir():
        # Check for Xenium files
        xenium_files = ['cells.csv', 'cells.csv.gz', 'cell_feature_matrix']
        has_xenium = any((input_path / f).exists() for f in xenium_files)

        # Check for MERSCOPE files
        merscope_files = ['cell_metadata.csv', 'cell_by_gene.csv']
        has_merscope = any((input_path / f).exists() for f in merscope_files)

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
        features_file = None
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
            gene_names = features_df[0].tolist()
            print(f"  ✓ Loaded {len(gene_names)} genes from features file")

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
        for key in adata.obsm.keys():
            if key not in ['X_spatial', 'spatial'] and key.startswith('X_'):
                embedding_name = key[2:].lower()  # Remove 'X_' prefix
                embeddings[embedding_name] = adata.obsm[key]

    elif data_format == 'xenium':
        cells_df, expr_matrix, gene_names = load_xenium_data(input_path)
        dataset_name = input_path.name
        dataset_type = 'xenium'
        embeddings = {}  # No embeddings for Xenium

        # Extract spatial coordinates
        coord_candidates_x = ['x_centroid', 'cell_centroid_x', 'center_x']
        coord_candidates_y = ['y_centroid', 'cell_centroid_y', 'center_y']
        coord_candidates_z = ['z_centroid', 'cell_centroid_z', 'center_z']

        x_col = next((c for c in coord_candidates_x if c in cells_df.columns), None)
        y_col = next((c for c in coord_candidates_y if c in cells_df.columns), None)
        z_col = next((c for c in coord_candidates_z if c in cells_df.columns), None)

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
        exclude_cols = [x_col, y_col, z_col, 'cell_id', 'id']
        for col in cells_df.columns:
            if col not in exclude_cols and col not in gene_names:
                obs_columns[col] = cells_df[col].values

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

    # 2. Process observation columns (clusters)
    print("\n🏷️  Processing observation columns...")
    obs_dir = output_dir / 'obs'
    palettes_dir = output_dir / 'palettes'
    obs_metadata = {}
    cluster_count = 0

    for col_name, col_values in obs_columns.items():
        categorical = is_categorical(col_values)

        obs_metadata[col_name] = {
            'type': 'categorical' if categorical else 'numerical',
            'unique_values': int(len(np.unique(col_values)))
        }

        # Write column values
        values_list = [str(v) if categorical else float(v) for v in col_values]
        obs_file = obs_dir / f'{col_name}.json.gz'
        obs_file.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(obs_file, 'wt', encoding='utf-8') as f:
            json.dump(values_list, f)

        print(f"  ✓ {col_name}: {obs_metadata[col_name]['type']} ({obs_metadata[col_name]['unique_values']} unique)")

        # Generate palette for categorical
        if categorical:
            cluster_count += 1
            unique_values = [str(v) for v in np.unique(col_values)]
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
            'available_embeddings': [],
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
