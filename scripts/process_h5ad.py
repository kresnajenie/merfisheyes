#!/usr/bin/env python3
"""
H5AD to Chunked Binary Format Processor
Converts H5AD files to the binary chunk format expected by ChunkedDataAdapter.ts

Usage:
    python process_h5ad.py input.h5ad output_folder/ [--chunk-size GENES_PER_CHUNK]

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
from typing import Dict, List, Tuple, Any, Optional

import anndata as ad
import numpy as np
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
    num_cells = len(indices)  # This should be total cells in dataset
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

    Format:
    - Header: [version: uint32, num_genes: uint32, chunk_id: uint32, total_cells: uint32]
    - Gene table: [gene_index: uint32, data_offset: uint32, data_size: uint32,
                   uncompressed_size: uint32, num_non_zero: uint32] * num_genes
    - Gene data: sparse format for each gene
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


def process_h5ad(input_path: Path, output_dir: Path, custom_chunk_size: Optional[int] = None):
    """Main processing function"""

    print(f"\n{'='*60}")
    print(f"Processing H5AD: {input_path.name}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}\n")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load H5AD file
    print("📂 Loading H5AD file...")
    adata = ad.read_h5ad(input_path)
    print(f"  ✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes\n")

    # Extract metadata
    num_cells = adata.n_obs
    num_genes = adata.n_vars
    gene_names = list(adata.var_names)

    # 1. Process spatial coordinates
    print("📍 Processing spatial coordinates...")
    spatial_coords = None
    spatial_key = None

    # Try to find spatial coordinates
    for key in ['X_spatial', 'spatial']:
        if key in adata.obsm:
            spatial_coords = adata.obsm[key]
            spatial_key = key
            break

    if spatial_coords is None:
        print("  ⚠ No spatial coordinates found in obsm, checking obs columns...")
        # Try to extract from obs columns (matching TypeScript H5adAdapter logic)
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
            spatial_key = 'from_obs'

    if spatial_coords is None:
        raise ValueError("No spatial coordinates found! Check your H5AD file structure.")

    # Normalize coordinates
    normalized_spatial, scaling_factor = normalize_coordinates(spatial_coords)
    spatial_dims = normalized_spatial.shape[1]

    # Write spatial coordinates
    coords_dir = output_dir / 'coords'
    write_coordinate_binary(normalized_spatial, coords_dir / 'spatial.bin.gz')

    # 2. Process embeddings (UMAP, etc.)
    print("\n🗺️  Processing embeddings...")
    available_embeddings = []
    MAX_EMBEDDING_DIMS = 3  # Maximum dimensions to save for embeddings
    for key in adata.obsm.keys():
        if key not in ['X_spatial', 'spatial'] and key.startswith('X_'):
            embedding_name = key[2:].lower()  # Remove 'X_' prefix
            embedding_coords = adata.obsm[key]

            # Limit to first 3 dimensions
            if embedding_coords.shape[1] > MAX_EMBEDDING_DIMS:
                embedding_coords = embedding_coords[:, :MAX_EMBEDDING_DIMS]
                print(f"  ℹ {embedding_name}: Limiting to first {MAX_EMBEDDING_DIMS} dimensions (from {adata.obsm[key].shape[1]})")

            normalized_emb, _ = normalize_coordinates(embedding_coords)
            write_coordinate_binary(normalized_emb, coords_dir / f'{embedding_name}.bin.gz')
            available_embeddings.append(embedding_name)

    if not available_embeddings:
        print("  ℹ No embeddings found")

    # 3. Process clusters
    print("\n🏷️  Processing cluster columns...")
    obs_dir = output_dir / 'obs'
    palettes_dir = output_dir / 'palettes'

    obs_metadata = {}
    cluster_count = 0

    # Get all columns from obs (excluding index)
    obs_columns = [col for col in adata.obs.columns if col != 'index']

    for col in obs_columns:
        values = adata.obs[col].values
        categorical = is_categorical(values)

        # Store metadata
        obs_metadata[col] = {
            'type': 'categorical' if categorical else 'numerical',
            'unique_values': int(len(np.unique(values)))
        }

        # Write column values as JSON
        values_list = [str(v) if categorical else float(v) for v in values]
        obs_file = obs_dir / f'{col}.json.gz'
        obs_file.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(obs_file, 'wt', encoding='utf-8') as f:
            json.dump(values_list, f)

        print(f"  ✓ {col}: {obs_metadata[col]['type']} ({obs_metadata[col]['unique_values']} unique)")

        # Generate palette for categorical columns
        if categorical:
            cluster_count += 1
            unique_values = [str(v) for v in np.unique(values)]
            palette = generate_color_palette(unique_values)

            palette_file = palettes_dir / f'{col}.json'
            palette_file.parent.mkdir(parents=True, exist_ok=True)
            with open(palette_file, 'w') as f:
                json.dump(palette, f, indent=2)

            print(f"    → Palette generated with {len(unique_values)} colors")

    # Write obs metadata
    with open(obs_dir / 'metadata.json', 'w') as f:
        json.dump(obs_metadata, f, indent=2)

    # 4. Process expression matrix
    print("\n🧬 Processing expression matrix...")

    # Determine chunk size
    chunk_size = determine_chunk_size(num_genes, custom_chunk_size)
    num_chunks = (num_genes + chunk_size - 1) // chunk_size  # Ceiling division

    print(f"  Chunk size: {chunk_size} genes/chunk")
    print(f"  Total chunks: {num_chunks}")

    # Get expression matrix
    if sparse.issparse(adata.X):
        expr_matrix = adata.X.tocsc()  # Column-sparse for efficient gene access
        print(f"  Matrix format: sparse CSC")
    else:
        expr_matrix = adata.X
        print(f"  Matrix format: dense")

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
        chunk_gene_count = end_gene - start_gene

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

    # 5. Create manifest
    print("\n📋 Creating manifest...")
    manifest = {
        'version': '1.0',
        'dataset_id': 'local_dataset',  # Placeholder
        'name': input_path.stem,
        'type': 'h5ad',
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
            'created_by': 'process_h5ad.py'
        },
        'files': {
            'coordinates': ['spatial'] + available_embeddings,
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
    print(f"{'='*60}")
    print(f"\nOutput structure:")
    print(f"  {output_dir}/")
    print(f"  ├── manifest.json")
    print(f"  ├── coords/")
    print(f"  │   ├── spatial.bin.gz")
    for emb in available_embeddings:
        print(f"  │   └── {emb}.bin.gz")
    print(f"  ├── expr/")
    print(f"  │   ├── index.json")
    print(f"  │   └── chunk_*.bin.gz ({num_chunks} files)")
    print(f"  ├── obs/")
    print(f"  │   ├── metadata.json")
    print(f"  │   └── *.json.gz ({len(obs_metadata)} files)")
    print(f"  └── palettes/")
    print(f"      └── *.json ({cluster_count} files)")
    print(f"\n📊 Statistics:")
    print(f"  Cells: {num_cells:,}")
    print(f"  Genes: {num_genes:,}")
    print(f"  Spatial dimensions: {spatial_dims}D")
    print(f"  Chunk size: {chunk_size} genes/chunk")
    print(f"  Total chunks: {num_chunks}")
    print(f"  Categorical columns: {cluster_count}")
    print(f"\n")


def main():
    parser = argparse.ArgumentParser(
        description='Convert H5AD files to chunked binary format for ChunkedDataAdapter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-determine chunk size based on gene count
  python process_h5ad.py data.h5ad output/

  # Use single gene per chunk
  python process_h5ad.py data.h5ad output/ --chunk-size 1

  # Use 100 genes per chunk
  python process_h5ad.py data.h5ad output/ --chunk-size 100
        """
    )

    parser.add_argument('input', type=Path, help='Input H5AD file')
    parser.add_argument('output', type=Path, help='Output directory for chunked data')
    parser.add_argument('--chunk-size', type=int, default=None,
                        help='Genes per chunk (default: auto-determined based on gene count, use 1 for single gene chunks)')

    args = parser.parse_args()

    # Validate input
    if not args.input.exists():
        print(f"❌ Error: Input file not found: {args.input}")
        sys.exit(1)

    if not args.input.suffix == '.h5ad':
        print(f"⚠ Warning: Input file doesn't have .h5ad extension")

    # Process
    try:
        process_h5ad(args.input, args.output, args.chunk_size)
    except Exception as e:
        print(f"\n❌ Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
