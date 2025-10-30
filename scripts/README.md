# Spatial Transcriptomics Data Processing Scripts

This directory contains Python scripts for preprocessing spatial transcriptomics datasets that are too large for browser processing.

## Overview

The `process_spatial_data.py` script converts **H5AD**, **Xenium**, and **MERSCOPE** datasets into a chunked binary format that can be loaded directly into the MERFISH Eyes viewer without requiring browser-based processing. This is essential for very large datasets (>500K cells) that would exceed browser memory limits.

## Why Use This Script?

### Problems with Large Datasets in Browser:
- **Memory Limits**: Browsers typically have 2-4GB memory limits
- **UI Freezing**: Even with web workers, very large files can cause issues
- **Slow Processing**: Parsing multi-GB H5AD files in JavaScript is slow

### Benefits of Python Preprocessing:
- ✅ **No Browser Memory Limits**: Process datasets of any size
- ✅ **Faster Processing**: Python with numpy/scipy is optimized for large data
- ✅ **Instant Loading**: Browser just reads pre-processed chunks
- ✅ **On-Demand Gene Loading**: Only load expression data for selected genes

## Installation

### Requirements:
```bash
pip install anndata numpy scipy pandas
```

### Dependencies:
- `anndata` - Reading H5AD files (required for H5AD format)
- `pandas` - Reading CSV files (required for Xenium/MERSCOPE)
- `numpy` - Numerical operations
- `scipy` - Sparse matrix support
- Standard library: `gzip`, `json`, `struct`, `pathlib`

## Usage

### H5AD Files:
```bash
python scripts/process_spatial_data.py input.h5ad output_folder/
```

### Xenium Folders:
```bash
python scripts/process_spatial_data.py xenium_output/ output_folder/
```

### MERSCOPE Folders:
```bash
python scripts/process_spatial_data.py merscope_output/ output_folder/
```

The script **automatically detects** the input format based on file/folder structure.

### Chunk Size Options:

**Auto-determined (default):**
- <100 genes: 50 genes/chunk
- <500 genes: 100 genes/chunk
- <2000 genes: 200 genes/chunk
- <10000 genes: 500 genes/chunk
- ≥10000 genes: 1000 genes/chunk

**Single Gene Per Chunk (For Very Large Datasets):**
```bash
python scripts/process_spatial_data.py large_data.h5ad output_folder/ --chunk-size 1
```

**Custom Chunk Size:**
```bash
python scripts/process_spatial_data.py data.h5ad output_folder/ --chunk-size 100
```

## Output Structure

The script creates a folder with this structure:

```
output_folder/
├── manifest.json                # Dataset metadata
├── coords/
│   ├── spatial.bin.gz          # Spatial coordinates (normalized to [-1, 1])
│   ├── umap.bin.gz             # UMAP embedding (max 3D)
│   └── pca.bin.gz              # PCA embedding (max 3D, truncated from 50)
├── expr/
│   ├── index.json              # Gene → chunk mapping
│   ├── chunk_00000.bin.gz      # Expression data chunks (sparse format)
│   ├── chunk_00001.bin.gz
│   └── ...
├── obs/
│   ├── metadata.json           # Column types (categorical/numerical)
│   ├── leiden.json.gz          # Cluster assignments
│   └── ...                     # Other observation columns
└── palettes/
    └── leiden.json             # Color mappings (categorical columns only)
```

**Note on Embeddings**: All embeddings (PCA, UMAP, etc.) are automatically limited to **first 3 dimensions** to reduce file size:
- PCA with 50 components → saves only first 3 (~94% size reduction)
- UMAP typically already 2D or 3D (no change)
- Embeddings >3D are truncated during processing

## Loading in MERFISH Eyes

### Step 1: Process Your Dataset
```bash
# H5AD file
python scripts/process_spatial_data.py my_large_dataset.h5ad processed_data/

# Xenium folder
python scripts/process_spatial_data.py xenium_output/ processed_data/

# MERSCOPE folder
python scripts/process_spatial_data.py merscope_output/ processed_data/
```

### Step 2: Upload Folder to MERFISH Eyes
1. Go to the MERFISH Eyes homepage
2. Click on **"Chunked Folder"** upload box (4th box on homepage)
3. **Select the entire `processed_data` folder**
4. The app will load it instantly (1-3 seconds) - no browser processing!
5. Click "Upload & Save" button
6. Upload modal shows:
   - ✅ Chunk settings hidden (already processed)
   - ✅ "Pre-chunked and ready for upload" message
   - Just enter email and dataset name
7. Upload to S3 and receive email notification
8. View at `/viewer/[datasetId]`

**Benefits of Pre-Chunked Upload**:
- **No browser processing** - instant loading
- **No memory limits** - Python handles all processing
- **Faster uploads** - no client-side chunking needed
- **Large dataset support** - handle 500K+ cells easily

## Binary Format Details

### Coordinate Files (`.bin.gz`)
```
Header:
  - num_points: uint32
  - dimensions: uint32
Data:
  - coordinates: float32[] (flattened row-major)
```

### Expression Chunks (`.bin.gz`)
```
Header (16 bytes):
  - version: uint32 (always 1)
  - num_genes: uint32
  - chunk_id: uint32
  - total_cells: uint32

Gene Table (24 bytes per gene):
  - gene_index: uint32
  - data_offset: uint32
  - data_size: uint32
  - uncompressed_size: uint32
  - num_non_zero: uint32
  - reserved: uint32

Gene Data (per gene, sparse format):
  - num_cells: uint32
  - num_non_zero: uint32
  - indices: uint32[] (cell indices with non-zero values)
  - values: float32[] (expression values)
```

### Cluster Files (`.json.gz`)
Simple JSON array with one value per cell:
```json
["Neuron", "Astrocyte", "Neuron", ...]
```

### Palette Files (`.json`)
Mapping of cluster values to hex colors:
```json
{
  "Neuron": "#1f77b4",
  "Astrocyte": "#ff7f0e",
  ...
}
```

## Color Palette

The script uses the exact same 40-color palette as the browser visualization:
- Bright, distinct colors optimized for black backgrounds
- Consistent across browser and Python processing
- Cycles through palette for >40 categories

## Compatibility

The output format is **100% compatible** with `ChunkedDataAdapter.ts` which supports both:
- **Remote mode**: Loading from S3 via presigned URLs
- **Local mode**: Loading from local File objects (this script's output)

Both modes use identical binary parsing logic, ensuring consistency.

## Performance

### Example Benchmarks:

**100K cells, 2K genes (H5AD ~800MB):**
- Python processing: ~30-60 seconds
- Browser loading: ~2-3 seconds (just reads metadata)
- Gene selection: <100ms (loads chunk on-demand)

**1M cells, 5K genes (H5AD ~8GB):**
- Python processing: ~5-10 minutes
- Browser loading: ~3-5 seconds
- Gene selection: ~500ms-1s (first time), <50ms (cached)

## Troubleshooting

### Error: "No spatial coordinates found"
- Your H5AD must have `obsm['X_spatial']` or spatial coordinates in `obs` columns
- Check your H5AD structure: `adata.obsm.keys()` and `adata.obs.columns`

### Error: "File missing webkitRelativePath"
- When uploading to browser, make sure to select the **entire folder**, not individual files
- Use browser's folder selection feature (available in Chrome/Edge)

### Large Memory Usage During Processing
- The Python script loads the entire expression matrix into memory
- For >10M cells, you may need 32GB+ RAM
- Consider using a server or cloud instance for very large datasets

## Advanced Usage

### Processing Multiple Datasets
```bash
for file in *.h5ad; do
    output_dir="${file%.h5ad}_processed"
    python scripts/process_h5ad.py "$file" "$output_dir/"
done
```

### Custom Chunk Size Based on Dataset Size
```python
# In Python
import anndata as ad

adata = ad.read_h5ad('data.h5ad')
num_genes = adata.n_vars

# Determine chunk size
if num_genes > 10000:
    chunk_size = 1  # Single gene per chunk for huge datasets
elif num_genes > 5000:
    chunk_size = 100
else:
    chunk_size = 500

# Run with custom chunk size
# python scripts/process_h5ad.py data.h5ad output/ --chunk-size {chunk_size}
```

## Future Enhancements

Potential improvements:
- Streaming processing to reduce memory usage
- Parallel chunk generation
- Support for other input formats (Seurat, Loom, etc.)
- Direct S3 upload after processing
- Progress bar for long-running operations

## Questions?

For issues or questions about the Python processing script, please check:
1. The main `CLAUDE.md` documentation
2. The `DATA_PROCESSING_ARCHITECTURE.md` for technical details
3. GitHub issues
