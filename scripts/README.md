# Spatial Transcriptomics Data Processing Scripts

This directory contains Python scripts for preprocessing spatial transcriptomics datasets that are too large for browser processing.

## Overview

This directory contains two types of processing scripts:

1. **`process_spatial_data.py`** - Converts **single cell data** (H5AD, Xenium, MERSCOPE) into chunked binary format
2. **`process_single_molecule.py`** - Converts **single molecule data** (parquet, CSV) into S3-compatible binary format

### Single Cell Data Processing

The `process_spatial_data.py` script converts **H5AD**, **Xenium**, and **MERSCOPE** datasets into a chunked binary format that can be loaded directly into the MERFISH Eyes viewer without requiring browser-based processing. This is essential for very large datasets (>500K cells) that would exceed browser memory limits.

### Single Molecule Data Processing

The `process_single_molecule.py` script converts **parquet** and **CSV** single molecule datasets into the same binary format used by the S3 storage. This allows processing large molecule datasets (millions of molecules) in Python and uploading directly to S3, bypassing browser memory limits.

> **These same two scripts run the "Upload & process on server" path.** The
> server-side ingestion worker ([`worker/entrypoint.py`](../worker/entrypoint.py))
> copies them into its Docker image and invokes them exactly as documented here —
> single-cell with `--chunk-size`/`--mmc-csv`, single-molecule with the
> user-confirmed `--dataset-type custom --gene-col …` mapping. So the flags below
> are also the contract the server uses. See
> [docs/server-side-ingestion/](../docs/server-side-ingestion/README.md).

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

### Requirements for Single Cell Processing:
```bash
pip install anndata numpy scipy pandas
```

### Requirements for Single Molecule Processing:
```bash
pip install numpy pandas
# pyarrow is optional (only needed for .parquet input, not for .csv)
pip install pyarrow
```

### Dependencies:
**Single Cell Scripts:**
- `anndata` - Reading H5AD files (required for H5AD format)
- `pandas` - Reading CSV files (required for Xenium/MERSCOPE)
- `numpy` - Numerical operations
- `scipy` - Sparse matrix support
- Standard library: `gzip`, `json`, `struct`, `pathlib`

**Single Molecule Scripts:**
- `numpy` - Numerical operations and coordinate normalization
- `pandas` - Reading CSV files
- `pyarrow` - Reading parquet files
- Standard library: `gzip`, `json`, `re`, `pathlib`

## Usage

### Single Cell Data Processing

#### H5AD Files:
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

#### Chunk Size Options:

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

### Single Molecule Data Processing

#### Basic Usage:

**Parquet files (MERSCOPE - default):**
```bash
python scripts/process_single_molecule.py molecules.parquet output_folder/
```

**CSV files (MERSCOPE):**
```bash
python scripts/process_single_molecule.py molecules.csv output_folder/
```

**Xenium format:**
```bash
python scripts/process_single_molecule.py transcripts.parquet output_folder/ --dataset-type xenium
```

#### Custom Column Names:

If your dataset has non-standard column names:
```bash
python scripts/process_single_molecule.py data.parquet output_folder/ \
  --gene-col "my_gene_column" \
  --x-col "my_x_column" \
  --y-col "my_y_column" \
  --z-col "my_z_column"
```

For 2D datasets (no Z coordinate):
```bash
python scripts/process_single_molecule.py data.parquet output_folder/ \
  --gene-col "gene_name" \
  --x-col "x_pos" \
  --y-col "y_pos"
```

#### Splitting molecules by cell assignment (`--cell-id-col`):

Off by default — every molecule for a gene goes into one `GENE.bin.gz`. Opt in
with `--cell-id-col` to split each gene's **unassigned** (not-in-a-cell) molecules
into a separate `GENE_uuuuuuuuuu.bin.gz` sibling, so the viewer can show
assigned vs. unassigned separately:

```bash
python scripts/process_single_molecule.py transcripts.parquet output_folder/ \
  --dataset-type xenium \
  --cell-id-col cell_id
```

**"Unassigned" is detected across platforms**, not just a single sentinel value:

| Platform | Unassigned cell-id | 
|----------|--------------------|
| MERSCOPE (Vizgen) | integer `-1` |
| Xenium (10x) | string `UNASSIGNED` |
| CosMx (NanoString) | integer `0` (assigned IDs are ≥ 1) |

The rule: any numeric value `≤ 0` (covers `-1` and `0`), plus `NaN`/blank and the
common string tokens (`unassigned`, `none`, `nan`, `na`, `null`). Any other
non-numeric string (e.g. a real Xenium cell id like `aaabbbcc-1`) counts as
**assigned**.

#### Dataset Types:

The script supports three preset column mappings:

**MERSCOPE (default):**
- Gene: `gene`
- X: `global_x`
- Y: `global_y`
- Z: `global_z`

**Xenium:**
- Gene: `feature_name`
- X: `x_location`
- Y: `y_location`
- Z: `z_location`

**Custom:**
- Same as Xenium, but use `--gene-col`, `--x-col`, etc. to override

## Output Structure

### Single Cell Output Structure

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

### Single Molecule Output Structure

The script creates a folder matching the S3 storage structure:

```
output_folder/
├── manifest.json.gz            # Dataset metadata
└── genes/
    ├── GENE1.bin.gz           # Gene coordinates (Float32Array, gzipped)
    ├── GENE2.bin.gz
    └── ...
```

**Manifest Format** (`manifest.json.gz`):
```json
{
  "uniqueGenes": ["GENE1", "GENE2", ...],
  "dimensions": 3,
  "scalingFactor": 1234.56,
  "totalMolecules": 21000000,
  "geneCount": 500
}
```

**Gene File Format** (`genes/{geneName}.bin.gz`):
- Gzipped Float32Array of normalized coordinates
- Flat array: `[x1, y1, z1, x2, y2, z2, ...]` for 3D
- Or: `[x1, y1, 0, x2, y2, 0, ...]` for 2D (z filled with zeros)
- Coordinates already normalized to [-1, 1] range
- Each coordinate is 4 bytes (Float32)

**Gene Filtering**: The script automatically filters out:
- Negative controls (`NegativeControl`, `NegCtrl`, etc.)
- Unassigned genes
- Deprecated genes
- Codewords
- Blanks

**Filename Sanitization**: Gene names with special characters are sanitized:
- `GENE-1` → `GENE_1.bin.gz`
- `CD8+` → `CD8_.bin.gz`

## Loading in MERFISH Eyes

### Single Cell Data

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

### Single Molecule Data

#### Step 1: Process Your Dataset
```bash
# MERSCOPE parquet
python scripts/process_single_molecule.py molecules.parquet processed_molecules/

# Xenium parquet
python scripts/process_single_molecule.py transcripts.parquet processed_transcripts/ --dataset-type xenium

# CSV file
python scripts/process_single_molecule.py molecules.csv processed_molecules/
```

#### Step 2: Upload to S3

**Option A: Manual S3 Upload**
```bash
# Upload the entire output folder to S3
aws s3 cp processed_molecules/ s3://your-bucket/datasets/{datasetId}/ --recursive

# The structure should be:
# s3://your-bucket/datasets/{datasetId}/manifest.json.gz
# s3://your-bucket/datasets/{datasetId}/genes/GENE1.bin.gz
# s3://your-bucket/datasets/{datasetId}/genes/GENE2.bin.gz
# ...
```

**Option B: Use the MERFISH Eyes Upload Flow**
1. Go to MERFISH Eyes homepage
2. Select the processed folder for upload
3. The app will upload it to S3 automatically
4. View at `/sm-viewer/[datasetId]`

**Benefits of Python Processing**:
- **No browser memory limits** - handle millions of molecules
- **Faster processing** - Python is optimized for large data
- **Lazy loading** - genes loaded on-demand when selected
- **Consistent format** - exact same structure as browser processing

## Preparing a folder of parts (`prepare_folder.py`)

An AnnData table (`.h5ad`, or `write_zarr` output in Zarr v2 or v3), a
molecule folder from `process_single_molecule.py` and an OME-Zarr image,
side by side in one folder:

```bash
python scripts/prepare_folder.py sample/              # --cluster: Leiden + DEGs when the table has no categorical column
python scripts/prepare_folder.py sample/ --url https://<bucket>.s3.<region>.amazonaws.com/<path>/sample
```

Writes one SpatialData store `sample.zarr` (Zarr v3): the table in
`tables/table` with `layers/X_csc` when `X` is not CSC, the image as
`images/<name>` written by the `spatialdata` library from the OME-Zarr
folder, the molecules as `points/molecules` (a tiled points element,
`legacy_to_points.py`); the image and molecule folders move to `sources/`;
`mapping.json` names the store; with `--cluster` the Leiden clusters go
into the table as `obs/leiden` and `annotations/` holds the DEG stats.
Then `validate_store.py`. `--expanded` keeps the image and molecule
folders beside the store, named in `mapping.json` (relative names, absolute
with `--url`), and checks that they sit where the cells are (a `WARN`
otherwise). A store that is there already is repaired in place where the
viewer cannot read it as written: shapes stored as a directory of parquet
parts or in the geoarrow encoding are rewritten as one WKB
`shapes.parquet` (the original kept as `shapes.parquet.parts`), image
levels in chunks over 1024² are rewritten in 1024² chunks with the same
codecs, and `--molecules-from-points` turns the store's parquet points
element into the tiled one. Dependencies for all of this: `pip install -r
scripts/requirements-prepare.txt spatialdata` (`geopandas` comes with it).
Tests: `tests/python/test_prepare_folder.py`.

### A Xenium export (`prepare_xenium.py`)

```bash
pip install spatialdata-io          # Python 3.12 or newer, on top of requirements-prepare.txt
python scripts/prepare_xenium.py Xenium_outs/ sample/ --workers 8    # --cluster, --labels, --keep-z, --url
```

Reads the export with `spatialdata_io.xenium` (boundaries, transcripts,
the morphology focus image with 1024² chunks and a pyramid, the cell
table; no raster labels unless `--labels`, the table then annotating the
polygons), writes the store, converts `points/transcripts` into
`sample/molecules/` (flat; `--keep-z` keeps z), then runs
`prepare_folder.py`'s checks. Test: `MERFISHEYES_XENIUM=/path/to/outs
python -m pytest tests/python/test_prepare_xenium.py` (skipped without it).

### A MERSCOPE region (`prepare_merscope.py`)

```bash
pip install spatialdata                 # on top of requirements-prepare.txt
python scripts/prepare_merscope.py region_R1/ sample/ --workers 8    # --z 3, --cluster, --keep-z, --url
```

Builds the store with the `spatialdata` library from the instrument's
files: `cell_metadata.csv` + `cell_by_gene.csv` (blanks dropped, the counts
aligned to the metadata's cell ids, written CSC; `cell_categories.csv` as
categorical obs columns) and, when present, the vpt `cell_boundaries.parquet`
as `shapes/cell_boundaries` (each cell's polygon on one z-plane, `--z`, the
largest part of a MultiPolygon). `detected_transcripts.csv` (or `.parquet`)
becomes `sample/molecules/` through `process_single_molecule.py`. Then
`prepare_folder.py`'s checks. Not through spatialdata-io: its MERSCOPE
reader needs the `images/` folder most exports lack. The mosaic images are
not written (tens of GB; `image_to_ome_zarr.py` reads an image whole).
Test: `tests/python/test_prepare_merscope.py` (a MERSCOPE-shaped export made
in the test).

### A Visium or Visium HD output (`prepare_visium.py`)

```bash
pip install spatialdata-io
python scripts/prepare_visium.py outs/ sample/ --fullres image.tif     # --bins 8 16 (HD), --cluster, --url
```

Reads the Space Ranger output with `spatialdata_io.visium` (or `visium_hd`
when `binned_outputs/` is there): the spots or bins as shapes, the tables
(the CSC copy added), the hires and lowres images and, with `--fullres`,
the microscope image in 1024² chunks with a pyramid. A Visium store has
no shared `global` system; the viewer targets the spots' own. For Visium
HD the 8 µm table is made the default (`metadata_general["Default Table"]`)
and 2 µm bins are not read unless asked. Test:
`tests/python/test_prepare_visium.py` (a Space Ranger-shaped output made in
the test; the HD case runs with `MERFISHEYES_VISIUM_HD=/path/to/outs`).

### A TIFF as an OME-Zarr image (`image_to_ome_zarr.py`)

```bash
python scripts/image_to_ome_zarr.py dapi.tif dapi.zarr --um-per-pixel 0.108 --origin 0 0 --channels DAPI
```

One plane, a `(c, y, x)` stack or an RGB image, uint8 or uint16, written as
OME-Zarr 0.5 in Zarr v3 with axes `c, y, x`, 1024² chunks and a pyramid
halved until a level fits one chunk. `--um-per-pixel` and `--origin` (the
centre of pixel 0, 0, in the cells' units) are what place it under the cells.
The image is held in memory while it is written.

## Binary Format Details

### Single Cell Binary Formats

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

### Single Molecule Binary Formats

#### Manifest File (`manifest.json.gz`)
Gzipped JSON with dataset metadata:
```json
{
  "uniqueGenes": ["GENE1", "GENE2", ...],
  "dimensions": 2 | 3,
  "scalingFactor": 1234.56,
  "totalMolecules": 21000000,
  "geneCount": 500
}
```

#### Gene Files (`genes/{geneName}.bin.gz`)
Gzipped binary file containing normalized coordinates as Float32Array:
```
Binary Format:
  - Flat Float32 array: [x1, y1, z1, x2, y2, z2, ...]
  - Each coordinate is 4 bytes (Float32)
  - Coordinates normalized to [-1, 1] range
  - For 2D: z coordinates are 0
  - Little-endian byte order
```

**Example**: Gene with 3 molecules in 3D:
```
[x1, y1, z1, x2, y2, z2, x3, y3, z3]
```
Byte size: `3 molecules × 3 dimensions × 4 bytes = 36 bytes` (before gzip)

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

#### Single Cell Processing

**100K cells, 2K genes (H5AD ~800MB):**
- Python processing: ~30-60 seconds
- Browser loading: ~2-3 seconds (just reads metadata)
- Gene selection: <100ms (loads chunk on-demand)

**1M cells, 5K genes (H5AD ~8GB):**
- Python processing: ~5-10 minutes
- Browser loading: ~3-5 seconds
- Gene selection: ~500ms-1s (first time), <50ms (cached)

#### Single Molecule Processing

**21M molecules, 500 genes (Parquet ~2GB):**
- Python processing: ~2-3 minutes
- Output size: ~1.5GB (compressed)
- Browser loading: ~1-2 seconds (manifest only)
- Gene selection: ~200-500ms (first time), instant (cached)

**100M molecules, 1000 genes (Parquet ~10GB):**
- Python processing: ~10-15 minutes
- Output size: ~7GB (compressed)
- Browser loading: ~2-3 seconds (manifest only)
- Gene selection: ~500ms-1s (first time), instant (cached)

## Troubleshooting

### Single Cell Processing

#### Error: "No spatial coordinates found"
- Your H5AD must have `obsm['X_spatial']` or spatial coordinates in `obs` columns
- Check your H5AD structure: `adata.obsm.keys()` and `adata.obs.columns`

#### Error: "File missing webkitRelativePath"
- When uploading to browser, make sure to select the **entire folder**, not individual files
- Use browser's folder selection feature (available in Chrome/Edge)

#### Large Memory Usage During Processing
- The Python script loads the entire expression matrix into memory
- For >10M cells, you may need 32GB+ RAM
- Consider using a server or cloud instance for very large datasets

### Single Molecule Processing

#### Error: "Missing required columns"
- Check your column names match the dataset type
- Use `--gene-col`, `--x-col`, `--y-col`, `--z-col` for custom columns
- For CSV: ensure headers are present in first row
- For Parquet: check schema with `pyarrow.parquet.read_schema()`

#### Error: Memory issues with large files
- Parquet files are processed in memory
- For >50M molecules, you may need 16GB+ RAM
- Consider splitting file into smaller chunks or using a server

#### Empty gene files / Missing genes
- Check if genes are being filtered (negative controls, blanks, etc.)
- Review filtered genes in console output
- Use original gene names (before sanitization) when selecting in viewer

## Advanced Usage

### Processing Multiple Datasets
```bash
for file in *.h5ad; do
    output_dir="${file%.h5ad}_processed"
    python scripts/process_spatial_data.py "$file" "$output_dir/"
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
# python scripts/process_spatial_data.py data.h5ad output/ --chunk-size {chunk_size}
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
