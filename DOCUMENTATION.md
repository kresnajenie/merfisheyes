# MERFISH Eyes Documentation

Complete guide to dataset formats, requirements, and upload workflows for MERFISH Eyes - a web-based 3D visualization platform for spatial transcriptomics data.

---

## Overview

MERFISH Eyes is a web-based 3D visualization platform designed for spatial transcriptomics data. The platform supports two main data types:

- **Single Cell Data** - Cell-level aggregated data with expression matrices (.h5ad, MERSCOPE, Xenium)
- **Single Molecule Data** - Individual molecule coordinates with gene labels (.parquet, .csv)

The platform provides interactive 3D visualization using Three.js with comprehensive gene expression analysis capabilities.

---

## Single Cell Data

Single cell data represents cell-level aggregated information, including spatial coordinates, gene expression matrices, and cluster annotations.

### Supported Formats

MERFISH Eyes supports three main single cell data formats:

#### H5AD Format

**Description**: AnnData format (.h5ad) is a standard file format for annotated data matrices in Python.

**Requirements**:
- Must contain `obsm['X_spatial']` for spatial coordinates
- Optionally supports `obsm['X_umap']` for UMAP embeddings
- Cell type annotations can be stored in `obs` columns
- Gene expression matrix in `X` (sparse or dense)

**Example Structure**:
```python
adata.obs          # Cell metadata (cell types, clusters)
adata.obsm         # Spatial coordinates, embeddings
adata.X            # Gene expression matrix
adata.var          # Gene names and metadata
```

**File Size**: Typically 100MB - 5GB depending on number of cells and genes

**Upload Method**: Direct file upload through the web interface

---

#### Xenium Format

**Description**: 10x Genomics Xenium platform output format.

**Requirements**:
- Must have `cells.csv` with centroid columns (`x_centroid`, `y_centroid`)
- `transcripts.csv` for single molecule data (optional)
- Auto-detects cell type columns from metadata
- Folder structure with multiple CSV files

**Expected Files**:
```
xenium_folder/
├── cells.csv              # Cell centroids and metadata
├── transcripts.csv        # Single molecule coordinates (optional)
└── cell_feature_matrix/   # Gene expression matrix (optional)
```

**Key Columns in cells.csv**:
- `cell_id` - Unique cell identifier
- `x_centroid`, `y_centroid` - Spatial coordinates
- Cluster/cell type columns (auto-detected)

**Upload Method**: Upload entire folder as a ZIP file or select folder in file picker

---

#### MERSCOPE Format

**Description**: Vizgen MERSCOPE platform output format.

**Requirements**:
- Must have `cell_metadata.csv` with coordinate columns
- Supports `cell_by_gene.csv` for expression matrix
- Auto-detects spatial coordinate columns

**Expected Files**:
```
merscope_folder/
├── cell_metadata.csv      # Cell coordinates and metadata
├── cell_by_gene.csv       # Gene expression matrix
└── detected_transcripts.csv  # Single molecule data (optional)
```

**Key Columns in cell_metadata.csv**:
- `EntityID` or `cell_id` - Cell identifier
- `center_x`, `center_y` - Spatial coordinates
- Cell type and cluster annotations

**Upload Method**: Upload entire folder as a ZIP file or select folder in file picker

---

#### Pre-chunked Data (Python Preprocessing)

**Description**: For very large datasets (>500K cells) that exceed browser memory limits, use Python preprocessing to create chunked format.

**When to Use**:
- Datasets with >500,000 cells
- Browser runs out of memory during processing
- Need reproducible preprocessing pipeline
- Want to process data offline

**How to Preprocess**:

1. Install dependencies:
```bash
pip install anndata numpy pandas
```

2. Run preprocessing script:
```bash
# Auto-detect format (H5AD, Xenium, or MERSCOPE)
python scripts/process_spatial_data.py path/to/data output_folder

# Or process H5AD specifically
python scripts/process_h5ad.py path/to/data.h5ad output_folder
```

3. Upload the generated folder through the web interface

**Output Structure**:
```
output_folder/
├── manifest.json              # Dataset metadata
├── coords/
│   ├── spatial.bin.gz        # Normalized spatial coordinates
│   ├── umap.bin.gz           # UMAP embedding (max 3D)
│   └── pca.bin.gz            # PCA embedding (max 3D, truncated from 50)
├── expr/
│   ├── index.json            # Expression matrix index
│   └── chunk_00000.bin.gz    # Gene expression chunks
├── obs/
│   ├── metadata.json         # Cluster metadata
│   └── leiden.json.gz        # Cluster values
└── palettes/
    └── leiden.json           # Color palettes for categorical clusters
```

**Benefits**:
- ✅ No browser memory pressure (processes in Python with unlimited memory)
- ✅ Faster upload initialization (no client-side chunking)
- ✅ Handles very large datasets (>500K cells)
- ✅ Reproducible preprocessing outside browser environment

**See Also**: [Python Scripts](#python-scripts) section for detailed documentation

---

## Single Molecule Data

Single molecule data consists of individual RNA molecule coordinates with gene labels, providing subcellular resolution.

### Parquet/CSV Format

**Description**: Individual molecule coordinates stored in columnar format.

**Supported File Types**:
- `.parquet` - Columnar binary format (recommended for large datasets)
- `.csv` - Comma-separated values (simpler but less efficient)

**Requirements**:
- Must have columns for gene names and x/y coordinates
- z coordinates optional (for 3D datasets)
- Column names configurable via mappings

**Example Data Structure**:

| feature_name | x_location | y_location | z_location |
|--------------|------------|------------|------------|
| GAPDH        | 1234.5     | 5678.9     | 12.3       |
| ACTB         | 2345.6     | 6789.0     | 15.7       |
| MYC          | 3456.7     | 7890.1     | 18.2       |

**File Size Considerations**:
- Small datasets (<1M molecules): <50MB
- Medium datasets (1-10M molecules): 50MB - 500MB
- Large datasets (>10M molecules): >500MB
- Warning shown for files >2GB (potential memory issues)

**Upload Method**: Direct file upload through single molecule viewer

---

### Column Mappings

Different platforms use different column names. MERFISH Eyes supports configurable column mappings.

**Default Mappings**:

**Xenium**:
```javascript
{
  gene: 'feature_name',
  x: 'x_location',
  y: 'y_location',
  z: 'z_location'
}
```

**MERSCOPE**:
```javascript
{
  gene: 'gene',
  x: 'global_x',
  y: 'global_y',
  z: 'global_z'
}
```

**Custom**:
You can define your own column names in the upload settings.

**Configuration File**: `lib/config/moleculeColumnMappings.ts`

---

### Python Preprocessing for Single Molecule

For very large single molecule datasets, you can preprocess using Python to create optimized binary files.

**When to Use**:
- Datasets with >20M molecules
- Browser processing takes too long
- Need to filter/preprocess molecules before upload

**Future Feature**: Python preprocessing scripts for single molecule data are planned but not yet implemented. Currently, all single molecule processing happens in the browser using web workers.

---

## Column Type Detection Logic

MERFISH Eyes automatically detects whether metadata columns are categorical or numerical.

### Detection Rules

**Categorical Columns**:
- ≤100 unique values
- Displayed with discrete color palette
- Supports checkbox filtering in UI
- Examples: cell types, clusters, tissue regions

**Numerical Columns**:
- \>100 unique values
- Displayed with coolwarm gradient (blue → white → red)
- Supports range filtering with number scrubbers
- Examples: gene expression, quality metrics, confidence scores

### Visualization Behavior

**Categorical**:
- Each category gets a unique color from the palette
- Click cells to toggle category selection
- Hover to see category name and color

**Numerical**:
- Gradient coloring based on value range
- Interactive scalebar with min/max controls
- Auto-scaling to 95th percentile
- Manual override via drag controls

### Implementation

The type detection happens in multiple places:

- `H5adAdapter.isCategoricalData()` - For H5AD files
- `standardized-dataset.worker.ts` - For Xenium/MERSCOPE
- `ChunkedDataAdapter` - Reads type from S3 metadata

**Best Column Selection**: The `selectBestClusterColumn()` utility automatically picks the best cluster column with priority:
1. "leiden" column
2. Any column containing "celltype"
3. First categorical column
4. First available column

---

## Upload Workflow & S3 Storage

### Single Cell Upload Flow

**Standard Flow (H5AD, Xenium, MERSCOPE)**:

1. **Client-side Processing**:
   - Parse file in web worker (keeps UI responsive)
   - Calculate dataset fingerprint (SHA-256 hash)
   - Check for duplicates via API
   - Process data using `GeneChunkProcessor` to create chunks

2. **Initiate Upload**:
   - POST to `/api/datasets/initiate` with metadata and file list
   - Creates `Dataset`, `UploadSession`, and `UploadFile` records in PostgreSQL
   - Batch creates file records (`createMany`) for performance
   - Generates presigned S3 URLs in batches of 50
   - Returns presigned URLs for direct upload

3. **Upload Files**:
   - Client uploads directly to S3 using presigned URLs
   - No data passes through the server (efficient)
   - Progress tracking for each file

4. **Complete Upload**:
   - POST to `/api/datasets/[datasetId]/complete` to finalize
   - Updates dataset status to COMPLETE

5. **Email Notification**:
   - POST to `/api/send-email` with dataset details
   - Subject: `"{DatasetName} - Dataset Ready - MERFISHeyes"`
   - Includes: total cells, genes, platform, clusters, shareable link

6. **View Dataset**:
   - Navigate to `/viewer/[datasetId]`
   - Loads data via `StandardizedDataset.fromS3()`

**Pre-chunked Flow (Python Preprocessed)**:

1. **Preprocessing**:
   - Run Python script locally to create chunked folder
   ```bash
   python scripts/process_spatial_data.py data.h5ad output_folder
   ```

2. **Client-side**:
   - Select chunked folder on homepage
   - Reads `manifest.json` only (no browser processing)
   - Dataset marked with `isPreChunked = true`
   - Fingerprint generated from manifest hash

3. **Upload**:
   - Chunk size settings hidden (already processed)
   - Shows "Pre-chunked and ready for upload"
   - Skips `GeneChunkProcessor` entirely
   - Direct S3 upload of pre-chunked files

4. **Complete & View**:
   - Same as standard flow

---

### Single Molecule Upload Flow

1. **Client-side Processing**:
   - `SingleMoleculeProcessor.processAndUpload()` runs in web worker
   - Parses parquet/CSV file with progress callbacks
   - Generates fingerprint (gene names + per-gene molecule counts)
   - Checks for duplicates via `/api/single-molecule/check-duplicate/[fingerprint]`

2. **File Preparation**:
   - Creates `manifest.json.gz` with metadata (genes, dimensions, scaling factor, statistics)
   - Creates `genes/{geneName}.bin.gz` for each gene (gzipped Float32Array of normalized coordinates)

3. **Initiate Upload**:
   - POST to `/api/single-molecule/initiate` with fingerprint, metadata, manifest, and file list
   - Creates `Dataset` (datasetType="single_molecule"), `UploadSession`, and `UploadFile` records
   - Stores manifest JSON in database `manifestJson` field
   - Returns presigned S3 URLs

4. **Upload Files**:
   - Client uploads directly to S3
   - S3 path: `datasets/{datasetId}/manifest.json.gz` and `datasets/{datasetId}/genes/{gene}.bin.gz`
   - Marks each file complete via `/api/single-molecule/[id]/files/[key]/complete`

5. **Complete Upload**:
   - POST to `/api/single-molecule/[id]/complete`
   - Sends email notification

6. **View Dataset**:
   - Navigate to `/sm-viewer/[datasetId]`
   - Lazy-loads data via `SingleMoleculeDataset.fromS3()`

---

### Database Schema

The database tracks upload state using Prisma ORM with PostgreSQL.

**Key Models**:

**Dataset**:
- `datasetType` - "single_cell" or "single_molecule"
- `fingerprint` - SHA-256 hash for deduplication
- `status` - UPLOADING → PROCESSING → COMPLETE/FAILED
- `manifestJson` - Stores manifest for single molecule datasets
- `numCells` - Cell count (or molecule count for single molecule)

**UploadSession**:
- Tracks multi-file upload progress
- Expiration timestamp for cleanup

**UploadFile**:
- Individual file upload status within a session
- S3 key and upload URL

All models use cascade deletion for referential integrity.

---

## Load from Custom S3

MERFISH Eyes supports loading datasets from custom S3 buckets without uploading through the interface.

### Single Cell S3 Loading

**Requirements**:
- Pre-chunked format (use Python preprocessing scripts)
- Files uploaded to your own S3 bucket
- Proper CORS configuration on S3 bucket

**Loading Process**:

1. Navigate to homepage
2. Click "Load from S3" button
3. Enter S3 URL to `manifest.json` file
4. Platform fetches manifest and chunk metadata
5. Loads chunks on-demand as needed

**S3 URL Format**:
```
https://your-bucket.s3.amazonaws.com/path/to/dataset/manifest.json
```

**CORS Configuration**:
Your S3 bucket must allow cross-origin requests:
```json
{
  "AllowedOrigins": ["https://merfisheyes.com"],
  "AllowedMethods": ["GET"],
  "AllowedHeaders": ["*"]
}
```

---

### Single Molecule S3 Loading

**Lazy Loading Architecture**:

Single molecule datasets use lazy loading to minimize initial load time and memory usage.

**Loading Process**:

1. **Fetch Metadata** (10%):
   - GET `/api/single-molecule/[id]` to get manifest URL

2. **Download Manifest** (30%):
   - Fetch and ungzip `manifest.json.gz` from S3
   - Contains gene list and dataset metadata

3. **Initialize Dataset** (60%):
   - Create dataset with gene list
   - Empty geneIndex (genes not loaded yet)

4. **Lazy Load Genes** (on-demand):
   - When user selects a gene, fetch presigned URL
   - Download `genes/{geneName}.bin.gz` from S3
   - Ungzip → Float32Array → number[] conversion
   - Cache in memory for instant re-selection

**Gene File Binary Format** (.bin.gz):
- Gzipped Float32Array of normalized coordinates
- Flat array: `[x1, y1, z1, x2, y2, z2, ...]` for 3D
- Each coordinate is 4 bytes (Float32)
- Coordinates already normalized to [-1, 1] range
- Gene names sanitized in filenames (special chars → underscores)

**Benefits**:
- ✅ Fast initialization (~1 second for manifest only)
- ✅ Memory efficient (only loads selected genes)
- ✅ Instant re-selection (genes cached after first load)
- ✅ Scales to datasets with hundreds of genes

**API Endpoints**:
- `/api/single-molecule/[id]` - Returns dataset metadata and manifest URL
- `/api/single-molecule/[id]/gene/[geneName]` - Returns presigned URL for specific gene file

---

## Python Scripts

For very large datasets (>500K cells) that exceed browser memory limits, Python scripts can preprocess data into chunked format.

### Installation

```bash
pip install anndata numpy pandas
```

### Available Scripts

**1. process_spatial_data.py** - Auto-detect and process any format

```bash
python scripts/process_spatial_data.py path/to/data output_folder
```

Automatically detects:
- H5AD files (.h5ad extension)
- Xenium folders (contains cells.csv)
- MERSCOPE folders (contains cell_metadata.csv)

**2. process_h5ad.py** - Process H5AD files specifically

```bash
python scripts/process_h5ad.py path/to/data.h5ad output_folder
```

### Output Structure

```
output_folder/
├── manifest.json              # Dataset metadata
├── coords/
│   ├── spatial.bin.gz        # Normalized spatial coordinates
│   ├── umap.bin.gz           # UMAP embedding (max 3D)
│   └── pca.bin.gz            # PCA embedding (max 3D, truncated from 50)
├── expr/
│   ├── index.json            # Expression matrix index
│   └── chunk_00000.bin.gz    # Gene expression chunks
├── obs/
│   ├── metadata.json         # Cluster metadata
│   └── leiden.json.gz        # Cluster values
└── palettes/
    └── leiden.json           # Color palettes for categorical clusters
```

### Key Features

**Embedding Dimension Limit**:
- Automatically truncates embeddings (PCA, UMAP, etc.) to first 3 dimensions
- Example: PCA with 50 components → saves only first 3 (94% size reduction)

**Coordinate Normalization**:
- Scales to [-1, 1] range with saved scaling factor
- Ensures consistent visualization

**Sparse Expression Matrix**:
- Stores only non-zero values with indices
- Significant size reduction for sparse data

**Binary Compression**:
- Uses gzip compression for all binary files
- Typical 3-5x size reduction

**Automatic Column Detection**:
- Detects categorical vs numerical metadata columns
- Creates appropriate visualizations

**Color Palette Generation**:
- Creates consistent color palettes matching browser visualization
- Saves to `palettes/` directory

### Documentation

See detailed documentation:
- `scripts/README.md` - Script usage and options
- `scripts/UPLOAD_GUIDE.md` - Step-by-step upload guide

---

## Troubleshooting

### Common Issues

#### "Out of Memory" Error During Upload

**Symptoms**:
- Browser tab crashes during file processing
- "Aw, Snap!" error in Chrome
- Processing stuck at certain percentage

**Solutions**:
1. Use Python preprocessing for large datasets (>500K cells)
2. Close other browser tabs to free memory
3. Try uploading in smaller batches
4. Use a computer with more RAM

**Prevention**:
- Use Python preprocessing for datasets >500K cells
- Pre-chunk data before uploading

---

#### File Format Not Recognized

**Symptoms**:
- "Invalid file format" error
- Upload fails immediately
- No processing progress

**Solutions**:
1. Check file format requirements above
2. Ensure required columns are present
3. Verify file is not corrupted (try opening in Python/Excel)
4. Check file extension matches content (.h5ad for H5AD, .parquet for Parquet)

**Common Causes**:
- Missing required columns (e.g., `X_spatial` in H5AD)
- Wrong file extension
- Corrupted file download

---

#### Slow Upload Speed

**Symptoms**:
- Upload takes very long time
- Progress bar stuck
- Network errors

**Solutions**:
1. Check internet connection speed
2. Use pre-chunked format to reduce upload size
3. Try uploading during off-peak hours
4. Contact support if persistent

**Optimization**:
- Pre-chunking can reduce upload size by 50-70%
- Parallel uploads happen automatically (50 files at a time)

---

#### Visualization Not Loading

**Symptoms**:
- Blank screen after upload complete
- Loading spinner never stops
- Console errors

**Solutions**:
1. Refresh the page
2. Check browser console for errors (F12)
3. Try different browser (Chrome recommended)
4. Clear browser cache
5. Check dataset status in database

**Common Causes**:
- Incomplete upload (check S3 for missing files)
- Browser compatibility issues
- Corrupted manifest file

---

#### Colors Don't Match Expected Cell Types

**Symptoms**:
- Cell type colors seem random
- Colors change between sessions

**Explanation**:
- Colors are assigned based on category order in the data
- First category gets first palette color, etc.
- Consistent within a dataset, but may differ between datasets

**Solution**:
- Colors are cosmetic and don't affect analysis
- Use the legend to identify which color represents which cell type

---

#### Gene Expression Values Seem Wrong

**Symptoms**:
- All values near zero
- Unexpected expression patterns

**Checks**:
1. Verify data is log-normalized (not raw counts)
2. Check if expression matrix is in correct format
3. Confirm gene names match between matrix and metadata
4. Look at scale bar range (may need manual adjustment)

**Scale Bar**:
- Auto-scales to 95th percentile by default
- Can manually adjust min/max with number scrubbers
- Try adjusting range if visualization looks washed out

---

### Getting Help

**Documentation**:
- Read this documentation thoroughly
- Check `CLAUDE.md` for technical details
- Review Python script documentation in `scripts/README.md`

**GitHub Issues**:
- Search existing issues: https://github.com/kresnajenie/merfisheyes/issues
- Create new issue with:
  - Dataset format and size
  - Browser and version
  - Error messages (with screenshots)
  - Steps to reproduce

**Discord Community**:
- Join: https://discord.gg/BRp6C2EVHU
- Ask questions in #support channel
- Share datasets and visualizations

**Email Support**:
- For sensitive data issues
- Contact: support@merfisheyes.com (if available)

---

## Additional Resources

### Links

- **GitHub Repository**: https://github.com/kresnajenie/merfisheyes
- **Discord Community**: https://discord.gg/BRp6C2EVHU
- **Twitter**: https://x.com/merfisheyes
- **Patreon**: https://www.patreon.com/cw/MERFISHEYES

### Example Datasets

Coming soon - we'll provide example datasets for testing:
- Small H5AD file (<100MB)
- Xenium sample data
- MERSCOPE sample data
- Single molecule parquet file

### Citation

If you use MERFISH Eyes in your research, please cite:

```
MERFISH Eyes - Web-based 3D Visualization for Spatial Transcriptomics
https://merfisheyes.com
```

(Formal citation to be added upon publication)

---

## Appendix

### Supported File Formats Summary

| Format | Type | Extension | Max Size | Processing |
|--------|------|-----------|----------|------------|
| H5AD | Single Cell | .h5ad | 5GB | Browser or Python |
| Xenium | Single Cell | folder | 10GB | Browser or Python |
| MERSCOPE | Single Cell | folder | 10GB | Browser or Python |
| Pre-chunked | Single Cell | folder | Unlimited | Python only |
| Parquet | Single Molecule | .parquet | 2GB | Browser |
| CSV | Single Molecule | .csv | 500MB | Browser |

### Browser Requirements

**Recommended**:
- Google Chrome 100+
- 16GB RAM minimum
- Modern GPU for WebGL
- Fast internet connection (10+ Mbps)

**Minimum**:
- Chrome/Firefox/Safari (latest version)
- 8GB RAM
- WebGL 2.0 support
- 5+ Mbps internet

### API Rate Limits

- Upload initiation: 10 requests per minute
- S3 presigned URL generation: 100 URLs per request
- Gene data fetching: No limit (cached after first load)

### Storage Limits

- Free tier: 10GB total storage
- Pro tier: 100GB total storage (coming soon)
- Enterprise: Custom limits

---

*Last updated: February 2026*
*Version: 1.0.0*
