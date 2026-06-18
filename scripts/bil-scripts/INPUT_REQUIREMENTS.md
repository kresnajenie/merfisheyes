# Input Requirements

This document describes the exact file formats, column names, and structure required for each pipeline.

## Table of Contents

- [Case 1 & 2: MERSCOPE CSV Files](#case-1--2-merscope-csv-files)
- [Case 3: H5AD Files](#case-3-h5ad-files)
- [Column Reference Tables](#column-reference-tables)

---

## Case 1 & 2: MERSCOPE CSV Files

### cell_metadata.csv

Contains spatial coordinates and metadata for each cell.

**Required columns:**

| Column | Description | Notes |
|--------|-------------|-------|
| Cell ID | Unique identifier per cell | Auto-detected: searches for `EntityID`, `id`, `cell_id` in order. Falls back to first column. |
| `center_x` | X coordinate of cell centroid | Numeric (microns) |
| `center_y` | Y coordinate of cell centroid | Numeric (microns) |

**Optional columns:**

| Column | Description |
|--------|-------------|
| `center_z` | Z coordinate (for 3D datasets) |
| Any other columns | Preserved as metadata (cluster annotations, QC metrics, etc.) |

**Example:**
```csv
EntityID,center_x,center_y,volume,class_name
cell_0001,1234.56,5678.90,150.3,Neuron
cell_0002,1235.12,5680.45,142.1,Astrocyte
cell_0003,1236.78,5682.10,165.7,Neuron
```

**Format:** UTF-8 CSV with header row. Comma-delimited.

---

### cell_by_gene.csv

Expression matrix — one row per cell, one column per gene.

**Required columns:**

| Column | Description | Notes |
|--------|-------------|-------|
| First column | Cell ID | Must match cell IDs in `cell_metadata.csv`. Column name is auto-detected (commonly `cell`). |
| Remaining columns | Gene expression values | Column names = gene names. Values = raw counts (integer or float). |

**Example:**
```csv
cell,Slc17a7,Gad1,Th,Drd1,Camk2a
cell_0001,5,0,0,2,12
cell_0002,0,3,0,0,1
cell_0003,8,0,1,0,15
```

**Important:**
- Cell IDs must match between `cell_metadata.csv` and `cell_by_gene.csv`
- Rows in different order are OK — they are aligned by cell ID during processing
- Gene names become the searchable gene list in the viewer

---

### detected_transcripts.csv

Individual molecule coordinates with gene labels. Used by the single molecule pipeline.

**Required columns (MERSCOPE default):**

| Column | Description | Notes |
|--------|-------------|-------|
| `gene` | Gene name | String |
| `global_x` | X coordinate | Numeric (microns) |
| `global_y` | Y coordinate | Numeric (microns) |

**Optional columns:**

| Column | Description | Notes |
|--------|-------------|-------|
| `global_z` | Z coordinate | Numeric. If present, dataset is treated as 3D. |
| `cell_id` | Cell assignment | Integer. `-1` = unassigned molecule. Unassigned molecules are written to separate files. |

**Example:**
```csv
gene,global_x,global_y,global_z,cell_id
Slc17a7,1234.56,5678.90,2.5,42
Gad1,1235.12,5680.45,1.8,-1
Th,1236.78,5682.10,3.1,42
Drd1,2345.67,6789.01,2.0,105
```

**Column name matching:** Fuzzy — the file name is matched by keywords ("detected" + "transcript"). Column names must match exactly (or be overridden with `--gene-col`, `--x-col`, etc.).

**Other dataset types:** Use command-line overrides for non-MERSCOPE column names:

| Dataset Type | Gene | X | Y | Z |
|-------------|------|---|---|---|
| MERSCOPE (default) | `gene` | `global_x` | `global_y` | `global_z` |
| Xenium | `feature_name` | `x_location` | `y_location` | `z_location` |
| Custom | Use `--gene-col`, `--x-col`, `--y-col`, `--z-col` |

---

### segmented_spot_table.csv (H5AD pipeline)

Used by `launch_h5ad_pipeline.sh` for single molecule processing alongside H5AD files.

**Required columns:**

| Column | Description |
|--------|-------------|
| `gene_names` | Gene name |
| `x` | X coordinate |
| `y` | Y coordinate |
| `z` | Z coordinate |
| `cell_ids` | Cell assignment (integer, `-1` = unassigned) |

**Example:**
```csv
gene_names,x,y,z,cell_ids
Slc17a7,1234.56,5678.90,2.5,42
Gad1,1235.12,5680.45,1.8,-1
Th,1236.78,5682.10,3.1,42
```

**Note:** These column names differ from MERSCOPE's `detected_transcripts.csv`. The H5AD pipeline passes custom column mappings automatically (`--gene-col gene_names --x-col x --y-col y --z-col z --cell-id-col cell_ids`).

---

## Case 3: H5AD Files

### cell_by_gene.h5ad

AnnData H5AD file containing single cell expression data.

**Required structure:**

| Key | Description | Notes |
|-----|-------------|-------|
| `X` or `obsm['X_raw']` | Expression matrix | Raw counts preferred. If `X` appears normalized, falls back to `obsm['X_raw']`. |
| `obsm['X_spatial']` or `obsm['spatial']` | Spatial coordinates | Shape: `(n_cells, 2)` or `(n_cells, 3)`. First choice: `X_spatial`, fallback: `spatial`. |
| `var` index | Gene names | Searched in order: `_index`, `gene`, `genes` |

**Optional structure:**

| Key | Description | Notes |
|-----|-------------|-------|
| `obsm['X_umap']` | UMAP embedding | Any `X_*` key (except `X_spatial`) is loaded as an embedding. Truncated to first 3 dimensions. |
| `obsm['X_pca']` | PCA embedding | Same as above. |
| `obs` columns | Cell metadata | All non-coordinate columns loaded as cluster/metadata. Auto-detected as categorical or numerical. |
| `uns['{column}_colors']` | Color palette | Hex color array matching unique values in the corresponding obs column. |

**Spatial coordinate fallback:** If `obsm['X_spatial']` is not found, the adapter searches `obs` columns in this order:

| Priority | X column | Y column | Z column |
|----------|----------|----------|----------|
| 1st | `center_x` | `center_y` | `center_z` |
| 2nd | `centerX` | `centerY` | `centerZ` |
| 3rd | `x` | `y` | `z` |
| 4th | `centroid_x` | `centroid_y` | `centroid_z` |

**Coordinate validation:**
- X and Y: at least 90% of values must be valid numbers
- Z: at least 90% valid AND at least 90% of valid values must be non-zero (otherwise treated as 2D)
- Rows with invalid X or Y coordinates are filtered out

**Expression data:**
- Raw counts expected (integers or near-integers)
- If the matrix appears normalized (no zeros or all floats), `map_my_cell` will look for `obsm['X_raw']` as fallback

**Example (Python):**
```python
import anndata as ad
import numpy as np
import pandas as pd

# Create minimal valid H5AD
n_cells = 1000
n_genes = 500

adata = ad.AnnData(
    X=np.random.poisson(2, (n_cells, n_genes)),  # raw counts
    obs=pd.DataFrame({
        'cell_type': np.random.choice(['Neuron', 'Astrocyte', 'Microglia'], n_cells),
    }),
    var=pd.DataFrame(index=[f'Gene_{i}' for i in range(n_genes)]),
)

# Add spatial coordinates
adata.obsm['X_spatial'] = np.random.uniform(0, 10000, (n_cells, 2))

# Save
adata.write_h5ad('cell_by_gene.h5ad')
```

---

## Column Reference Tables

### Single Cell Pipeline (combine_slices → process_spatial)

| File | Column | Required | Auto-detected | Fallback |
|------|--------|----------|---------------|----------|
| cell_metadata.csv | Cell ID | Yes | `EntityID` → `id` → `cell_id` → first column | |
| cell_metadata.csv | `center_x` | Yes | Exact match | |
| cell_metadata.csv | `center_y` | Yes | Exact match | |
| cell_metadata.csv | `center_z` | No | Exact match | Treated as 2D |
| cell_by_gene.csv | Cell ID | Yes | First column | |
| cell_by_gene.csv | Gene columns | Yes | All remaining columns | |

### Single Molecule Pipeline (process_single_molecule)

| Dataset Type | Gene Column | X Column | Y Column | Z Column | Cell ID Column |
|-------------|-------------|----------|----------|----------|----------------|
| MERSCOPE (default) | `gene` | `global_x` | `global_y` | `global_z` | `cell_id` |
| Xenium | `feature_name` | `x_location` | `y_location` | `z_location` | — |
| H5AD spot table | `gene_names` | `x` | `y` | `z` | `cell_ids` |
| Custom | `--gene-col` | `--x-col` | `--y-col` | `--z-col` | `--cell-id-col` |

### H5AD Structure (process_spatial_data / H5adAdapter)

| Key | Required | Description |
|-----|----------|-------------|
| `X` | Yes | Expression matrix (n_cells x n_genes), raw counts |
| `obsm['X_spatial']` | Yes* | Spatial coordinates (n_cells x 2 or 3) |
| `var._index` | Yes | Gene names |
| `obs.*` | No | Metadata columns (auto-loaded as clusters) |
| `obsm['X_umap']` | No | UMAP embedding |
| `obsm['X_pca']` | No | PCA embedding |
| `uns['{col}_colors']` | No | Color palettes for categorical obs columns |

*Falls back to `obs` coordinate columns if `obsm['X_spatial']` is missing.

### File Discovery (Fuzzy Matching)

The pipeline uses fuzzy keyword matching to find files in directories:

| Target File | Keywords (all must match) | Excludes | Extension |
|-------------|--------------------------|----------|-----------|
| cell_metadata | "metadata" | "gene", "detected", "transcript" | .csv |
| cell_by_gene | "cell" + "gene" | "metadata", "detected", "transcript" | .csv |
| detected_transcripts | "detected" + "transcript" | "metadata", "cell_by_gene" | .csv |

**Examples of matched filenames:**
- `cell_metadata.csv` -> cell_metadata
- `cellpose_metadata.csv` -> cell_metadata (contains "metadata")
- `cell_by_gene_exons.csv` -> cell_by_gene (contains "cell" + "gene")
- `detected_transcripts.csv` -> detected_transcripts
- `detected_transcripts_v2.csv` -> detected_transcripts
