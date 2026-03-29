# Data Input Requirements

This document defines the exact input requirements for each supported data format, derived from the adapter implementations. All column name matching is **case-insensitive** unless noted otherwise.

---

## Single Cell Formats

### 1. H5AD (.h5ad)

**Source**: `lib/adapters/H5adAdapter.ts`

#### Spatial Coordinates (REQUIRED)

Checked in this order — first match wins:

1. `obsm/X_spatial` — shape `[n_cells, 2]` or `[n_cells, 3]`
2. `obsm/spatial` — same shape
3. `obs` columns — tries these naming patterns in order:
   - `center_x`, `center_y`, `center_z`
   - `centerX`, `centerY`, `centerZ`
   - `x`, `y`, `z`
   - `centroid_x`, `centroid_y`, `centroid_z`

**Validation**: X and Y each require >=90% valid finite floats. Z (if present) requires >=90% valid AND >=90% non-zero to be treated as 3D; otherwise falls back to 2D.

**Error**: Throws if no coordinates found from any source.

#### Gene Names (OPTIONAL)

Checked in this order:

1. `var/_index`
2. `var/gene`
3. `var/genes`

Returns empty array if all fail (non-fatal).

#### Expression Matrix (REQUIRED for gene visualization)

- Reads `X` dataset from HDF5 root
- Shape: `[n_cells, n_genes]`
- Dense array loaded fully into memory

> **TODO: H5AD Sparse Matrix Support**
> Currently only dense `X` matrices are supported. Many H5AD files store expression data as sparse matrices (CSR/CSC format via `scipy.sparse`). Support for reading sparse `X` matrices needs to be implemented.

> **TODO: H5AD Zarr Backend**
> Some AnnData files use Zarr as a backing store instead of HDF5. Support for reading `.h5ad` files backed by Zarr (or standalone `.zarr` directories) is not yet implemented.

#### Embeddings (OPTIONAL)

- Scans all keys in `obsm` starting with `X_` (except `X_spatial`)
- `X_umap` stored as `"umap"`, `X_pca` as `"pca"`, etc.

#### Cluster/Observation Columns (OPTIONAL)

- Reads all columns from `obs` group that don't start with `_`
- Each column classified as categorical or numerical:
  - **Always categorical**: column name contains `leiden` or `louvain`
  - **Always numerical**: column name matches `n_genes`, `total_counts`, `n_counts`, `pct_counts`, `log1p_*`, `n_cells`, `doublet_score`
  - **Auto-detected**: >50% float values → numerical; integers with >=80% uniqueness → numerical; otherwise categorical

#### Color Palette

- First tries `uns/{columnName}_colors` (from Scanpy/AnnData)
- Validates palette matches unique cluster count and has >1 unique color
- Falls back to built-in `DEFAULT_COLOR_PALETTE` (40 colors, cycles)

---

### 2. Xenium (folder)

**Source**: `lib/adapters/XeniumAdapter.ts`

#### Required Files

- `cells.csv` or `cells.csv.gz` — **REQUIRED**

#### Optional Files

- `features.tsv` (or `.tsv.gz`, `.csv`, `.csv.gz`) — gene names
- `cell_feature_matrix/features.tsv` (same variants) — gene names (alternate location)
- `cell_feature_matrix/matrix.mtx` (or `.mtx.gz`) — expression matrix (Matrix Market format)
- `cell_feature_matrix/barcodes.tsv` (or `.tsv.gz`) — needed with matrix.mtx
- `analysis/clustering/*.csv` or `*.tsv` — cluster labels

#### Cell ID Column (from cells.csv)

Checked in order:

1. `cell_id`
2. `id`
3. `barcode`
4. `cell`
5. `cellid`
6. `cell_id_original`
7. Falls back to first column

#### Spatial Coordinates (from cells.csv)

**X coordinate** — checked in order:

1. `cell_centroid_x`
2. `x_centroid`
3. `centroid_x`
4. `center_x`
5. `centerx`
6. `cx`
7. `x`
8. `x_um`
9. `x_px`
10. `x_position`
11. `xpos`
12. Heuristic: any numeric column matching `/(^|_)x($|_|[a-z])|centroid.*x|x.*centroid|cx/i`

**Y coordinate** — same pattern with `y` substituted.

**Validation**: >=90% of rows must have valid numeric X,Y. Rows with invalid coords are dropped. Error if <90% valid.

**Dimensions**: Always 2D.

#### Cluster Column (from cells.csv)

Checked in order (40 candidate names):

`cell_type`, `celltype`, `cell type`, `cell-type`, `celltypes`, `cell types`, `predicted_celltype`, `predicted celltype`, `predicted cell type`, `cluster`, `clusters`, `leiden`, `louvain`, `cell.types`, `cell_types`, `annotation`, `annotations`, `class`, `class_label`, `class label`, `subclass`, `subclass_label`, `subclass label`, `label`, `labels`, `cell_annotation`, `cell annotation`, `cell_annotations`, `cell annotations`, `cell_class`, `cell class`, `cell_subclass`, `cell subclass`, `cluster_label`, `cluster label`, `cluster_id`, `cluster id`, `group`, `group_label`, `group label`

**Guard**: If selected column has >max(50, 10% of cells) unique values, it's rejected (likely an ID column). Falls back to single "All" group.

#### Cluster Labels from Analysis Folder (OPTIONAL)

- Searches `analysis/clustering/*.csv` or `*.tsv`
- Prefers `analysis/clustering/gene_expression_graphclust/`
- Picks file matching the most cells with reasonable unique label count

#### Gene Names (from features file)

- Format: `gene_name\tfeature_type` (tab or comma delimited)
- Type column is optional

---

### 3. MERSCOPE (folder)

**Source**: `lib/adapters/MerscopeAdapter.ts`

#### Required Files

- `cell_metadata.csv` (or `.csv.gz`) — **REQUIRED**

#### Optional Files

- `cell_categories.csv` (or `.csv.gz`) — cluster labels
- `cell_numeric_categories.csv` (or `.csv.gz`) — UMAP embeddings
- `cell_by_gene.csv` (or `.csv.gz`) — gene expression

#### Cell ID Column

Checked in order:

1. `EntityID`
2. `cell`
3. `cell_id`
4. `id`

Falls back to first column with warning.

#### Spatial Coordinates (from cell_metadata.csv)

**X coordinate** — checked in order:

1. `center_x`
2. `centroid_x`
3. `centerX`
4. `centerx`
5. `x`
6. `x_centroid`
7. `x_px`
8. `x_position`

**Y coordinate** — same pattern with `y` substituted.

**Validation**: >=90% valid. Error if <90%.

**Dimensions**: Always 2D.

#### Cluster Column

**Stage 1** — from `cell_categories.csv`:

`leiden`, `cluster`, `clusters`, `cell_type`, `celltype`, `annotation`, `annotations`, `class`, `subclass`, `label`, `labels`

**Stage 2** — from obs keys (24+ candidates, same as Xenium list above)

**Stage 3** — heuristic detection:
- Skips columns with prefixes: `x`, `y`, `z`, `row`, `col`, `column`, `coord`, `n_`, `sum`, `total`, `count`, `umi`, `reads`, `intensity`, `area`, `volume`
- Skips columns containing: `gene`, `umi`
- Requires: 2+ unique values, <=200 unique, <=20% of rows, <25% numeric-looking
- Selects column with lowest score (unique count + numeric ratio)

Falls back to single "All" group if nothing found.

#### UMAP Embeddings (from cell_numeric_categories.csv)

- ID column: `EntityID`, `cell`, `cell_id`, `id`
- UMAP X: `umap_X`, `umap_x`, `UMAP_X`, `umap1`, `x_umap`
- UMAP Y: `umap_Y`, `umap_y`, `UMAP_Y`, `umap2`, `y_umap`

Only created if all three columns found. 2D only.

#### Gene Expression (from cell_by_gene.csv)

**Long format** — detected if all of these columns exist:
- `cell` (cell ID)
- `gene` (gene name)
- One of: `count`, `expression`, `value` (expression value, priority in that order)

**Wide format** — default fallback:
- First column is cell ID (detected from: `cell`, `EntityID`, `cell_id`, `id`)
- All other columns (not starting with `_`) are gene names
- One row per cell, one column per gene

---

### 4. Pre-Chunked (Python-preprocessed folder)

**Source**: `lib/adapters/ChunkedDataAdapter.ts`, `scripts/process_h5ad.py`

#### Required Folder Structure

```
output_folder/
├── manifest.json                    # REQUIRED
├── coords/
│   └── spatial.bin.gz              # REQUIRED
├── expr/
│   ├── index.json                  # REQUIRED
│   └── chunk_00000.bin.gz          # REQUIRED (at least one)
├── obs/
│   ├── metadata.json               # OPTIONAL (needed for cluster columns)
│   └── {column_name}.json.gz       # OPTIONAL (one per cluster column)
└── palettes/
    └── {column_name}.json          # OPTIONAL (one per categorical column)
```

#### manifest.json

```json
{
  "version": "1.0",
  "statistics": {
    "total_cells": <int>,
    "total_genes": <int>,
    "spatial_dimensions": 2 | 3,
    "available_embeddings": ["umap", "pca", ...],
    "cluster_count": <int>
  },
  "processing": {
    "spatial_scaling_factor": <float>,
    "chunk_size": <int>,
    "num_chunks": <int>
  },
  "files": {
    "coordinates": ["spatial", "umap", "pca", ...],
    "expression_chunks": <int>,
    "observation_columns": ["leiden", "celltype", ...]
  }
}
```

#### coords/*.bin.gz (Binary, gzipped, little-endian)

| Offset | Type | Field |
|--------|------|-------|
| 0 | uint32 | num_points |
| 4 | uint32 | dimensions (2 or 3) |
| 8 | float32[num_points * dimensions] | coordinates (row-major: x1,y1,z1,x2,y2,z2,...) |

Coordinates are normalized to [-1, 1].

Embeddings truncated to first 3 dimensions (e.g. PCA 50 dims → 3).

#### expr/index.json

```json
{
  "total_genes": <int>,
  "num_chunks": <int>,
  "chunk_size": <int>,
  "genes": [
    { "name": "GENE_A", "chunk_id": 0, "position_in_chunk": 0 },
    { "name": "GENE_B", "chunk_id": 0, "position_in_chunk": 1 }
  ]
}
```

#### expr/chunk_XXXXX.bin.gz (Binary, gzipped, little-endian)

| Offset | Type | Field |
|--------|------|-------|
| 0 | uint32 | version (always 1) |
| 4 | uint32 | num_genes |
| 8 | uint32 | chunk_id |
| 12 | uint32 | total_cells |
| 16 | Gene table: 24 bytes per gene (gene_index, data_offset, data_size, uncompressed_size, num_non_zero, reserved) |
| varies | Per-gene sparse data: num_cells(u32), num_non_zero(u32), indices(u32[]), values(f32[]) |

#### obs/metadata.json

```json
{
  "column_name": { "type": "categorical" | "numerical" }
}
```

#### obs/{column_name}.json.gz

JSON array, one value per cell: `["Neuron", "Astrocyte", "Neuron", ...]`

#### palettes/{column_name}.json

Only for categorical columns: `{ "Neuron": "#1f77b4", "Astrocyte": "#ff7f0e" }`

#### Chunk Size Auto-Determination (Python script)

| Gene Count | Chunk Size |
|------------|------------|
| < 100 | 50 |
| < 500 | 100 |
| < 2000 | 200 |
| < 10000 | 500 |
| >= 10000 | 1000 |

---

## Single Molecule Formats

### 5. Parquet (.parquet)

**Source**: `lib/SingleMoleculeDataset.ts`, `lib/config/moleculeColumnMappings.ts`

Parsed via `hyparquet` (pure JavaScript, column-oriented reading).

#### Column Mappings by Dataset Type

| Dataset Type | Gene Column | X Column | Y Column | Z Column (optional) |
|-------------|-------------|----------|----------|---------------------|
| **xenium** | `feature_name` | `x_location` | `y_location` | `z_location` |
| **merscope** | `gene` | `global_x` | `global_y` | `global_z` |
| **custom** | `feature_name` | `x_location` | `y_location` | `z_location` |

Z column is optional — if present and non-empty, dataset is 3D; otherwise 2D.

#### Gene Filtering

Automatically filters control probes and unassigned genes via `shouldFilterGene()`:
- Negative controls
- "unassigned" genes
- Deprecated markers
- Codewords
- Blanks

---

### 6. CSV (.csv)

**Source**: `lib/SingleMoleculeDataset.ts`

Parsed via `papaparse` with `dynamicTyping: true` and `header: true`.

Same column mappings as Parquet (see table above). Same gene filtering.

---

### 7. S3 Single Molecule (uploaded dataset)

**Source**: `lib/SingleMoleculeDataset.ts`, `lib/utils/SingleMoleculeProcessor.ts`

#### S3 Structure

```
datasets/{datasetId}/
├── manifest.json.gz        # Dataset metadata (gzipped JSON)
└── genes/
    └── {gene_name}.bin.gz  # Per-gene coordinate file (gzipped Float32Array)
```

#### manifest.json.gz

```json
{
  "version": "1.0",
  "statistics": {
    "total_molecules": <int>,
    "unique_genes": <int>,
    "spatial_dimensions": 2 | 3
  },
  "genes": {
    "unique_gene_names": ["GENE1", "GENE2", ...]
  },
  "processing": {
    "compression": "gzip",
    "coordinate_format": "float32_flat_array",
    "coordinate_range": "normalized_[-1,1]",
    "scaling_factor": <float>
  }
}
```

#### genes/{gene_name}.bin.gz

- Gzipped Float32Array of normalized coordinates
- Flat array: `[x1, y1, z1, x2, y2, z2, ...]` (3D) or `[x1, y1, x2, y2, ...]` (2D)
- Each coordinate: 4 bytes (float32)
- Coordinates already normalized to [-1, 1]

**Gene name sanitization for filenames**: non-alphanumeric chars replaced with `_`, collapsed, stripped from edges. Example: `DAPI-488` → `DAPI_488`

#### Loading (lazy)

1. Fetch manifest only on initial load
2. Gene data loaded on demand when user selects a gene
3. Once loaded, cached in memory for instant re-selection

---

## Common Invariants

These hold for all formats after processing:

| Property | Requirement |
|----------|-------------|
| Spatial coordinates | Normalized to [-1, 1] range |
| Coordinate dimensions | 2 or 3 |
| Cluster uniqueValues | Sorted with `localeCompare({ numeric: true })` |
| Cluster valueIndices | Length equals point count |
| Cluster valueIndices type | Uint16Array if <65536 unique values, Uint32Array otherwise |
| Color palette keys | Match uniqueValues entries |
| Gene names | Non-empty strings |

---

## TODOs / Planned Changes

### Coordinate Normalization

Currently all spatial coordinates are normalized to [-1, 1] during loading (both in browser adapters and Python preprocessing). This adds processing overhead and loses the original coordinate scale.

**Planned**: Move normalization out of the data loading pipeline. Store raw coordinates and handle scaling at the visualization layer instead. This would:
- Preserve original coordinate values for analysis
- Allow toggling between raw and normalized views
- Reduce processing time during loading
- Simplify the data pipeline (one less transformation step)

### H5AD Sparse Matrix

See TODO in H5AD section above. Required for many real-world datasets that store expression as sparse.

### H5AD Zarr Backend

See TODO in H5AD section above. Required for large datasets that exceed HDF5 file size limits or use cloud-native storage.
