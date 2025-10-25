# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MERFISH Eyes is a web-based 3D visualization platform for spatial transcriptomics data. It supports both:
- **Single Cell Data** (.h5ad, MERSCOPE, Xenium) - Cell-level aggregated data with expression matrices
- **Single Molecule Data** (.parquet, .csv) - Individual molecule coordinates with gene labels

The platform provides interactive 3D visualization using Three.js with gene expression analysis capabilities.

## Development Commands

### Basic Commands
```bash
npm install              # Install dependencies
npm run dev             # Start dev server with Turbopack (http://localhost:3000)
npm run build           # Production build (4GB memory)
npm run build:low-memory # Production build for low-memory servers (2GB)
npm run build:strict    # Production build with ESLint enabled
npm start               # Start production server
npm run lint            # Run ESLint with auto-fix
```

### Database Commands
```bash
npx prisma generate     # Generate Prisma client after schema changes
npx prisma migrate dev  # Create and apply database migrations
npx prisma studio       # Open Prisma Studio GUI
npm run test-prisma     # Test Prisma connection
```

### Environment Setup
Required environment variables in `.env` or `.env.local`:
- `DATABASE_URL` - PostgreSQL connection string
- `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_REGION`, `AWS_S3_BUCKET` - S3 storage
- `NEXT_PUBLIC_BASE_URL` - Base URL for the application

## Architecture

### Data Processing Pipeline

The application supports two distinct data types with separate processing pipelines:

#### Single Cell Data Pipeline

Uses an **Adapter Pattern** to standardize different formats into a unified `StandardizedDataset`:

1. **Format-Specific Adapters** ([lib/adapters/](lib/adapters/)):
   - `H5adAdapter` - Parses AnnData .h5ad files using h5wasm
   - `XeniumAdapter` - Parses Xenium folder structure (cells.csv, transcripts.csv)
   - `MerscopeAdapter` - Parses MERSCOPE folder structure (cell_metadata.csv, cell_by_gene.csv)
   - `ChunkedDataAdapter` - Loads pre-processed data from S3 storage

2. **StandardizedDataset** ([lib/StandardizedDataset.ts](lib/StandardizedDataset.ts)):
   - Normalizes spatial coordinates to [-1, 1] range
   - Provides consistent interface: `spatial`, `embeddings`, `genes`, `clusters`
   - Factory methods: `fromH5ad()`, `fromXenium()`, `fromMerscope()`, `fromS3()`
   - Handles gene expression queries via adapter's matrix interface

#### Single Molecule Data Pipeline

Uses **Web Workers with Comlink** for non-blocking data processing:

1. **File Parsing** ([lib/services/hyparquetService.ts](lib/services/hyparquetService.ts)):
   - `hyparquet` - Pure JavaScript parquet reader with `onPage` callback for column-oriented data
   - `hyparquet-compressors` - Compression codec support (gzip, snappy, etc.)
   - `papaparse` - CSV parsing fallback
   - Memory optimization: Stores column chunks during reading, flattens once at end to avoid repeated array concatenation
   - Warning for files >2GB to alert about potential memory pressure

2. **Web Worker Processing** ([lib/workers/](lib/workers/)):
   - `single-molecule.worker.ts` - Background worker for parsing parquet/CSV files
   - `singleMoleculeWorkerManager.ts` - Singleton pattern for worker lifecycle management
   - Uses Comlink for seamless function proxying across worker boundary
   - Progress callbacks proxied via `Comlink.proxy()` for real-time UI updates
   - Serializable dataset format transfers geneIndex as array of entries
   - Prevents UI freezing during processing of large files (e.g., 26+ seconds for 21M molecules)

3. **SingleMoleculeDataset** ([lib/SingleMoleculeDataset.ts](lib/SingleMoleculeDataset.ts)):
   - Stores molecules as: `uniqueGenes: string[]`, `geneIndex: Map<string, number[]>`, `dimensions: 2|3`
   - Pre-computes normalized coordinates ([-1, 1]) during initialization
   - Gene index stores flattened coordinate arrays: `[x1,y1,z1, x2,y2,z2, ...]`
   - `getCoordinatesByGene()` provides O(1) lookup of pre-normalized coordinates (async for S3 lazy loading)
   - Saves `scalingFactor` for potential denormalization
   - Factory methods:
     - `fromParquet()`, `fromCSV()` - Local file parsing with async progress callbacks (called from worker)
     - `fromS3()` - **Lazy loading from S3** (loads manifest only, genes on-demand)
     - `fromSerializedData()` - Reconstructs dataset from worker-serialized data
   - **Progress Tracking**: Reports progress during file reading (10%), parsing (30%), extraction (50-60%), normalization (60%), indexing (70-85%), filtering (85%), and finalization (90-100%)
   - **Gene Filtering**: Automatically filters control probes and unassigned genes using shared `shouldFilterGene()` utility
   - **Granular Progress**: Indexing loop reports every 5% for real-time feedback on large datasets (millions of molecules)
   - **Responsive UI**: Progress callbacks with `await` yield to event loop for UI updates
   - **Performance Timer**: Displays elapsed time in progress messages and final completion time (e.g., "2.45s", "1m 32.5s")
   - **S3 Lazy Loading**: Downloads gene files on-demand when selected, caches in memory for instant re-selection

4. **Column Mapping Configuration** ([lib/config/moleculeColumnMappings.ts](lib/config/moleculeColumnMappings.ts)):
   - Configurable column names for different dataset types (xenium, merscope, custom)
   - Default mappings: `feature_name`, `x_location`, `y_location`, `z_location` (Xenium)
   - Alternative: `gene`, `global_x`, `global_y`, `global_z` (MERSCOPE)

5. **Gene Filtering Utilities** ([lib/utils/gene-filters.ts](lib/utils/gene-filters.ts)):
   - Shared `shouldFilterGene()` function filters control probes and unassigned genes
   - Used by both single cell (XeniumAdapter) and single molecule (SingleMoleculeDataset) pipelines
   - Filters patterns: negative controls, unassigned, deprecated, codewords, blanks
   - Reduces clutter in gene selection UI and improves performance

#### State Management

Separate stores for each data type ([lib/stores/](lib/stores/)):
- `datasetStore` - Manages single cell datasets (StandardizedDataset)
- `singleMoleculeStore` - Manages single molecule datasets (SingleMoleculeDataset)
- `visualizationStore` - Controls 3D scene state for single cell viewer (camera, colors, filters, gene/celltype selection)
  - **Separate mode states**: `mode` (actual visualization) and `panelMode` (which panel is open)
  - Allows browsing genes/celltypes without changing visualization
  - Visualization updates only when user selects a gene or toggles a celltype
- `singleMoleculeVisualizationStore` - Controls visualization state for single molecule viewer (gene selection with colors, local/global scaling, view mode)
- Uses Zustand for client-side state management

**Routing**: Components check `pathname.startsWith("/sm-viewer")` to determine which store to use:
- `/viewer` - Single cell viewer using `datasetStore` and `visualizationStore`
- `/sm-viewer` - Single molecule viewer using `singleMoleculeStore` and `singleMoleculeVisualizationStore`

### Upload & Storage Flow

#### Single Cell Upload Flow

1. **Client-side**: Calculate dataset fingerprint, check for duplicates via `/api/datasets/check-duplicate/[fingerprint]`
2. **Initiate Upload**: POST to `/api/datasets/initiate` with metadata and file list
   - Creates `Dataset`, `UploadSession`, and `UploadFile` records in PostgreSQL
   - Returns presigned S3 URLs for direct upload
3. **Upload Files**: Client uploads directly to S3 using presigned URLs
4. **Complete**: POST to `/api/datasets/[datasetId]/complete` to finalize upload
5. **Email Notification**: POST to `/api/send-email` with dataset name and metadata
   - Subject line includes dataset name: `"{DatasetName} - Dataset Ready - MERFISHeyes"`
   - Email body includes: total cells, unique genes, platform, cluster columns count, shareable link
6. **View**: Navigate to `/viewer/[datasetId]` which loads data via `StandardizedDataset.fromS3()`

#### Single Molecule Upload Flow

1. **Client-side Processing**:
   - `SingleMoleculeProcessor.processAndUpload()` processes dataset locally
   - Generates fingerprint via `generateSingleMoleculeFingerprint()` (gene names + per-gene molecule counts)
   - Checks for duplicates via `/api/single-molecule/check-duplicate/[fingerprint]`
2. **File Preparation**:
   - Creates `manifest.json.gz` with metadata (genes, dimensions, scaling factor, statistics)
   - Creates `genes/{geneName}.bin.gz` for each gene (gzipped Float32Array of normalized coordinates)
3. **Initiate Upload**: POST to `/api/single-molecule/initiate` with fingerprint, metadata, manifest, and file list
   - Creates `Dataset` (datasetType="single_molecule"), `UploadSession`, and `UploadFile` records
   - Stores manifest JSON in database `manifestJson` field
   - Returns presigned S3 URLs for direct upload
4. **Upload Files**: Client uploads directly to S3 using presigned URLs
   - S3 path: `datasets/{datasetId}/manifest.json.gz` and `datasets/{datasetId}/genes/{gene}.bin.gz`
   - Marks each file complete via `/api/single-molecule/[id]/files/[key]/complete`
5. **Complete**: POST to `/api/single-molecule/[id]/complete` to finalize upload and send email
6. **View**: Navigate to `/sm-viewer/[datasetId]` which lazy-loads data via `SingleMoleculeDataset.fromS3()`

#### Database Schema

Database schema ([prisma/schema.prisma](prisma/schema.prisma)) tracks upload state:
- `datasetType` field distinguishes "single_cell" vs "single_molecule"
- `manifestJson` field stores manifest for single molecule datasets
- `numCells` field stores molecule count for single molecule datasets
- `DatasetStatus` enum: UPLOADING → PROCESSING → COMPLETE/FAILED
- `FileStatus` enum tracks individual file upload status

### 3D Visualization Architecture

#### Single Cell Visualization

**WebGL Rendering** ([lib/webgl/](lib/webgl/)):
- `point-cloud.ts` - Creates Three.js point cloud with custom shaders
- `shaders.ts` - Vertex/fragment shaders for point rendering with color/size control
- `scene-manager.ts` - Manages Three.js scene, camera, controls, and rendering loop
- `visualization-utils.ts` - Color mapping utilities with support for:
  - Gene expression → coolwarm gradient (95th percentile normalization)
  - Categorical clusters → discrete palette colors
  - **Numerical clusters** → coolwarm gradient (same as gene expression)

**Key Features**:
- GPU-accelerated point rendering using BufferGeometry and custom shaders
- Dynamic attribute updates for color (gene expression or cell type) and size
- Coordinate normalization ensures consistent visualization across datasets
- Supports both 2D and 3D spatial coordinates
- **Automatic column type detection**: Columns with ≤100 unique values are categorical, >100 are numerical
- **Numerical cluster visualization**: Numerical metadata columns use gradient coloring instead of discrete categories
- **Interactive hover tooltips**: Mouse over points to see cluster/gene information
  - Celltype mode: Shows original cluster palette color + cluster name (even when filtered/greyed)
  - Gene mode: Shows cluster palette color + cluster name, plus gene gradient color + expression value
  - Numerical clusters: Shows numerical value without color circle
  - Fine-tuned adaptive raycaster threshold (0.1-2.0 based on camera distance) for accurate selection
  - Throttled intersection checking (50ms) for performance
- **Double-click interaction**: Double-click points to toggle cluster selection
  - Automatically switches to celltype mode
  - Toggles cluster in/out of selectedCelltypes for filtering
  - Only works for categorical clusters (numerical clusters excluded)
- **Ref-based state management**: Uses refs to avoid JavaScript closure issues with event handlers

#### Single Molecule Visualization

**Three.js Point Clouds** ([components/single-molecule-three-scene.tsx](components/single-molecule-three-scene.tsx)):
- One point cloud per gene for independent control (easier alpha/visibility management)
- Random HSL color assignment per gene (70-100% saturation, 50-70% lightness for visibility on black background)
- Coordinates pre-scaled by 100x (from normalized [-1,1] to [-100,100] range)
- Point size: `localScale × globalScale × 2.0`
- **Async lazy loading**: `getCoordinatesByGene()` is async to support on-demand S3 loading
- **Gene caching**: Once loaded from S3, genes remain in memory for instant re-selection

**View Modes**:
- **2D Mode**: Top-down orthographic-like view, camera at (0, 0, 200), rotation disabled
- **3D Mode**: Perspective view, camera at (150, 150, 150), rotation enabled
- Camera resets to center when toggling between modes

**Rendering Strategy**:
- Uses `fixed inset-0` positioning to fill entire viewport below navbar
- Black background for high contrast with bright gene colors
- Three.js OrbitControls for pan/zoom/rotate (rotate disabled in 2D mode)
- Multiple point clouds rendered simultaneously (one per selected gene)
- Point cloud update loop wrapped in async function to handle lazy S3 loading

### Component Structure

#### Shared Components
- `components/three-scene.tsx` - Three.js scene for single cell visualization (uses `useEffect` for scene lifecycle)
- `components/visualization-controls.tsx` - UI for gene selection, cluster filtering (single cell viewer)
- `components/dataset-card.tsx` - Dataset preview cards on explore page
- `components/navbar-wrapper.tsx` - Smart routing to appropriate store/modal based on pathname
  - Detects `/sm-viewer` path and shows `SingleMoleculeUploadModal`
  - Detects `/viewer` path and shows `UploadSettingsModal`
  - Manages "Upload & Save" button visibility based on dataset presence

#### Single Cell Components
- `app/viewer/[id]/page.tsx` - Single cell viewer page with dynamic dataset loading from S3
- `components/upload-settings-modal.tsx` - Upload settings for cell datasets (shows point count, genes, clusters)

#### Single Molecule Components
- `app/sm-viewer/page.tsx` - Single molecule viewer page (local file upload, auto-selects 5 random genes on load)
- `app/sm-viewer/[id]/page.tsx` - Single molecule viewer page with S3 lazy loading (auto-selects first 5 genes)
- `components/single-molecule-three-scene.tsx` - Three.js scene for molecule visualization with multiple point clouds
- `components/single-molecule-controls.tsx` - Gene selection UI with search, checkboxes, and colored chips
- `components/view-mode-toggle.tsx` - Toggle button for 2D/3D camera views
- `components/single-molecule-upload-modal.tsx` - Upload settings for molecule datasets (shows molecule count, genes, dimensions, handles S3 upload)

#### File Upload
- `components/file-upload.tsx` - Unified upload component with `singleMolecule` prop
  - **Single Molecule**: Worker-based parsing with granular progress tracking (.parquet/.csv)
  - **Single Cell**: Worker-based parsing with progress callbacks (h5ad, xenium, merscope)
  - Routes to appropriate store based on data type
  - Navigates to `/viewer` for cell data, `/sm-viewer` for molecule data
  - Real-time progress bar (0-100%) and status messages
  - Uses singleton worker managers to keep UI responsive during heavy processing

## Key Technical Details

### Data Format Requirements

#### Single Cell Data

**H5AD**: Requires `obsm['X_spatial']` for coordinates. Optionally supports `obsm['X_umap']` and cell type annotations in `obs`.

**Xenium**: Must have `cells.csv` with centroid columns. Auto-detects cell type columns.

**MERSCOPE**: Must have `cell_metadata.csv` with coordinate columns. Supports `cell_by_gene.csv` for expression matrix.

#### Single Molecule Data

**Parquet**: Columnar format read via `hyparquet` (pure JavaScript, no WASM)
- Requires columns for gene names and x/y/z coordinates
- Column names configurable via `MOLECULE_COLUMN_MAPPINGS`
- Supports both 2D and 3D datasets (z_location optional)
- Processed in web worker to keep UI responsive

**CSV**: Parsed via `papaparse` with same column requirements as parquet
- Automatically infers dimensions based on z column presence
- Less efficient than parquet for large datasets (millions of molecules)

**Dataset Types**:
- **Xenium**: `feature_name`, `x_location`, `y_location`, `z_location`
- **MERSCOPE**: `gene`, `global_x`, `global_y`, `global_z`
- **Custom**: User-defined column mappings

### S3 Loading Architecture

#### Single Molecule S3 Lazy Loading

**API Endpoints**:
- `/api/single-molecule/[id]` - Returns dataset metadata and manifest URL
- `/api/single-molecule/[id]/gene/[geneName]` - Returns presigned URL for specific gene file

**Loading Process** (`SingleMoleculeDataset.fromS3()`):
1. **Fetch Metadata** (10%): GET `/api/single-molecule/[id]` to get manifest URL
2. **Download Manifest** (30%): Fetch and ungzip `manifest.json.gz` from S3
3. **Initialize Dataset** (60%): Create dataset with gene list, empty geneIndex
4. **Lazy Load Genes** (on-demand): Override `getCoordinatesByGene()` to:
   - Check cache first (return immediately if already loaded)
   - Fetch presigned URL from `/api/single-molecule/[id]/gene/[geneName]`
   - Download `genes/{geneName}.bin.gz` from S3
   - Ungzip → Float32Array → number[] conversion
   - Cache in geneIndex for future use

**Gene File Binary Format** (.bin.gz):
- Gzipped Float32Array of normalized coordinates
- Flat array: `[x1, y1, z1, x2, y2, z2, ...]` for 3D or `[x1, y1, x2, y2, ...]` for 2D
- Each coordinate is 4 bytes (Float32)
- Coordinates already normalized to [-1, 1] range
- Gene names sanitized in filenames (special chars → underscores)

**Benefits**:
- ✅ Fast initialization (~1 second for manifest only)
- ✅ Memory efficient (only loads selected genes)
- ✅ Instant re-selection (genes cached after first load)
- ✅ Scales to datasets with hundreds of genes

### Cluster Column Type Detection & Selection

**Type Detection** - All adapters and workers automatically classify columns:
- **Categorical**: ≤100 unique values → discrete color palette, checkbox filtering
- **Numerical**: >100 unique values → coolwarm gradient, no filtering UI

**Implementation**:
- `H5adAdapter.isCategoricalData()` - Analyzes column values during H5AD parsing
- `standardized-dataset.worker.ts` - Uses `isCategoricalData()` helper for Xenium/Merscope
- `ChunkedDataAdapter` - Reads type from S3 metadata, skips palette loading for numerical columns

**Column Selection** - The `selectBestClusterColumn()` utility ([lib/utils/dataset-utils.ts](lib/utils/dataset-utils.ts)) automatically picks the best cluster column with priority:
1. "leiden" column
2. Any column containing "celltype"
3. First categorical column
4. First available column

### Memory Considerations

- Production builds require 4GB RAM by default (`--max-old-space-size=4096`)
- Use `npm run build:low-memory` for servers with <4GB RAM
- Large datasets may cause SIGBUS errors on low-memory systems

**Single Cell Data**:
- Gene expression matrices are cached in `StandardizedDataset.matrix` after first query
- Cell-level aggregation reduces memory footprint vs raw molecule data

**Single Molecule Data**:
- Pre-computed gene index trades memory for O(1) query speed
- For millions of molecules, index stores normalized coordinates per gene
- **Web worker processing** prevents UI freezing during large file parsing
- **Chunked storage** during parquet reading reduces memory pressure (avoids repeated concatenation)
- Parquet files more memory-efficient than CSV parsing
- Warning for files >2GB to alert about potential browser memory limits

### Threading Model

**ALL DATA PROCESSING NOW USES WEB WORKERS** to keep UI responsive:

**Single Cell Adapters** (H5adAdapter, XeniumAdapter, MerscopeAdapter, ChunkedDataAdapter):
- `standardized-dataset.worker.ts` - Background worker for parsing all single cell formats
- `standardizedDatasetWorkerManager.ts` - Singleton pattern for worker lifecycle management
- Workers parse files → load expression matrix → normalize coordinates → serialize
- Main thread receives serialized data → reconstructs `StandardizedDataset` with cached matrix
- Expression matrix pre-loaded in worker for non-S3 datasets
- For S3 datasets: ChunkedDataAdapter created in main thread after worker completes
  - Enables on-demand gene expression loading via `adapter.fetchGeneExpression()`
  - Loads gene data from S3 chunks as needed (avoids loading entire matrix)
- Prevents UI freezing during parsing of large H5AD/Xenium/MERSCOPE files
- ChunkedDataAdapter (S3 loading) also runs in worker with absolute URLs for fetch
- Progress callbacks use `Comlink.proxy()` to avoid DataCloneError in production

**Single Molecule Processing** (hyparquetService, SingleMoleculeDataset):
- `single-molecule.worker.ts` - Background worker for parsing parquet/CSV files
- `singleMoleculeWorkerManager.ts` - Singleton pattern for worker lifecycle management
- Uses Comlink for seamless function proxying across worker boundary
- Progress callbacks proxied via `Comlink.proxy()` for real-time UI updates
- Prevents UI freezing for large files (e.g., 21M+ molecules taking 26+ seconds)
- Worker serializes dataset as JSON-compatible structure for main thread transfer

**Worker Compatibility Fixes**:
- `gzip.ts`: Uses `typeof DecompressionStream !== "undefined"` instead of `window` check
- `ChunkedDataAdapter.ts`: Uses `self.location.origin` for absolute URLs (works in both main thread and workers)

### TypeScript Configuration

- Path alias `@/*` maps to project root for clean imports
- Strict mode enabled
- Target ES2015 with downlevel iteration support
- Next.js plugin for type generation

## Database Schema

Three main entities:
- **Dataset** - Stores dataset metadata, fingerprint for deduplication, status tracking
- **UploadSession** - Tracks multi-file upload progress with expiration
- **UploadFile** - Individual file upload status within a session

All use cascade deletion to maintain referential integrity.

## Tech Stack

- **Framework**: Next.js 15 (App Router) with Turbopack
- **UI**: HeroUI v2 + Tailwind CSS
- **3D Graphics**: Three.js with custom WebGL shaders
- **Data Processing**:
  - Single Cell: h5wasm (WebAssembly H5AD parser), pako (gzip)
  - Single Molecule: hyparquet (pure JS parquet reader), hyparquet-compressors, papaparse, comlink (web workers)
- **State**: Zustand (separate stores for cell vs molecule datasets)
- **Database**: PostgreSQL + Prisma ORM
- **Storage**: AWS S3 with presigned URLs
- **Animation**: Framer Motion, GSAP
- **Background Processing**: Web Workers with Comlink for all data processing

## Performance Optimizations

### Single Cell Data
- **Web Worker Processing**: All parsing (H5AD/Xenium/MERSCOPE) runs in background workers
- **Pre-loaded Expression Matrix**: Matrix cached during worker processing for instant gene queries
- **On-Demand S3 Loading**: ChunkedDataAdapter fetches gene data from S3 as needed
- **Coordinate Validation**: Ensures ≥90% valid coordinates, filters invalid rows automatically
- **Progress Tracking**: Real-time progress updates via Comlink proxying
- **Comlink Proxy Safety**: Progress callbacks properly proxied to avoid DataCloneError in production

### Single Molecule Data
- **Web Worker Processing**: All parsing happens in background worker to prevent UI freezing
- **Comlink Proxying**: Seamless function calls and progress updates across worker boundary
- **Chunked Storage**: During parquet reading, stores column chunks then flattens once at end (avoids repeated array concatenation)
- **Async Progress Callbacks**: Yields to event loop during heavy processing to keep UI responsive
- **Pre-computed Gene Index**: Normalized coordinates computed once during initialization
- **O(1) Gene Lookup**: `getCoordinatesByGene()` is direct Map lookup with no computation
- **Hyparquet Column-Oriented Reading**: Only processes requested columns via `onPage` callback
- **Typed Arrays**: Float32Array for memory-efficient coordinate storage
- **Coordinate Normalization**: [-1, 1] range with saved `scalingFactor` for consistent visualization
- **Granular Progress Updates**: Every 5% during indexing for responsive UI feedback on large datasets
- **Memory Warnings**: Alerts for files >2GB that may cause browser memory pressure
