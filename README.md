# MERFISHEYES

Web-based 3D visualization platform for spatial transcriptomics data. Supports both **single cell** and **single molecule** datasets.

## Features

### Single Cell Visualization
- Multiple format support: .h5ad (AnnData), MERSCOPE, and Xenium
- 3D visualization of cell-level spatial data using Three.js
- Color cells by gene expression or cell type annotations
- **Numerical metadata support**: Continuous metadata (e.g., QC metrics) visualized with gradient coloring
- Interactive filtering and selection of cell populations
- Automatic column type detection (categorical vs numerical)

### Single Molecule Visualization
- Parquet and CSV file support for molecule coordinates
- 3D point cloud visualization with one cloud per gene
- Lazy loading from S3 for efficient memory usage
- Gene-based filtering and multi-gene overlay
- **Automatic control gene filtering**: Removes negative controls, unassigned probes, and codewords
- 2D/3D view mode toggle

### General
- **Web worker processing**: Non-blocking background processing for all data parsing
- **Cloud storage** with AWS S3 integration and lazy loading
- **Duplicate detection** via dataset fingerprinting
- **Email notifications** with shareable links and dataset metadata (cell count, gene count, platform)
- **Dark mode**
- Works on desktop and tablet

## Tech Stack

- [Next.js 15](https://nextjs.org/) - React framework with Turbopack
- [HeroUI v2](https://heroui.com/) - UI components
- [Three.js](https://threejs.org/) - 3D visualization
- [TypeScript](https://www.typescriptlang.org/) - Type safety
- [Tailwind CSS](https://tailwindcss.com/) - Styling
- [Zustand](https://zustand-demo.pmnd.rs/) - State management
- [Prisma](https://www.prisma.io/) - Database ORM (PostgreSQL)
- [AWS S3](https://aws.amazon.com/s3/) - Cloud file storage
- [h5wasm](https://github.com/usnistgov/h5wasm) - HDF5/H5AD file reading
- [Hyparquet](https://github.com/hyparam/hyparquet) - Pure JavaScript parquet parsing
- [Comlink](https://github.com/GoogleChromeLabs/comlink) - Web worker communication
- [Pako](https://github.com/nodeca/pako) - Gzip compression/decompression
- [PapaParse](https://www.papaparse.com/) - CSV parsing
- [SendGrid](https://sendgrid.com/) - Email notifications

## Getting Started

Requires Node.js 18+

### Installation

```bash
git clone <repository-url>
cd merfisheyes-heroui
npm install
```

### Environment Setup

Copy `.env.example` to `.env.local` and configure:

```bash
cp .env.example .env.local
```

Required environment variables:
- `DATABASE_URL` - PostgreSQL connection string
- `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_REGION`, `AWS_S3_BUCKET` - S3 storage
- `NEXT_PUBLIC_BASE_URL` - Base URL for the application
- `SENDGRID_API_KEY`, `SENDGRID_FROM_EMAIL` - Email notifications

See `.env.example` for full list and examples.

### Database Setup

```bash
npx prisma generate    # Generate Prisma client
npx prisma migrate dev # Run database migrations
```

### Development

```bash
npm run dev            # Start dev server (http://localhost:3000)
```

### Production Build

```bash
npm run build          # Build for production (requires 4GB RAM)
npm start              # Start production server
```

**Low memory servers:** If you encounter SIGBUS errors on servers with limited RAM:
```bash
npm run build:low-memory  # Build with 2GB memory limit
```

Or manually set memory limit:
```bash
NODE_OPTIONS='--max-old-space-size=2048' npm run build
```

### Deployment

For deploying to remote servers:

1. **Build locally** on a high-memory machine
2. **Deploy** using the included script:
   ```bash
   ./deploy.sh
   ```

The deploy script:
- Builds the project locally
- Transfers `.next`, `public`, `package.json`, and `prisma` to remote server
- Runs `npx prisma generate` on production
- Restarts the application with PM2

**Note**: The script does NOT transfer `.env.local` - production server should have its own environment variables configured.

## Usage

### Uploading Data

#### Single Cell Data
Navigate to `/viewer` and upload:
- **H5AD**: Single `.h5ad` file
- **Xenium**: Folder with `cells.csv` and related files
- **MERSCOPE**: Folder with `cell_metadata.csv` and related files

#### Single Molecule Data
Navigate to `/sm-viewer` and upload:
- **Parquet**: `.parquet` file with columns for gene names and x/y/z coordinates
- **CSV**: `.csv` file with same column structure
- Supports Xenium and MERSCOPE column naming conventions

Drag and drop or click to upload. After processing, you'll receive an email with a shareable link to view your dataset.

### Viewer Controls

#### Single Cell Viewer (`/viewer/[id]`)
- **Rotate**: Left click + drag
- **Pan**: Right click + drag or middle click + drag
- **Zoom**: Mouse wheel
- **Filter**: Use side panel to filter by cell type (categorical columns only)
- **Color**: Select gene from dropdown to color by expression, or choose cluster column:
  - **Categorical columns** (≤100 unique values): Discrete colors with checkbox filtering
  - **Numerical columns** (>100 unique values): Coolwarm gradient, no filtering UI

#### Single Molecule Viewer (`/sm-viewer/[id]`)
- **Rotate**: Left click + drag (disabled in 2D mode)
- **Pan**: Right click + drag
- **Zoom**: Mouse wheel
- **Select Genes**: Search and check genes to display
- **2D/3D Toggle**: Switch between top-down and perspective views
- **Scale**: Adjust point size with global and per-gene local scales

## Data Format Requirements

### Single Cell Formats

#### H5AD
- Standard AnnData format
- Requires `obsm['X_spatial']` for coordinates
- Optional: `obsm['X_umap']`, cell type annotations in `obs`

#### Xenium (Cell-level)
- Required: `cells.csv` with centroids
- Optional: `transcripts.csv`, `features.tsv`
- Detects cell type columns automatically

#### MERSCOPE (Cell-level)
- Required: `cell_metadata.csv` with coordinates
- Optional: `cell_categories.csv`, `cell_numeric_categories.csv`, `cell_by_gene.csv`

### Single Molecule Formats

#### Parquet
- Columnar binary format (most efficient for large datasets)
- Required columns (configurable):
  - Gene name: `feature_name` (Xenium) or `gene` (MERSCOPE)
  - X coordinate: `x_location` (Xenium) or `global_x` (MERSCOPE)
  - Y coordinate: `y_location` (Xenium) or `global_y` (MERSCOPE)
  - Z coordinate: `z_location` (Xenium) or `global_z` (MERSCOPE) - optional for 2D data

#### CSV
- Text format with same column requirements as parquet
- Automatically infers 2D vs 3D based on z column presence
- Less memory-efficient than parquet for large datasets (millions of molecules)

## API Routes

The application provides RESTful API endpoints for dataset upload and management:

### Single Cell Endpoints

| Route | Method | Purpose |
|-------|--------|---------|
| `/api/datasets/check-duplicate/{fingerprint}` | GET | Check if dataset already exists |
| `/api/datasets/initiate` | POST | Start upload, get presigned S3 URLs |
| `/api/datasets/{datasetId}/complete` | POST | Mark upload as complete |
| `/api/datasets/{datasetId}` | GET | Get dataset info + download URLs |

### Single Molecule Endpoints

| Route | Method | Purpose |
|-------|--------|---------|
| `/api/single-molecule/check-duplicate/{fingerprint}` | GET | Check if single molecule dataset exists |
| `/api/single-molecule/initiate` | POST | Start upload, get presigned S3 URLs |
| `/api/single-molecule/{id}/files/{key}/complete` | POST | Mark individual file as uploaded |
| `/api/single-molecule/{id}/complete` | POST | Finalize upload, send email |
| `/api/single-molecule/{id}` | GET | Get dataset metadata and manifest URL |
| `/api/single-molecule/{id}/gene/{geneName}` | GET | Get presigned URL for specific gene file |

### Upload Flow

#### Single Cell Upload
1. **Check for duplicates** - `GET /api/datasets/check-duplicate/{fingerprint}`
2. **Initiate upload** - `POST /api/datasets/initiate` with metadata and file list
   - Creates database records (Dataset, UploadSession, UploadFile)
   - Returns presigned S3 URLs for file upload
3. **Upload files** - Use presigned URLs to upload directly to S3
4. **Complete upload** - `POST /api/datasets/{datasetId}/complete`
   - Finalizes the upload session, sends email notification

#### Single Molecule Upload
1. **Process locally** - Client processes dataset into manifest + gene files
2. **Check for duplicates** - `GET /api/single-molecule/check-duplicate/{fingerprint}`
3. **Initiate upload** - `POST /api/single-molecule/initiate`
   - Returns presigned S3 URLs for manifest and all gene files
4. **Upload files** - Upload `manifest.json.gz` and `genes/{gene}.bin.gz` files to S3
5. **Mark files complete** - `POST /api/single-molecule/{id}/files/{key}/complete` for each file
6. **Complete upload** - `POST /api/single-molecule/{id}/complete`
   - Sends email with link to `/sm-viewer/{id}`

## Project Structure

```
├── app/                           # Next.js app directory
│   ├── api/                      # API routes
│   │   ├── datasets/             # Single cell endpoints
│   │   ├── single-molecule/      # Single molecule endpoints
│   │   └── send-email*/          # Email notification services
│   ├── viewer/                   # Single cell viewer
│   │   └── [id]/                 # Dynamic dataset routes
│   ├── sm-viewer/                # Single molecule viewer
│   │   └── [id]/                 # Dynamic dataset routes with S3 lazy loading
│   ├── explore/                  # Example datasets page
│   └── about/                    # About page
├── components/                   # React components
│   ├── three-scene.tsx           # Single cell Three.js scene
│   ├── single-molecule-three-scene.tsx  # Single molecule Three.js scene
│   ├── visualization-controls.tsx       # Single cell controls
│   ├── single-molecule-controls.tsx     # Single molecule controls
│   └── file-upload.tsx           # Unified upload component
├── lib/                          # Core logic
│   ├── adapters/                # Single cell format adapters
│   │   ├── H5adAdapter.ts
│   │   ├── XeniumAdapter.ts
│   │   ├── MerscopeAdapter.ts
│   │   └── ChunkedDataAdapter.ts  # S3 loading adapter
│   ├── workers/                 # Web workers for background processing
│   │   ├── standardized-dataset.worker.ts  # Single cell parsing worker (H5AD/Xenium/MERSCOPE)
│   │   ├── standardizedDatasetWorkerManager.ts  # Single cell worker manager
│   │   ├── single-molecule.worker.ts  # Parquet/CSV parsing worker
│   │   └── singleMoleculeWorkerManager.ts  # Single molecule worker manager
│   ├── stores/                  # Zustand state stores
│   │   ├── datasetStore.ts      # Single cell datasets
│   │   ├── singleMoleculeStore.ts  # Single molecule datasets
│   │   ├── visualizationStore.ts   # Single cell viz state
│   │   └── singleMoleculeVisualizationStore.ts  # Single molecule viz state
│   ├── services/                # Data processing services
│   │   └── hyparquetService.ts  # Hyparquet parquet reader
│   ├── utils/
│   │   ├── SingleMoleculeProcessor.ts  # S3 upload processing
│   │   ├── fingerprint.ts       # Dataset fingerprinting
│   │   └── gene-filters.ts      # Shared gene filtering (control probes, etc.)
│   ├── webgl/                   # WebGL/Three.js utilities (single cell)
│   ├── s3.ts                    # S3 client utilities
│   ├── prisma.ts                # Database client
│   ├── StandardizedDataset.ts   # Single cell dataset class
│   └── SingleMoleculeDataset.ts # Single molecule dataset class (lazy loading)
├── prisma/                      # Database schema
│   └── schema.prisma            # Supports both dataset types
└── public/                      # Static assets
```

## Contributing

Pull requests welcome.

## Acknowledgments

Built at the Bogdan Bintu Lab, UCSD.

## License

MIT
