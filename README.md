# MERFISH Eyes

Web-based 3D visualization for spatial transcriptomics data. Supports .h5ad, MERSCOPE, and Xenium formats.

## Features

- Multiple format support: .h5ad (AnnData), MERSCOPE, and Xenium
- 3D visualization of spatial data using Three.js
- Color cells by gene expression or cell type annotations
- Interactive filtering and selection of cell populations
- Dark mode
- Works on desktop and tablet

## Tech Stack

- [Next.js 15](https://nextjs.org/)
- [HeroUI v2](https://heroui.com/)
- [Three.js](https://threejs.org/)
- [TypeScript](https://www.typescriptlang.org/)
- [Tailwind CSS](https://tailwindcss.com/)
- [Zustand](https://zustand-demo.pmnd.rs/)
- [React Toastify](https://fkhadra.github.io/react-toastify/)
- [Prisma](https://www.prisma.io/) - Database ORM
- [AWS S3](https://aws.amazon.com/s3/) - File storage

## Getting Started

Requires Node.js 18+

```bash
git clone <repository-url>
cd merfisheyes-heroui
npm install
npm run dev
```

Open http://localhost:3000

### Production Build

```bash
npm run build
npm start
```

**Low memory servers:** If you encounter SIGBUS errors on servers with limited RAM:
```bash
npm run build:low-memory
```

Or manually set memory limit:
```bash
NODE_OPTIONS='--max-old-space-size=2048' npm run build
```

## Usage

### Uploading Data

Upload from the home page. Supports:
- H5AD: Single `.h5ad` file
- Xenium: Folder with `cells.csv` and related files
- MERSCOPE: Folder with `cell_metadata.csv` and related files

Drag and drop or click to upload. The viewer loads automatically after processing.

### Controls

- Rotate: Left click + drag
- Pan: Right click + drag or middle click + drag
- Zoom: Mouse wheel
- Filter cell types using the side panel
- Color by gene expression using the gene dropdown

## Data Format Requirements

### H5AD
- Standard AnnData format
- Requires `obsm['X_spatial']` for coordinates
- Optional: `obsm['X_umap']`, cell type annotations in `obs`

### Xenium
- Required: `cells.csv` with centroids
- Optional: `transcripts.csv`, `features.tsv`
- Detects cell type columns automatically

### MERSCOPE
- Required: `cell_metadata.csv` with coordinates
- Optional: `cell_categories.csv`, `cell_numeric_categories.csv`, `cell_by_gene.csv`

## API Routes

The application provides RESTful API endpoints for dataset upload and management:

| Route | Method | Purpose |
|-------|--------|---------|
| `/api/datasets/check-duplicate/{fingerprint}` | GET | Check if dataset already exists |
| `/api/datasets/initiate` | POST | Start upload, get presigned S3 URLs |
| `/api/datasets/{datasetId}/complete` | POST | Mark upload as complete |
| `/api/datasets/{datasetId}` | GET | Get dataset info + download URLs |

### Upload Flow

1. **Check for duplicates** - `GET /api/datasets/check-duplicate/{fingerprint}`
2. **Initiate upload** - `POST /api/datasets/initiate` with metadata and file list
   - Creates database records (Dataset, UploadSession, UploadFile)
   - Returns presigned S3 URLs for file upload
3. **Upload files** - Use presigned URLs to upload directly to S3
4. **Complete upload** - `POST /api/datasets/{datasetId}/complete`
   - Finalizes the upload session
   - Sets `completedAt` timestamp and `manifestUrl`

## Project Structure

```
├── app/                    # Next.js app directory
│   ├── api/               # API routes
│   │   └── datasets/      # Dataset management endpoints
│   ├── viewer/            # 3D visualization page
│   ├── explore/           # Example datasets page
│   └── about/             # About page
├── components/            # React components
│   ├── three-scene.tsx    # Three.js scene component
│   ├── visualization-controls.tsx
│   └── dataset-card.tsx
├── lib/                   # Core logic
│   ├── adapters/         # Data format adapters
│   │   ├── H5adAdapter.ts
│   │   ├── XeniumAdapter.ts
│   │   └── MerscopeAdapter.ts
│   ├── stores/           # Zustand state stores
│   ├── webgl/            # WebGL/Three.js utilities
│   ├── s3.ts             # S3 client utilities
│   ├── prisma.ts         # Database client
│   └── StandardizedDataset.ts
├── prisma/                # Database schema
└── public/               # Static assets

```

## Contributing

Pull requests welcome.

## Acknowledgments

Built at the Bintu Lab, Stanford University.

## License

MIT
