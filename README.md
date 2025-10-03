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

## Project Structure

```
├── app/                    # Next.js app directory
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
│   └── StandardizedDataset.ts
└── public/               # Static assets

```

## Contributing

Pull requests welcome.

## Acknowledgments

Built at the Bintu Lab, Stanford University.

## License

MIT
