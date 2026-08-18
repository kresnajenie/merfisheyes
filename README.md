# MERFISHEYES

Web-based 3D visualization platform for spatial transcriptomics data. Supports both **single cell** and **single molecule** datasets, viewable directly in the browser or processed on the server for datasets too large for a browser to handle.

> Live at **[merfisheyes.com](https://www.merfisheyes.com)**.

---

## What files can I bring? (input requirements)

MERFISHEYES reads the **raw outputs** of the common spatial platforms — you don't need to pre-format anything. Below is exactly which files each type needs and what they should look like. Drop a whole folder or a single file on the homepage; the app figures out the format from the file names.

### Single cell

<table>
<tr><th>Format</th><th>What to drop</th><th>Required</th><th>Also read if present</th></tr>
<tr>
<td><b>H5AD</b><br>(AnnData)</td>
<td>A single <code>.h5ad</code> file</td>
<td>Spatial coordinates in <code>obsm['X_spatial']</code></td>
<td><code>obsm['X_umap']</code> / other embeddings, and any cell-type / metadata columns in <code>obs</code></td>
</tr>
<tr>
<td><b>Xenium</b><br>(10x)</td>
<td>The Xenium output <b>folder</b></td>
<td><code>cells.csv</code> or <code>cells.csv.gz</code> (cell centroids)</td>
<td>An expression matrix — <code>cell_feature_matrix.h5</code>, or <code>cell_feature_matrix/</code> (<code>matrix.mtx.gz</code> + <code>features.tsv.gz</code> + <code>barcodes.tsv.gz</code>), or <code>cell_by_gene.csv</code>. A <code>transcripts.parquet</code>/<code>.csv</code> is picked up as single-molecule data (see combined upload below)</td>
</tr>
<tr>
<td><b>MERSCOPE</b><br>(Vizgen)</td>
<td>The MERSCOPE output <b>folder</b></td>
<td><code>cell_metadata.csv</code> (per-cell coordinates)</td>
<td><code>cell_by_gene.csv</code> (expression matrix)</td>
</tr>
<tr>
<td><b>Pre-chunked</b></td>
<td>A folder produced by <a href="scripts/README.md">the Python scripts</a></td>
<td><code>manifest.json</code> + <code>coords/</code>, <code>expr/</code>, <code>obs/</code>, <code>palettes/</code></td>
<td>—</td>
</tr>
</table>

**Cell-type annotations (optional).** Cell types don't have to live in the file. Run [MapMyCells](https://knowledge.brain-map.org/mapmycells/process) (or any tool) yourself and drop the resulting **one-row-per-cell CSV** alongside the dataset; it's merged in as another cluster column, with its own palette and differential-expression stats. There is no server-side cell-typing — you stay in control of the labels.

### Single molecule

<table>
<tr><th>Format</th><th>What to drop</th><th>Required columns</th></tr>
<tr>
<td><b>Parquet</b> (preferred for large data)<br>or <b>CSV</b></td>
<td>One <code>.parquet</code> / <code>.csv</code> of per-molecule rows</td>
<td>A <b>gene</b> column, an <b>x</b> and <b>y</b> column, and optionally a <b>z</b> column (omit for 2D data). A <b>cell-id</b> column is optional.</td>
</tr>
</table>

Default column names are auto-detected per platform, and you can remap them at upload time:

| | Gene | X | Y | Z | Cell ID |
|---|---|---|---|---|---|
| **Xenium** | `feature_name` | `x_location` | `y_location` | `z_location` | `cell_id` (string, `UNASSIGNED` = no cell) |
| **MERSCOPE** | `gene` | `global_x` | `global_y` | `global_z` | `cell_id` (`-1` = no cell) |
| **Custom** | *(you pick each column from a dropdown)* | | | | |

Parquet is strongly preferred for millions of molecules (columnar, far less memory than CSV).

### Combined single-cell + single-molecule

Drop a **single-cell folder that also contains a transcripts file** (a typical Xenium export) and MERFISHEYES offers to upload **both** at once: the cells and the molecules become two datasets grouped into one project, with the molecules overlaid on the cells.

---

## Two ways to load a dataset

On the homepage, pick **Single Cell** or **Single Molecule**, then choose one of:

### 1. Preview in browser
The file is parsed **entirely in your browser** (h5wasm / hyparquet in web workers) and opens instantly — nothing leaves your computer. Best for quick looks and datasets that fit in browser memory (roughly < 500K cells / tens of millions of molecules).

### 2. Upload & process on server
For large datasets, upload the **raw bytes** and let a background worker do the processing. Requires signing in (so we know where to email you). The flow:

1. Drop your file/folder → confirm the dataset name.
2. **Single-molecule uploads add a "Confirm columns" step** — a preview of the first rows with dropdowns to map gene / x / y / z / cell-id (auto-filled from the detected platform). This is what guarantees the server reads your columns correctly.
3. A full-screen progress bar tracks the upload; you can leave once it says "Uploaded — we'll email you."
4. An AWS worker runs the same Python processors used for pre-chunking, writes the result to S3, and **emails you a link** to the viewer when it's done.

This path has no browser memory limits and produces the exact same chunked format as everything else. See **[docs/server-side-ingestion/](docs/server-side-ingestion/README.md)** for the full architecture.

---

## Accounts, projects & community

Signing in (**Google**, via NextAuth) unlocks:

- **Your datasets** (`/account`) — everything you've uploaded, with the same cards as Explore, plus edit controls for titles, descriptions, thumbnails, and visibility.
- **Projects** (`/account/projects/[id]`) — group related datasets (e.g. the cells + molecules from one experiment) into a folder with its own page and thumbnail.
- **Community** — submit a dataset or project for review; once an admin approves it, it appears in the **Community** view on **[Explore](https://www.merfisheyes.com/explore)** for everyone.

Explore (`/explore`) lists curated and community datasets with a persistent **Featured** row on top; community projects open a detail page (`/explore/[id]`).

---

## Viewer

Datasets open at `/viewer/[id]` (single cell) or `/sm-viewer/[id]` (single molecule). You can also open local or S3 data directly via `/viewer/from-local`, `/viewer/from-s3`, and the `sm-viewer` equivalents.

### Single cell viewer

- **Loading progress**: real-time bar while the dataset streams from S3.
- **Navigate**: rotate (left-drag), pan (right/middle-drag), zoom (wheel).
- **Hover** a point for a tooltip (cluster + gene value, keeps original palette colors even when filtered); **double-click** a point to toggle that cluster's selection.
- **Color by gene** (coolwarm gradient, auto-scaled to the 95th percentile) or by a **cluster column**:
  - *Categorical* (≤100 unique values): discrete palette + checkbox filtering.
  - *Numerical* (>100 unique values): coolwarm gradient, no filtering.
- **Combined mode**: select a gene **and** toggle celltypes to show the gene gradient only on those populations (others greyed).
- **Gene / numerical scalebar**: draggable min/max scrubbers, auto-scaled per selection.
- **Legends panel** (top-right): active gene + celltype badges, embedded scalebar, one-click removal, "Clear All".
- **Plot panel**: a floating, draggable, resizable panel with quantification plots — boxplots and histograms of gene expression per celltype, celltype bar charts, and numerical histograms, with an optional secondary "group by" column, Top-N / axis caps, outlier and density toggles, and CSV export.

### Single molecule viewer

- **Auto-selection**: the first few genes are shown on load.
- **Select genes**: search + check; each gene is an independently colored point cloud, lazy-loaded from S3 and cached.
- **2D / 3D toggle**: top-down orthographic view vs. full 3D rotation.
- **Scale**: global and per-gene point sizes.

### Split screen

Compare two datasets side by side (single-cell or single-molecule, in any combination), with camera and selection kept in sync between the panels.

---

## Feature summary

- Single cell: H5AD / Xenium / MERSCOPE, gene + celltype + numerical coloring, combined mode, interactive legends & scalebar, quantification plots.
- Single molecule: parquet/CSV, one point cloud per gene, S3 lazy loading, 2D/3D, automatic control-probe/unassigned filtering.
- Two ingestion paths: **in-browser** (web workers) and **server-side** (AWS Batch worker) with email delivery.
- Column remap step for single-molecule server uploads.
- Accounts, projects, and an admin-reviewed community catalog.
- Split-screen comparison, duplicate detection via fingerprinting, dark mode, desktop + tablet.

---

## Tech Stack

- [Next.js 15](https://nextjs.org/) (App Router, Turbopack) + [TypeScript](https://www.typescriptlang.org/)
- [HeroUI v2](https://heroui.com/) + [Tailwind CSS](https://tailwindcss.com/) — UI
- [Three.js](https://threejs.org/) with custom WebGL shaders — 3D rendering
- [Zustand](https://zustand-demo.pmnd.rs/) — state; [Plotly](https://plotly.com/javascript/) + [react-rnd](https://github.com/bokuweb/react-rnd) — plot panel
- [Prisma](https://www.prisma.io/) + PostgreSQL (Supabase) — database
- [NextAuth v5](https://authjs.dev/) — Google sign-in
- [AWS S3](https://aws.amazon.com/s3/) — file storage; [AWS Batch](https://aws.amazon.com/batch/) on Fargate — server-side processing worker; [AWS SES](https://aws.amazon.com/ses/) — email
- [h5wasm](https://github.com/usnistgov/h5wasm) — H5AD; [hyparquet](https://github.com/hyparam/hyparquet) — parquet; [PapaParse](https://www.papaparse.com/) — CSV; [Comlink](https://github.com/GoogleChromeLabs/comlink) — web workers; [pako](https://github.com/nodeca/pako) — gzip

---

## Getting Started (development)

Requires Node.js 18+.

```bash
git clone <repository-url>
cd merfisheyes-heroui
npm install
```

### Environment

Copy `.env.example` to `.env.local` and fill in:

- `DATABASE_URL` — PostgreSQL connection string
- `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_REGION`, `AWS_S3_BUCKET` — S3 storage
- `NEXT_PUBLIC_BASE_URL` — base URL of the app
- Google OAuth (`AUTH_GOOGLE_ID`, `AUTH_GOOGLE_SECRET`) + `AUTH_SECRET` — sign-in
- SES sender config — email notifications
- Server-side ingestion (optional, for the Upload-&-process-on-server path): AWS Batch job queue / definition, `CALLBACK_SECRET`. See [docs/server-side-ingestion/](docs/server-side-ingestion/README.md).

### Database

```bash
npx prisma generate     # generate the Prisma client
npx prisma migrate dev  # apply migrations
```

### Run

```bash
npm run dev             # dev server at http://localhost:3000
```

### Production build

```bash
npm run build           # production build (needs ~4GB RAM)
npm start               # start the production server
```

Low-memory servers (avoids SIGBUS on limited RAM):

```bash
npm run build:low-memory   # 2GB memory limit
```

---

## Testing

Unit tests (Vitest), end-to-end browser tests (Playwright), and a performance regression suite. See **[TESTING.md](TESTING.md)** for the full guide.

```bash
npm test                 # unit tests
npm run e2e:install      # one-time: download the Playwright browser
npm run test:e2e:quick   # fast E2E smoke suite over the tiny test datasets
npm run test:perf        # performance suite (compares against perf/baselines/)
```

---

## Processing large datasets outside the browser

Two options for datasets too big to parse in a browser:

1. **Upload & process on server** (in-app) — described above; a background worker does it for you.
2. **Pre-process locally** with the Python scripts and upload the ready-made chunks:

   ```bash
   python scripts/process_spatial_data.py    data.h5ad   output/   # single cell
   python scripts/process_single_molecule.py data.parquet output/  # single molecule
   ```

   See **[scripts/README.md](scripts/README.md)** for formats, flags, and the output layout. (BIL HPC / SLURM pipelines live in [scripts/bil-scripts/](scripts/bil-scripts/README.md).)

---

## Project structure

```
├── app/                       # Next.js App Router
│   ├── api/                   # API routes
│   │   ├── datasets/          #   single-cell browser upload
│   │   ├── single-molecule/   #   single-molecule browser upload
│   │   ├── ingest/            #   server-side ingestion (initiate/complete/callback/status)
│   │   ├── projects/          #   projects (grouping)
│   │   └── admin/             #   admin: processing dashboard, community review, users
│   ├── viewer/ , sm-viewer/   # viewers (+ [id], from-local, from-s3)
│   ├── explore/               # curated + community catalog (+ [id], bil/[bilCode])
│   ├── account/               # your datasets & projects (+ edit)
│   ├── admin/                 # admin dashboard & catalog
│   └── auth/signin/
├── components/                # React components (viewers, controls, plots, upload, cards)
├── lib/
│   ├── adapters/              # H5ad / Xenium / Merscope / ChunkedData adapters
│   ├── ingest/                # server-ingestion helpers (classify-folder, error classification, …)
│   ├── stores/                # Zustand stores (+ per-panel factories for split screen)
│   ├── workers/               # web workers for parsing
│   ├── webgl/                 # Three.js / shader utilities
│   ├── config/                # visualization + column-mapping config
│   ├── ses.ts                 # email
│   ├── StandardizedDataset.ts # single-cell dataset class
│   └── SingleMoleculeDataset.ts
├── worker/                    # server-side ingestion worker (Dockerfile + entrypoint.py)
├── scripts/                   # Python preprocessing (+ bil-scripts/ HPC pipelines)
├── prisma/                    # schema + migrations
└── docs/                      # design & ops docs (server-side-ingestion/, …)
```

---

## Contributing

Pull requests welcome. PRs target the `develop` branch.

## Acknowledgments

Built at the Bogdan Bintu Lab, UCSD.

## License

MIT
