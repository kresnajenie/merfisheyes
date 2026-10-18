# MERFISHEYES

Web-based 3D visualization platform for spatial transcriptomics data. Supports both **single cell** and **single molecule** datasets, viewable directly in the browser or processed on the server for datasets too large for a browser to handle.

> Live at **[merfisheyes.com](https://www.merfisheyes.com)** — user guide at **[merfisheyes.com/docs](https://www.merfisheyes.com/docs)**.

---

## What files can I bring? (input requirements)

MERFISHEYES reads the **raw outputs** of the common spatial platforms — you don't need to pre-format anything. Below is exactly which files each type needs and what they should look like. Drop a whole folder or a single file on the homepage; the app figures out the format from the file names.

There are three kinds of data, and they can be shown together: **cells**, **molecules**, and an **image** (e.g. a DAPI mosaic) drawn under them.

| I have… | Go to |
|---|---|
| Raw platform output (H5AD, Xenium, MERSCOPE, a transcripts file) | [Single cell](#single-cell), [Single molecule](#single-molecule), [Combined upload](#combined-single-cell--single-molecule--image-raw-data) |
| An image | [Image](#image) |
| Data that is already in a format the viewer reads — a SpatialData store, or merfisheyes' own output folders | [Already-processed data](#already-processed-data-opens-as-it-is) |

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

### Image

An image is drawn under the cells or the molecules, or opened on its own. It must be an **OME-Zarr multiscale image**:

| | Requirement |
|---|---|
| **Format** | OME-Zarr 0.5 on **Zarr v3** (a folder with a `zarr.json` at its root). Zarr v2 and TIFF are not read. |
| **Axes** | Exactly `c, y, x` — one or more channels, no z and no time. |
| **Pixels** | `uint8` or `uint16`. |
| **Pyramid** | Several resolution levels, each half the one before. The coarsest should fit in a chunk or two: it is loaded whole as the backdrop. |
| **Placement** | Each level's `scale` and `translation` must put the image in the **same coordinates as the cells** (usually microns). Rotations and flips are not applied, so they must already be in the pixels. |
| **Channel names** | Optional, from `omero.channels[].label`. |

```
dapi.zarr/
├── zarr.json        the multiscale metadata (axes, scale + translation per level)
├── 0/               full resolution: zarr.json + c/<channel>/<row>/<column> chunks
├── 1/               half resolution
└── …
```

`docs/spatialdata/README.md` has a [short Python recipe](docs/spatialdata/README.md#writing-an-image) for writing one from a stitched array.

**Where an image can come from:**

- **A folder on your computer** — drop it on the single-cell **Folder** card ("Preview in browser"). Tiles are read from disk as you zoom; nothing is uploaded.
- **A public URL** — paste it into **Load from S3** on the homepage, or open `/image-viewer/from-s3?url=<the .zarr folder's URL>`. The bucket must allow the site in CORS.
- **Beside your cells**, in a [combined upload](#combined-single-cell--single-molecule--image-raw-data) or an [already-processed folder](#already-processed-data-opens-as-it-is).

**What you set on an image is kept in the link.** Which channels are on, their colours and their lower and upper bounds are written into the page's URL as you change them, so a copied link opens looking the same.

An image is only drawn under cells that kept their **original coordinates**. Folders written by older versions of the Python scripts normalized them to [-1, 1]; those load fine, but an image cannot be placed under them.

### Combined single-cell + single-molecule + image (raw data)

With **Upload & process on server**, drop a **single-cell folder that also contains a transcripts file** (a typical Xenium export) and/or an **OME-Zarr image folder** (any name ending `.zarr`), and MERFISHEYES offers to upload them together: each becomes its own dataset, grouped into one project, and the cells open with the molecules over them and the image under them.

```
sample/
├── cells.h5ad            or a Xenium / MERSCOPE export
├── transcripts.parquet   (optional) the molecules
└── dapi.zarr/            (optional) the image
```

In **Preview in browser**, raw cells are shown on their own: overlaying needs the other parts to be [already processed](#already-processed-data-opens-as-it-is).

### Already-processed data (opens as it is)

If your data is already in a format the viewer reads, nothing needs converting. Drop the folder on the single-cell **Folder** card in **Preview in browser** and it opens at once, whatever its size — each part is read from disk a piece at a time, as you look at it. Cells, molecules and an image in the same folder are overlaid.

Two layouts are accepted. **Folder names don't matter**: each part is recognised by what is inside it.

**A. A SpatialData store, with the molecules beside it**

```
sample/
├── anything.zarr/     SpatialData store (Zarr v3): cells in tables/, the image in images/
├── anything/          merfisheyes molecule folder: manifest.json.gz + genes/
└── annotations/       (optional) clusters + DEGs from scripts/spatialdata_cluster.py
```

**B. merfisheyes' own output folders**

```
sample/
├── anything/          cells: manifest.json + expr/ + coords/ (+ obs/, palettes/)
│                      from scripts/process_spatial_data.py
├── anything/          molecules: manifest.json.gz + genes/
│                      from scripts/process_single_molecule.py
└── anything.zarr/     the image: an OME-Zarr image (see above)
```

| Part | Recognised by |
|---|---|
| Cells (SpatialData) | a `zarr.json` plus `tables/`, `images/` or another SpatialData element group |
| Cells (AnnData table, Zarr v3) | a `zarr.json` with `obs/` and `var/` beside it — a table written by `write_zarr`, or copied out of a store's `tables/` |
| Cells (merfisheyes) | `manifest.json`, `expr/index.json` and `coords/spatial.bin.gz` |
| Molecules | `manifest.json.gz` (or `manifest.json`) and `genes/*.bin.gz` |
| Image | a `zarr.json` describing an OME-Zarr image with axes `c, y, x` |

- **Any part may be missing.** Cells + image, cells + molecules, molecules + image, or a single part on its own all open. Without cells the molecules are the scene, with the image under them; an image alone opens in the image viewer.
- **A part is a folder directly inside the one you drop**, or the dropped folder itself.
- **A table on its own works in either Zarr version.** A Zarr v2 AnnData folder (`.zgroup`) opens as before; a Zarr v3 one (`zarr.json`), including `tables/<table>/` or `tables/` taken out of a SpatialData store, opens with the same reader a store uses.
- **With several candidates for one part** (two images, say) the first by name is used, and a message says which were skipped.
- **A SpatialData store's molecules** are not read from its `points/`; convert them once with `scripts/spatialdata_points_to_sm.py` and put the result beside the store.
- **Upload & Save** from this preview uploads every part exactly as it is — nothing is processed — each as its own dataset in one project, and the saved cells open with the same overlays. (Molecules and an image saved without cells are grouped in a project but not overlaid.) Signing in is required; no single file may exceed 5 GB.

**The same data on a public bucket** opens from **Load from S3** on the homepage. Paste the URL of a merfisheyes folder, a SpatialData store, an OME-Zarr image, or the folder holding a store with its molecules beside it — the last opens with everything overlaid, and the link you land on reopens all of it. A bucket's folders can't be listed from a browser, so the folder's `mapping.json` names its store (`"cells": "name.zarr"`); without that, a store named after the folder (`sample/sample.zarr`) is tried. A store nobody owns yet shows **Add to my account**, like any other dataset opened by URL. See [docs/spatialdata](docs/spatialdata/README.md) for the layout.

What is and isn't converted:

| You drop | Preview in browser | Upload & process on server |
|---|---|---|
| Raw cells (H5AD, Xenium, MERSCOPE) | Parsed in your browser | Converted by a worker |
| Raw molecules (Parquet / CSV) | Parsed in your browser | Converted by a worker |
| Raw cells + transcripts and/or an image | Cells only | All converted or uploaded, and overlaid |
| Already-processed cells, molecules, image | **Opened as they are, overlaid** | — (use Upload & Save from the preview) |
| An OME-Zarr image | **Opened as it is** | Uploaded as it is, beside raw cells |

---

## Ways to load a dataset

On the homepage, pick **Single Cell** or **Single Molecule**, then choose one of:

### 1. Preview in browser
The file is parsed **entirely in your browser** (h5wasm / hyparquet in web workers) and opens instantly — nothing leaves your computer. Best for quick looks and datasets that fit in browser memory (roughly < 500K cells / tens of millions of molecules).

Data that is [already processed](#already-processed-data-opens-as-it-is) — merfisheyes output folders, a SpatialData store, an AnnData table in Zarr, an OME-Zarr image — is not parsed at all: it is read from disk as the view needs it, so it opens at once whatever its size.

### 2. Upload & process on server
For large datasets, upload the **raw bytes** and let a background worker do the processing. Requires signing in (so we know where to email you). The flow:

1. Drop your file/folder → confirm the dataset name.
2. **Single-molecule uploads add a "Confirm columns" step** — a preview of the first rows with dropdowns to map gene / x / y / z / cell-id (auto-filled from the detected platform). This is what guarantees the server reads your columns correctly.
3. A full-screen progress bar tracks the upload; you can leave once it says "Uploaded — we'll email you."
4. An AWS worker runs the same Python processors used for pre-chunking, writes the result to S3, and **emails you a link** to the viewer when it's done.

This path has no browser memory limits and produces the exact same chunked format as everything else. See **[docs/server-side-ingestion/](docs/server-side-ingestion/README.md)** for the full architecture.

### 3. Load from S3
Data already on a public bucket opens straight from it: click **Load from S3** and paste the URL of a merfisheyes folder, a SpatialData store, an OME-Zarr image, or the folder holding a store with its molecules beside it (which opens with everything overlaid). Nothing is copied. The bucket must allow the site in CORS. A dataset opened this way that nobody owns yet can be **added to your account**, and its link reopens it — including the image settings you chose.

---

## Accounts, projects & community

Signing in (**Google**, via NextAuth) unlocks:

- **Your datasets** (`/account`) — everything you've uploaded, with the same cards as Explore, plus edit controls for titles, descriptions, thumbnails, and visibility.
- **Projects** (`/account/projects/[id]`) — group related datasets (e.g. the cells + molecules from one experiment) into a folder with its own page and thumbnail.
- **Community** — submit a dataset or project for review; once an admin approves it, it appears in the **Community** view on **[Explore](https://www.merfisheyes.com/explore)** for everyone.

Explore (`/explore`) lists curated and community datasets with a persistent **Featured** row on top; community projects open a detail page (`/explore/[id]`).

---

## Viewer

Datasets open at `/viewer/[id]` (single cell) or `/sm-viewer/[id]` (single molecule). You can also open local or S3 data directly via `/viewer/from-local`, `/viewer/from-s3`, and the `sm-viewer` equivalents. SpatialData stores open at `/spatialdata-viewer/[id]` or `/spatialdata-viewer/from-s3?url=…` (the store, or the folder holding it), and an image on its own at `/image-viewer/[id]`, `/image-viewer/from-s3?url=…` or `/image-viewer/from-local` — see [docs/spatialdata](docs/spatialdata/README.md).

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

Requires Node.js 20 (`.nvmrc`; the browser test runner breaks on newer majors).

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
│   ├── spatialdata-viewer/    # SpatialData stores (+ [id], from-s3)
│   ├── image-viewer/          # an OME-Zarr image on its own (+ [id], from-local, from-s3)
│   ├── explore/               # curated + community catalog (+ [id], bil/[bilCode])
│   ├── account/               # your datasets & projects (+ edit)
│   ├── admin/                 # admin dashboard & catalog
│   └── auth/signin/
├── components/                # React components (viewers, controls, plots, upload, cards)
├── lib/
│   ├── adapters/              # H5ad / Xenium / Merscope / ChunkedData adapters
│   ├── ingest/                # server-ingestion helpers (classify-folder, error classification, …)
│   ├── spatialdata/           # Zarr v3 stores, tables and OME-Zarr images, read lazily
│   ├── local/                 # already-processed folders: recognise, open, save as they are
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

Pull requests welcome. PRs target the `develop` branch. Read
[docs/PHILOSOPHY.md](docs/PHILOSOPHY.md) first: it states the rules a
change is judged by, and [docs/specs/](docs/specs/) holds the budgets,
contracts and format profiles behind them.

## Acknowledgments

Built at the Bogdan Bintu Lab, UCSD.

## License

MIT
