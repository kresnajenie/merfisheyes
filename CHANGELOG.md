# Changelog

All notable changes to MERFISHEYES are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
See [docs/RELEASING.md](docs/RELEASING.md) for how a release is cut.

## [Unreleased]

_Nothing yet._

## [0.4.1] - 2026-09-04

Single-molecule loading becomes visible and cancellable, and pre-chunked
single-molecule uploads stop timing out on production.

### Fixed

- **Pre-chunked single-molecule uploads no longer 500 on production.** The
  initiate endpoint created one DB row and one presigned URL per file serially
  inside a single transaction (~33s for a 529-gene dataset — past Vercel's
  function timeout). It now batch-inserts the file rows and presigns outside
  the transaction in parallel batches (~2s).
- **Streamed default genes render again.** When every default gene streamed,
  the camera auto-fit ran against still-empty clouds and produced a broken
  (NaN) camera that never recovered, leaving the viewer permanently blank.
- **WebGL failures are no longer silent.** If the browser can't create a WebGL
  context (GPU blocklist, half-applied browser update, software GL), both
  viewers now show a "3D view unavailable" panel with remediation hints
  instead of a blank canvas.

### Changed

- **The single-molecule scene shows on first data.** The camera fits as soon
  as the first gene (or first streamed batch) has points, so molecules appear
  progressively instead of after every default gene finishes downloading.
- **Gene loads show progress and truly cancel.** Whole-file downloads report
  byte progress in the loading toast (e.g. `Drd1: 12.4 / 38.2 MB`); loading
  toasts are blue, and closing one aborts the download and deselects the gene.
  Deselecting a gene or leaving the view also aborts in-flight downloads.
- **Streaming kicks in at 100k molecules** (was 2M), so most genes stream in
  with visible progress.
- **Gene expression no longer fades low-expression points by default**
  (`EXPRESSION_ALPHA_MIN` 0.3 → 1.0; still adjustable in the advanced panel).

## [0.4.0] - 2026-09-01

Single-molecule overlays become first-class on single-cell datasets, the DEG
panel is rebuilt around a plain comparison of means, and signing in to upload
stops getting in the way.

### Added

- **Single-molecule overlays on single-cell datasets.** Link a molecule dataset
  to a cell dataset and its spots render on top of the cells. The overlay picker
  is searchable and paginated, works for both app-uploaded and S3-registered
  molecule datasets, and the selected molecule genes now persist in the viewer
  link (`ov=`) so a shared URL restores them.
- **Route DEG gene clicks to the overlay.** When a molecule overlay is present,
  a Cell / Molecule toggle beside the Expression control chooses whether clicking
  a gene colours the cells by expression or adds it as a molecule layer.
- **Expression-mode explainer.** An info button by the Expression control shows
  how counts-vs-log-normalized is auto-detected and how it affects the means and
  fold change.

### Changed

- **DEG panel shows mean1 / mean2 / linear fold change** instead of
  log2FC / mean_in / mean_pct — a direct comparison of means, with counts-vs-log
  auto-detected (and overridable) and log means un-logged before the ratio.
  Higher-contrast target / reference and mean columns.

### Fixed

- Signing in no longer blocks uploads: the email field is hidden when signed in
  and the account email is preserved, so "Process & Save" / "Process & Upload"
  enable correctly. Sign-in happens in place from the upload modals.
- Retroactive overlays no longer orphan a duplicate re-upload; the finalize step
  redirects to the existing dataset.
- Large single-cell `.h5ad` files with BigInt values no longer crash the loader,
  and the homepage featured datasets and server-fallback modal behave correctly.

## [0.3.0] - 2026-08-25

Large datasets stop being a wall: whole-transcriptome Xenium runs convert with
their full expression matrices, out-of-memory jobs retry themselves on bigger
machines, and big single-molecule genes stream into the viewer instead of
loading whole.

### Added

- **Xenium `cell_feature_matrix.h5` support.** Current Xenium exports ship
  expression only as the 10x HDF5 matrix, which neither pipeline could read —
  datasets loaded cells with zero genes. The server now reads it natively
  (a 717,576-cell x 27,104-gene run converts in 44 minutes), and the browser
  loads standard panels in full; whole-transcriptome matrices too large for a
  tab fall back to cells + gene list instead of crashing.
- **Automatic memory-tier escalation.** Processing jobs that run out of memory
  are retried on larger instances (16 -> 32 -> 64 GB) without user action, and
  admins can retry failures from the dashboard.
- **Progressive single-molecule streaming.** Gene files are written in
  bit-reversal order and streamed into the scene, so a half-loaded gene shows
  the whole slide densifying instead of one corner.
- **Account & Projects overhaul**: dataset cards, a project picker, optimistic
  toggles, and a shared dataset detail page.
- **Benchmark harness** (`bench/`, `scripts/bench/`): measures browser vs
  server processing across a real-dataset databank, with per-machine limits,
  costs, and a self-contained report page; 64-dataset results included.

### Changed

- The ingestion worker reads `.h5ad` selectively and streams the expression
  matrix: a 3-million-cell dataset that previously needed >64 GB now converts
  in a 16 GB instance (peak 13.3 GB).
- Browser and server pipelines produce byte-identical chunked output,
  enforced by a parity harness.
- The single-molecule gene picker paginates (250 genes per page) —
  thousand-gene panels no longer lag the panel open.

### Fixed

- Ingestion failures report the container's real reason instead of the
  generic task status.
- Community catalog entries resolve thumbnails from the live dataset.
- Unnamed single-molecule files upload correctly.

## [0.2.0] - 2026-08-19

Explore overhaul: the catalog is fast to search, gene search actually works
(and covers uploads), and datasets open identically everywhere.

### Added

- **Multi-gene search** on Explore: comma/space-separated terms are AND-matched
  with fuzzy per-term matching; cards show which genes matched; Enter creates
  one exact chip per term.
- **Gene search over uploaded datasets**: new `Dataset.genes` column filled
  automatically when processing completes (all three upload paths), a gene
  search box on "Your datasets", and gene copy-through on community
  submission. Backfill script for pre-existing rows
  (`scripts/backfill-dataset-genes.ts` — must be run once against production).
- **BIL → Projects organizer** (`scripts/ingest-bil-projects.ts`): groups the
  auto-registered BIL datasets into properly titled Projects and renames the
  "combined_output"-style rows (must be run once against production).

### Changed

- **Explore API is card-slim**: list payloads dropped from ~700 KB to
  ~17–41 KB (gene arrays never leave the database; featured/filters fetched
  once, not per keystroke). Search inputs are debounced with stale-response
  cancellation.
- **Dataset cards**: no more clipped titles/descriptions (fixed flexbox
  slicing), taller layout with investigator · age · publication-year byline,
  tooltips for full text, gradient fallback for broken thumbnails.

### Fixed

- **Split screen loads datasets like the full viewer**: one shared
  open-dataset pipeline (reset → owner config → URL-state overlay) for all
  viewer pages and both split panels — per-sample alignment, rotation, and
  saved colors now apply in the right panel, stores reset between datasets,
  and the loading screen clears only after coordinates are aligned (no more
  plot-then-jump). Same dataset in both panels renders identically.
- In-viewer catalog modal fetched 50 items under a "Showing 1–20" counter;
  now fetches exactly one page.

## [0.1.0] - 2026-08-18

First tracked release — a baseline snapshot of the production platform. Earlier
history is in the git log; changes from here on are recorded per release.

### Added

- **Single-cell viewer** (H5AD / Xenium / MERSCOPE): gene, celltype, and numerical
  coloring; combined gene + celltype mode; interactive legends and gene/numerical
  scalebar; a floating quantification plot panel (boxplots, histograms, celltype
  barplots) with secondary "group by", Top-N / axis caps, and CSV export.
- **Single-molecule viewer**: parquet/CSV, one point cloud per gene, S3 lazy
  loading, 2D/3D toggle, automatic control-probe / unassigned filtering.
- **Two ingestion paths**: in-browser (web workers) and **server-side** (AWS Batch
  worker) with email delivery. Server single-molecule uploads include a
  column-confirm/remap step; unassigned cell-ids are detected across MERSCOPE,
  Xenium, and CosMx conventions.
- **Combined single-cell + single-molecule upload**: a cells folder containing a
  transcripts file uploads both, auto-creates a project, and overlays molecules on
  cells.
- **Accounts** (Google sign-in), **Your datasets**, **Projects**, and an
  admin-reviewed **Community** catalog on Explore.
- **Split-screen** comparison of two datasets with synced camera and selection.
- **Admin processing dashboard** with failure classification and owner-notified
  failure emails.
- **Python preprocessing** (`process_spatial_data.py`, `process_single_molecule.py`)
  and BIL HPC / SLURM pipelines for very large datasets.

[Unreleased]: https://github.com/kresnajenie/merfisheyes/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/kresnajenie/merfisheyes/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kresnajenie/merfisheyes/releases/tag/v0.1.0
