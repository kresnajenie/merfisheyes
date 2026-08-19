# Changelog

All notable changes to MERFISHEYES are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
See [docs/RELEASING.md](docs/RELEASING.md) for how a release is cut.

## [Unreleased]

_Nothing yet._

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

[Unreleased]: https://github.com/kresnajenie/merfisheyes/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kresnajenie/merfisheyes/releases/tag/v0.1.0
