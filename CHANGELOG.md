# Changelog

All notable changes to MERFISHEYES are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
See [docs/RELEASING.md](docs/RELEASING.md) for how a release is cut.

## [Unreleased]

_Nothing yet._

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

[Unreleased]: https://github.com/kresnajenie/merfisheyes/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/kresnajenie/merfisheyes/releases/tag/v0.1.0
