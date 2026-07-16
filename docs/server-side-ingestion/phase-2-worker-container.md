# Phase 2 — Worker container (runs manually first)

← [Back to design doc](README.md)

## Goal

Build the Docker image that runs the Python processor against a dataset's raw S3
prefix and produces the exact chunked layout the viewer already reads. Prove it
works **run-by-hand** (local `docker run`, or against a real `QUEUED` dataset's
`raw/{id}/`) before wiring any AWS triggers. This phase can be built in parallel
with [Phase 1](phase-1-raw-upload.md).

## Depends on

Conceptually independent, but end-to-end testing wants a `QUEUED` dataset with raw
files in S3 (from Phase 1).

## Tasks

### Image

- [ ] `Dockerfile` (new, e.g. `worker/Dockerfile`): Python 3.11 base +
      `anndata, numpy, pandas, scipy, pyarrow` + copy `scripts/process_spatial_data.py`
      and `scripts/process_single_molecule.py` + the entrypoint.
- [ ] Keep the image lean — **no MapMyCells reference data baked in** (that's a later
      stage, mounted from S3/EFS).

### Entrypoint wrapper (thin Python or shell)

- [ ] Read job env: `DATASET_ID`, `DATASET_KIND` (`single_cell` | `single_molecule`),
      `PROCESSING_PARAMS` (JSON), `CALLBACK_URL`, `CALLBACK_SECRET`, bucket/prefix.
- [ ] Download `raw/{id}/` from S3 to local scratch.
- [ ] Run the **enabled stages** (README §5.8) — v1 = the single chunk stage:
  - single-cell → `process_spatial_data.py raw/ out/`
  - single-molecule → `process_single_molecule.py raw_file out/ [--dataset-type …]`
- [ ] Upload `out/` to `datasets/{id}/` (same layout as `ChunkedDataAdapter` /
      `SingleMoleculeDataset.fromS3()` expect).
- [ ] Compute the canonical fingerprint (dedup — README §6).
- [ ] Delete `raw/{id}/`.
- [ ] (Callback POSTs are stubbed/printed here; wired for real in Phase 3.)

## Files touched

- `worker/Dockerfile`, `worker/entrypoint.py` (or `.sh`) — new
- Reuses `scripts/process_*.py` verbatim (no changes expected)

## Verification

- [ ] `docker run` with a local raw folder produces `out/` with
      `manifest.json`, `coords/`, `expr/`, `obs/`, `palettes/` (SC) or
      `manifest.json.gz` + `genes/` (SM).
- [ ] Upload that `out/` to a test `datasets/{id}/` and open it in `/viewer/{id}`
      or `/sm-viewer/{id}` — renders correctly end-to-end.
- [ ] Coordinates match today's pre-chunked datasets (raw, rounded 2dp, `normalized:false`).
- [ ] Image size is reasonable; build is reproducible.

## Notes / risks

- Ephemeral scratch must hold raw + output; size the Fargate task storage accordingly
  in Phase 3 (up to 200 GB configurable).
- Pin dependency versions in the image for reproducibility.
- The processors already emit staged logs (`STEP 4.5`, `chunk 3/10`) — the entrypoint
  will map these to progress callbacks in Phase 3.
