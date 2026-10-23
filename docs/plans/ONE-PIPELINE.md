# Plan: one pipeline, run by you or by the server

Raw data — an `.h5ad`, a Xenium export, a MERSCOPE region, a transcripts
table — becomes one SpatialData store by **one script**, `scripts/prepare.py`.
Someone preparing their own data runs it on their machine; the server's
ingestion worker runs the same script, with the same pinned packages, on
what was uploaded. Same input, same bytes out. The browser preview of raw
exports is left as it is (decision 4).

This file is the brief for every session working on it: what was measured
on 2026-10-22, the decisions, the contract, the steps with their gates.
Keep it current: a step that changes a fact or a decision edits this file in
the same commit.

## 1. What is true today (measured on develop `ef15bf7`)

**Two pipelines that never meet.**

- The server worker (`worker/entrypoint.py`) runs `scripts/process_spatial_data.py`
  (single cell, line 651) and `scripts/process_single_molecule.py` (molecules,
  line 612). They write the legacy chunked format: `manifest.json` + `expr/`
  + `coords/` + `obs/`, and `manifest.json.gz` + `genes/*.bin.gz`. The
  image (`worker/Dockerfile`, `python:3.12-slim`) holds those two scripts
  only; `worker/requirements.txt` pins `anndata==0.12.10` and has no
  `spatialdata` / `spatialdata-io`.
- A person preparing their own data runs `scripts/prepare_folder.py`,
  `prepare_xenium.py`, `prepare_merscope.py` or `prepare_visium.py`. They
  write one SpatialData store (table with a CSC copy of `X`, images, the
  tiled points element, shapes) and end in `scripts/validate_store.py`.
  `scripts/requirements-prepare.txt` is unpinned: installed today it gives
  `anndata 0.13.4`, `spatialdata 0.8.0`, `spatialdata-io 0.7.1`, `zarr 3.4.1`.
- The new scripts import from the old one: `should_filter_gene` and
  `write_de_stats_binary` come from `process_spatial_data.py`
  (`prepare_merscope.py`, `transcripts_to_points.py`,
  `spatialdata_cluster.py`). The old file cannot simply be deleted.

**What the server receives** (`lib/ingest/classify-folder.ts`): for Xenium
`cells.csv[.gz]` and the expression matrix; for MERSCOPE `cell_metadata.csv`
and `cell_by_gene.csv`; a transcripts file; `.zarr` folders as candidate
images. Not enough for the new scripts:

- spatialdata-io's `xenium()` with every image and label flag off needs 7
  files — `experiment.xenium`, `cell_feature_matrix.h5`, `cells.parquet`,
  `cells.zarr.zip` (needed even with labels off), `cell_boundaries.parquet`,
  `nucleus_boundaries.parquet`, `transcripts.parquet`. Measured on the 10x
  2-FOV exports: 24.7 MB of 284 MB (lung), 28.8 MB of 385 MB (breast), both
  read to the same cells, boundaries and transcripts as the whole export.
  On the full mouse brain export the images are 9.2 GB of 13 GB.
- `prepare_merscope.py` reads `cell_metadata.csv`, `cell_by_gene.csv`, and
  when present `cell_boundaries.parquet` and `cell_categories.csv`; the
  transcripts go straight into the tiled points. spatialdata-io's
  `merscope()` would also need `images/micron_to_mosaic_pixel_transform.csv`
  even with images off — the reason the script builds the store itself.

**The worker's limits.** Three compute tiers (`COMPUTE_TIERS` in
`lib/batch.ts`): standard 4 vCPU / 16 GB, large 8 / 32 GB, xlarge 16 / 64 GB,
100 GiB ephemeral disk. A job that runs out of memory is retried one tier up
(`lib/ingest/mark-failed.ts`), so a 40 GB job fails twice before it runs.

**Measured on this machine, one run at a time** (wall clock, peak RSS;
`scripts/` at develop `ef15bf7`, packages from `requirements-prepare.txt`
installed 2026-10-22). "Old" is what the worker runs today, "new" the
`prepare_*` script that makes one store. Old MERSCOPE / h5ad rows are the
cells only (molecules are a second job); new rows include the molecules.

| Input | Size | Old: cells | Old: molecules | New: one store |
|---|---|---|---|---|
| Xenium lung 2-FOV (10x) | 11,898 cells, 826k tx | — | — | 0:06, 1.3 GB |
| Xenium V1 mouse brain | 13 GB export | 0:39, 0.5 GB | **fails** (bytes gene names) at 5.6 GB | 1:16, 9.0 GB (with images) |
| MERSCOPE Mouse_receptors | 198 MB counts, 3.4 GB tx | 0:24, 0.6 GB | — | 1:10, 1.9 GB |
| MERSCOPE Human_Ovarian_Cancer | 151 MB, 2.4 GB | 0:25, 0.5 GB | — | 0:43, 1.8 GB |
| MERSCOPE hum_brain2 | 266 MB, 5.7 GB | 1:00, 1.4 GB | — | 1:35, 3.2 GB |
| MERSCOPE human_uterine_cancer | 797 MB, 32 GB | 5:50, 6.5 GB | — | 7:58, 10.8 GB |
| MERSCOPE colon (704k cells) | 1.2 GB, 20 GB | 2:41, 8.8 GB | — | 2:51, 12.5 GB |
| BIL `ace-low-cop` h5ad | 785 MB | 0:43, 0.4 GB | — | **fails**: `uns['cell_polygons']` too large for Zarr v3 |
| Beaumont set1 h5ad (124k × 1,067) + set1 molecules | 457 MB + 6.1 GB | 0:57, 1.4 GB | (set1 below) | 17:41, 15.7 GB: 124,364 cells, 851,483,449 molecules after the gene filter |
| Beaumont set2 molecules | 86M rows, 548 MB | — | 2:15, 11.2 GB | 1:24, 3.3 GB |
| Beaumont set4 molecules | 367M rows, 2.3 GB | — | 9:54, **41.3 GB** | 6:09, 8.5 GB |
| Beaumont set3 molecules | 545M rows, 3.4 GB | — | 15:44, **62.4 GB** | 9:03, 12.8 GB |
| Beaumont set1 molecules | 975M rows, 6.1 GB | — | 28:55, **110.2 GB**: no tier fits | 17:16, 15.6 GB |

The old molecule processor needs about 0.11 GB per million molecules: the
32 GB tier from ~140M, the 64 GB tier from ~290M, and nothing past ~580M.
The new one needs about 0.016–0.024 GB per million; set1, at 975M, peaks at
15.6 GB — inside the 16 GB tier with no margin, so the gate in step 2 is on
the large tier for anything over ~500M molecules until the writer streams.

**Deterministic.** `prepare_xenium.py` twice on the lung export: every file
identical except spatialdata-io's parquet copy of the transcripts
(`points/transcripts/points.parquet`), which the viewer does not read.

**What the old pipeline does that the new one does not (yet):**

| Old (`process_spatial_data.py`) | New (`prepare_*`) |
|---|---|
| DE stats for every categorical column (`write_de_stats_binary`, line 2024) | only for `leiden`, with `--cluster` |
| Merges a user's MapMyCells / label CSV (`--mmc-csv`) | no such option |
| Drops control probes and blanks from every format (`should_filter_gene`) | MERSCOPE and molecules only; an `.h5ad` table keeps them, and `lib/spatialdata/` does not filter |
| One output per part (cells, molecules, image: three datasets) | one store; the store also keeps spatialdata-io's parquet points copy (`WARN … not read by the viewer`) |
| — | the store is named after the input folder (`outs.zarr` for a Xenium export) |
| Single molecule upload with no cells | `transcripts_to_points.py` writes a store with only points (opens in the molecule viewer since #165) |

Visium / Visium HD have a script (`prepare_visium.py`) but no server path;
adding one is not in this plan.

**What uploading a lab's own data to staging showed** (Beaumont4GBM, 2026-10-22,
`localhost:3000` on the staging `.env`, the real Batch queue):

- **Molecules with their own column names never process.** Server mode asks
  no columns (`startPreparedUpload` removed the step: a large CSV left the
  dialog blank), and the worker's `detect_molecule_type` knows only
  `global_x` (MERSCOPE) and `x_location` (Xenium). `Gene, X, Y, Z` →
  "Could not detect molecule schema"; the dataset only says "container
  exited with exit code 1" — the worker's reason never reaches it.
- **Molecules beside cells are found by file name.** `scoreTranscriptFile`
  wants "transcript" in the name; `…set1__single_molecules.parquet` beside
  the h5ad was dropped without a word, and only the cells uploaded (124,364
  cells × 1,067 genes, COMPLETE in 476 s).
- **The old molecule processor needs the 64 GB tier from ~290M molecules**
  (table above): three of the four Beaumont sets would fail at least once.
- **Prepared stores upload and open** through the server card's save-as-is
  path (`feat/processed-drops`), once two viewer bugs were fixed there: the
  molecule scene did not frame a points element (camera on an empty
  plane), and a table whose index has a name (`var/feature_name`, as this
  h5ad has) failed with "No gene list found" — the second fixed on develop
  meanwhile by #167 (`var`'s `_index` attribute).

| Uploaded to staging | Size | Upload | Opened (first draw) | Registered |
|---|---|---|---|---|
| set2 molecules store | 1.07 GB, 575 files | 40 s (26.5 MB/s) | 6.8 s | 74,229,043 molecules × 1,067 genes |
| set4 molecules store | 4.59 GB, 1,108 files | 157 s (29.3 MB/s) | 6.7 s | 312,120,318 |
| set3 molecules store | 6.70 GB, 1,537 files | 415 s (16.1 MB/s; machine busy) | 6.5 s | 469,710,146 |
| set1 molecules store | 12.03 GB, 2,503 files | 413 s (29.1 MB/s) | 5.0 s | 851,483,449 |
| set1 cells + molecules store | 12.86 GB, 2,726 files | 609 s (21.1 MB/s) | 13.3 s | 124,364 cells × 1,067 genes, molecules over them |
| raw set2 parquet → worker | 548 MB | 89 s | — | **FAILED**: column detection |
| raw h5ad + set1 parquet → worker | 457 MB sent (parquet not offered) | — | — | cells only, COMPLETE in 476 s |

These belong to the new pipeline's contract: molecules are recognised by
their columns, never their name; a person's own column names are accepted
(`--molecule-columns`, and the upload dialog asks when it cannot tell); a
failure reaches the dataset with the script's own message.

## 2. Decisions (user, 2026-10-22)

1. **One pipeline.** The worker runs the same entry point a user runs, with
   the same pinned packages. Output: one SpatialData store per upload.
2. **Images are optional when uploading, off by default.** Cells and
   molecules are always uploaded. The image files are not sent unless the
   person ticks the box, so they never reach the worker.
3. **Test data** is fetched from public sources: the 10x 2-FOV Xenium
   exports (lung, breast; CC BY 4.0) in
   `/data/merfisheyes-test-data/downloads/`, plus the real exports already
   on disk (Vizgen `Mouse_receptors` and `Human_Ovarian_Cancer`,
   `hum_brain2`, the colon and uterine regions, the Xenium V1 and Prime
   mouse brains, the BIL `.h5ad` files).
4. **The browser preview of raw exports stays as it is.**
5. **The main README says exactly which Python file to run and where it
   is**, and a person who has never seen the repo can follow it.

Second round (user, 2026-10-22, after the Beaumont staging test):

- **a. Images are processed when the input has them**, never otherwise.
  On the server that means only when the person ticked the image box (off
  by default), because unticked image files are never uploaded; locally,
  whatever is in the folder. Same script, same rule.
- **b. DE stats: a separate topic** (§8) — stored inside the table rather
  than an `annotations/` folder. Until that lands, step 1 keeps today's
  behaviour of each script.
- **c. Existing datasets stay as they are**, readable as today; no backfill.
- **d. spatialdata-io's parquet copy of the molecules is left out** of the
  store; the tiled points element carries them.
- **e. Molecule columns are recognised by content.** Xenium and MERSCOPE
  names first, then the obvious ones in any case (`gene`/`feature`/
  `target`, `x`/`y`/`z`, optional `cell`/`cell_id`); when still unsure the
  server path asks, with a column step built from the parquet schema or
  the CSV header (so it cannot come up blank). A failure reaches the
  dataset with the script's own message, not "exit code 1".
- **f. The compute tier is chosen up front** from the upload's size
  (molecules from the parquet footer or the CSV's size, cells from the
  table), not by failing on the smaller ones first.
- **g. An `.h5ad` entry the Zarr writer refuses** (the BIL files'
  `uns['cell_polygons']`) is left out with a warning naming it.
- **h. A molecule file beside the cells is found by content**: familiar
  names first (`transcripts.*`, `detected_transcripts.*`), then any
  `.parquet` / `.csv` whose columns are molecules (decision e), whatever its
  name (`…set1__single_molecules.parquet`).

## 3. The contract: `scripts/prepare.py`

```
python scripts/prepare.py INPUT OUT [--name NAME] [--no-images] [--cluster]
                          [--labels CSV] [--molecule-columns GENE,X,Y[,Z]]
                          [--workers N]
```

- **INPUT** is recognised by content, never by name: an `.h5ad` or AnnData
  Zarr; a Xenium export (`experiment.xenium`); a MERSCOPE region
  (`cell_by_gene.csv` + `cell_metadata.csv`); a Visium output
  (`spatial/` + a feature matrix); a transcripts table (`.parquet` /
  `.csv` with gene and coordinate columns); a folder of parts (what
  `prepare_folder.py` takes). A folder may also hold a transcripts table
  and an OME-Zarr image beside the cells.
- **OUT** gets exactly:
  ```
  OUT/
  ├── NAME.zarr/      tables/table (X + layers/X_csc unless X is CSC, obsm/spatial)
  │                   points/molecules (tiled; when there are transcripts)
  │                   shapes/cell_boundaries (when the input has polygons)
  │                   images/… (when the input has an image; --no-images skips it)
  ├── annotations/    de/<column>.bin.gz for every categorical column
  ├── mapping.json    {"cells": "NAME.zarr"}
  └── stats.json      {"cells": n, "genes": n, "molecules": n, "dims": 2|3,
                       "images": [...], "input": "xenium" | "merscope" | …,
                       "pipeline": "<git sha>", "packages": {…}}
  ```
- **Exit codes**: 0 written and validated; 2 input not recognised; 3 input
  incomplete (names the missing file); 4 `validate_store.py` failed. The
  worker reports the last line of output as the failure reason, as it does
  today (`run_processor`).
- **Packages**: `scripts/requirements.txt`, every package pinned. The worker
  image installs exactly that file; the README tells people to install
  exactly that file.
- **Shared helpers** (`should_filter_gene`, `write_de_stats_binary`, the
  palette) move to `scripts/meyes_common.py`; both the old and the new
  scripts import them from there until the old ones go.

## 4. Steps (each its own PR into develop, each with its gate)

### Step 1: one entry point (`feat/prepare-entry`)

`scripts/prepare.py` (dispatch by content to the existing `prepare_*`
functions), `scripts/meyes_common.py`, images read only when present (spatialdata-io's
`morphology_focus` / `morphology_mip` / `aligned_images` / labels flags
set from what is on disk; `--no-images` to skip them anyway), `--labels CSV` (the
MapMyCells merge, into `obs`), DE stats for every categorical column, the
gene filter on every table, `--name`, `stats.json`, pinned
`scripts/requirements.txt`; molecule tables recognised by their columns
(decision e), `uns` entries the writer refuses left out with a warning
(decision g).

Gate: `tests/python/test_prepare_*.py` green; a new `test_prepare_entry.py`
runs each input kind on a tiny fixture; `prepare.py` twice on the lung and
breast exports gives identical bytes; spatialdata `read_zarr` and anndata
`read_zarr` open every output.

**Done (2026-10-23, `feat/prepare-entry`).** As planned, with: obs
coordinate columns (`center_x/y/z`…) used when a table has no
`obsm/spatial`, as the legacy pipeline and the browser do (the BIL h5ad:
151,391 cells, 3-D); the README section of step 4 brought forward. DE stats
inside the table are their own PR (§8). Gate: Python suite 68 passed; the
Xenium test on the 10x lung export; lung and breast twice → identical;
anndata reads every table. **Found:** spatialdata 0.8's `read_zarr` stops on
the tiled points element (`KeyError: '0.3-zarr'` — our own format); with
`on_bad_files="warn"` it reads the table, images and shapes and skips the
molecules. The docs now say so; whether to also write a standard points
element is an open question for the user.

### Step 2: the worker runs it (`feat/worker-one-pipeline`)

`worker/Dockerfile` copies `scripts/prepare*.py`, the modules they import
and `scripts/requirements.txt`; `worker/entrypoint.py` runs
`prepare.py raw/ out/ --name <dataset id>` (plus `--labels` when a CSV
came with it; images are processed exactly when the upload included them), reads
`stats.json` instead of `manifest.json`, uploads `out/` under
`datasets/{id}/`, and reports `datasetType: "spatialdata"` in the callback.

Gate — **parity**: `scripts/test_datasets.py --pipeline=prepare` runs the
lung and breast exports through `prepare.py` on the host and through the
worker image (`docker run` against the local MinIO), and compares the two
stores file for file: identical. **Memory**: the app picks the tier from
the upload's size (molecule count from the parquet footer or the CSV's
size; cells from the table) instead of learning it from failed retries, and
every real export in §1 peaks under 75% of its tier, measured in the
container with `--memory` set to that tier. Today's numbers already need
the large tier for Beaumont set1 (15.6 GB); a streaming writer for the
points element is the fix if the user wants everything on standard.

**Done with step 3 (2026-10-23, `feat/one-pipeline-server`)**, behind
`NEXT_PUBLIC_INGEST_PIPELINE=prepare` (stamped on each upload's
`processingParams.pipeline`; unset, nothing changes). The image keeps the
legacy processors, so `:latest` can carry both; it was built and run
locally only, not pushed. Parity: the lung and breast exports through the
worker (download from MinIO, `prepare.py`, upload) and through `prepare.py`
run in the same image give identical folders, `stats.json` included (410
and 405 files). Memory, as the container's anonymous memory (its
`memory.peak` counts page cache and fills to the limit), on the tier
`lib/ingest/compute-tier.ts` picks:

| Upload | Tier | Time | Peak |
|---|---|---|---|
| Xenium lung 2-FOV, 7 files, 24.7 MB | standard | 0:10 | 0.8 GB |
| Xenium breast 2-FOV, 7 files, 28.8 MB | standard | 0:11 | 0.8 GB |
| Beaumont set2 molecules, 548 MB | standard | 1:34 | 2.7 GiB |
| Beaumont set1 h5ad + molecules, 6.5 GB | large | 16:11 | 16 GiB (50%) |

Found on the way: the image's build context left out `meyes_common.py`
(fixed in #171) and the schemas `validate_store.py` reads; `prepare.py`'s
output reached the log only at exit (`PYTHONUNBUFFERED`); numba cannot
cache as a non-root user (`NUMBA_CACHE_DIR`). Gate:
`e2e/tests/upload/prepare-upload.spec.ts` (run with the setting, and with
`E2E_WORKER_IMAGE` the worker image processes the upload and its signed
callback finishes it); the rest of the upload suite passes without the
setting. **Not yet run:** the staging pipeline (AWS Batch) — it needs the
image in ECR.

### Step 3: the app sends what the script reads (`feat/upload-one-store`)

- `classify-folder.ts`: per platform, the files `prepare.py` reads — the 7
  Xenium files above (CSV fallbacks when the parquet ones are absent);
  MERSCOPE's two CSVs, `cell_boundaries.parquet`, `cell_categories.csv`;
  the transcripts table, found by its columns (the parquet schema or the
  CSV header), not by "transcript" in its name; image files (`morphology*.ome.tif`,
  `morphology_focus/`, MERSCOPE `images/`, `.zarr` images) only when the
  image box is ticked.
- The name dialog: **cells** and **molecules** listed, always included;
  **"Also upload the images (N GB)"** unticked by default, with the size.
- `initiate` creates one dataset (`datasetType: "spatialdata"`), no linked
  molecule or image datasets; `finalize.ts` registers the store, its
  counts from the callback's stats, the email link from `viewerHref`.
- The single molecule card: a transcripts file becomes a store holding only
  points.

Gate: `e2e/tests/upload/combined-upload.spec.ts` changed to expect one
store with cells and molecules, images absent unless ticked; the lung
export uploaded through the real staging pipeline (AWS Batch, dev bucket)
opens with cells and molecules; the uploaded byte count equals the 7 files.

### Step 4: say where it is (`docs/prepare-entry`)

- **README.md**: a section near the top, "Prepare your data yourself"
  (text in §7), linked from the table of contents.
- `scripts/README.md`: starts with `prepare.py`; the old processors
  described below it as legacy until step 5.
- `/docs/prepare`, `docs/spatialdata/README.md`, CLAUDE.md's "Python
  Preprocessing Scripts" section, the `prepare_xenium.py` /
  `prepare_merscope.py` docstrings (they still say the transcripts go to a
  `molecules/` folder).

Gate: a fresh session, given only the README and a clean Python 3.12 venv,
prepares the breast export and drops the result in the app.

### Step 5: retire the old pipeline (`chore/retire-chunked-pipeline`)

After step 3 has run on production for two weeks without a rollback:
remove `process_spatial_data.py` and `process_single_molecule.py` from the
worker image and the repo, `scripts/parity/` and `npm run parity` (they
compare the browser with the old processor), and `scripts/test_datasets.py`'s
old mode. The viewer keeps reading existing chunked datasets (`/viewer/[id]`,
`ChunkedDataAdapter`); that code is not touched.

## 5. Rollout and rollback

- Step 2 ships a new image tag (`:one-pipeline`) and a new job definition
  revision; `:latest` and rev 3 stay until step 3 is live. Rollback =
  point the app's job submission back at rev 3 (old processor, old
  output), no data migration either way: each dataset records its
  `datasetType` and opens in its own viewer.
- Steps 2 and 3 go out together (the app must send the files the worker
  reads). A dataset uploaded under the old filter and processed by the new
  worker fails with exit 3, naming the missing file.
- Worker code reaches Fargate only by rebuilding and pushing to ECR
  (`docker build --platform linux/amd64 -f worker/Dockerfile .`).

## 6. What shows it is right

- `spatialdata.read_zarr` and `anndata.read_zarr` open every store the
  server writes; `scripts/validate_store.py` passes it.
- The host and the worker image write identical stores from the same input
  (step 2's parity gate).
- Every real export in §1 peaks under 75% of the tier the app sends it to,
  in a container of that size.
- Beaumont4GBM set1 — the h5ad and its 975M molecules, dropped as the lab
  has them — uploads through the server card and opens with cells and
  molecules on staging.
- The staging upload of the lung export opens with its 11,898 cells, 289
  genes (before the gene filter) and its molecules; images only when ticked.
- A fresh session follows the README alone to a prepared, opened dataset.

## 7. The README section (draft)

> ## Prepare your data yourself
>
> Everything the server does to raw data is one Python script in this
> repository, **[`scripts/prepare.py`](scripts/prepare.py)**. Run it on your
> own machine and you get exactly what the server would make, without
> uploading anything first.
>
> ```bash
> git clone https://github.com/kresnajenie/merfisheyes.git
> cd merfisheyes
> python3.12 -m venv .venv && source .venv/bin/activate
> pip install -r scripts/requirements.txt
>
> python scripts/prepare.py path/to/your_data my_sample/
> ```
>
> | You have | `path/to/your_data` is |
> |---|---|
> | An AnnData file | the `.h5ad` (or its folder, with a transcripts table and an image beside it) |
> | A Xenium export | the `outs/` folder |
> | A MERSCOPE region | the folder with `cell_by_gene.csv` and `cell_metadata.csv` |
> | A Visium / Visium HD output | Space Ranger's `outs/` |
> | Molecules only | the transcripts `.parquet` / `.csv` |
>
> You get `my_sample/my_sample.zarr`, one SpatialData store with the cells,
> the molecules and, when your folder has one, the image (`--no-images` to
> leave it out), checked before the script exits. Drag `my_sample/` onto the homepage to look at it, then
> **Upload & Save** to keep it. Options: `python scripts/prepare.py --help`;
> the full guide is [/docs/prepare](https://merfisheyes.com/docs/prepare).

## 8. Separate topic: DE stats inside the table

Today the DEG panel reads per-cluster statistics from an `annotations/`
folder beside the store (`annotations/de/<column>.bin.gz`, written by
`write_de_stats_binary`: for every gene and every category, the mean of the
log-normalized expression and the fraction of cells expressing it, plus the
category sizes). The user wants them inside the table instead, so a store
is one self-contained object.

Proposal, to plan on its own once the pipeline above is settled:

- For each categorical `obs` column, two dense float32 arrays in the
  table's `varm` (genes × categories), `varm["de_mean__<column>"]` and
  `varm["de_pct__<column>"]`, and a record in
  `uns["merfisheyes"]["de"][<column>]` of the category order, the cell count
  per category and how they were computed (normalization, version).
- `varm` and `uns` are ordinary AnnData elements: anndata and spatialdata
  read and write them, and they travel with the table. Size is small
  (Beaumont: 1,067 genes × ~20 leiden clusters × 2 arrays ≈ 170 KB).
- The viewer reads them through `SpatialDataTableAdapter` and falls back to
  `annotations/` for stores written before.

To decide then: every categorical column (as the old pipeline did), or only
the ones a person names; and whether the browser still computes them when a
table has none.

**Decided (user, 2026-10-23): every categorical column, the server path
first. Done on `feat/de-in-table`:** `scripts/de_stats.py` writes, for every
categorical obs column with 2–500 categories (not the index or join keys),
`varm/de_mean__<column>` and `varm/de_pct__<column>` (float32, genes ×
categories) and `uns/merfisheyes_de/<column>` (categories, cell counts,
method) — the record sits at `uns["merfisheyes_de"]`, flatter than proposed.
`normalize_total` + `log1p` are reimplemented in numpy (scanpy's defaults,
equal to 1.2e-7; no scanpy needed). `prepare_folder.py` (so `prepare.py`,
so the worker in step 2) runs it last on every table; `spatialdata_cluster.py`
writes `obs/leiden` and its stats into the table and no `annotations/`
(`--sidecar-only` keeps the old layout). `SpatialDataTableAdapter` offers
the table's columns and reads them, falling back to `annotations/`. The
browser computes nothing for a store (it never holds the whole matrix), so a
store without them has no DEG panel. Checked: the stats against a hand
computation (Python and the viewer), and the DEG panel listing genes for a
dropped store that has them and absent for one that has not (`deg-table.spec.ts`).
