# Plan: the preview is for looking, the server saves

**The change.** Data dropped on the homepage still opens in the browser right
away. That is the preview, and it stays. **Upload & Save** stops processing
in the browser. It sends the files you dropped to the server, where AWS
Fargate runs `scripts/prepare.py`, the same script you can run yourself, and
writes one SpatialData store (or keeps an OME-Zarr image as it is).

**The outcome.** No new dataset is written in the old merfisheyes formats:
- the chunked cell folder: `manifest.json`, `expr/`, `obs/`, `coords/`
- the molecule folder: `manifest.json.gz`, `genes/`
- `h5ad-zarr`

New data opens only in `/spatialdata-viewer` and `/image-viewer`. The old
formats stay readable, so existing datasets keep working.

**Relation to the one-pipeline plan.** This finishes
[ONE-PIPELINE.md](ONE-PIPELINE.md) and replaces its decision 4 ("the browser
preview of raw exports is left as it is"). Its step 5 (retire the legacy
processors) becomes steps 6–7 here.

**How to use this file.** It is the brief for every session working on this:
- what was measured on 2026-10-24, on develop `55cdbe2` (v0.14.0)
- the decisions, with the ones still open marked
- the contract
- the steps, each with its gate

A step that changes a fact or a decision edits this file in the same commit.

**Where the evidence is.** The research scripts and their outputs are in
`/data/preview-save-research/` on the lab machine:
- `parity/` compares the preview with the saved store
- `measure/` holds the cost measurements
- `misc/` covers fingerprints, sign-in and old cell folders

The scripts re-run as described in each folder.

## 1. What is true today

### 1.1 Every way a new dataset can still be written in an old format

| # | Where | Writes | `datasetType` | Opens in |
|---|---|---|---|---|
| 1 | Homepage, cells, **Preview in browser** (the default mode): `.h5ad`, Xenium, MERSCOPE, then **Upload & Save** (`upload-settings-modal.tsx`, `GeneChunkProcessor`) | chunked cells | `single_cell` | `/viewer` |
| 2 | Same, molecule file (`.parquet` / `.csv`) (`single-molecule-upload-modal.tsx`, `SingleMoleculeProcessor`) | molecule folder | `single_molecule` | `/sm-viewer` |
| 3 | Same, an AnnData **Zarr v2** folder (`uploadZarrToS3`, `/api/datasets/initiate-zarr`) | `data.zarr/` as is | `h5ad-zarr` | `/viewer` |
| 4 | A **pre-chunked** cell folder dropped on the Folder card (`fromLocalChunked`, then Upload & Save) | chunked cells, verbatim | `single_cell` | `/viewer` |
| 5 | The molecule **Chunked Folder** card (browser only: it has no server mode, `app/page.tsx:534`) | molecule folder, verbatim | `single_molecule` | `/sm-viewer` |
| 6 | Processed parts in a folder, on either card, with chunked cells or a molecule folder among them (`lib/local/save-bundle.ts`: `saveChunkedCells`, `saveMolecules`). On the server card this happens **even with the prepare flag on** (`file-upload.tsx:462`) | chunked cells, molecule folder | `single_cell`, `single_molecule` | `/viewer`, `/sm-viewer` |
| 7 | **Upload & process on server** without `NEXT_PUBLIC_INGEST_PIPELINE=prepare`, and admin Retry of such a dataset (`lib/ingest/resubmit.ts` reuses its stored params). The flag is on in dev and, since 2026-10-24, in production, so this now only reaches a Retry of an older dataset | legacy processors on Fargate | `single_cell`, `single_molecule` | `/viewer`, `/sm-viewer` |
| 8 | **Load from S3** / admin catalog import / seed scripts registering an existing old-format folder by URL (`ensureDatasetForS3Url`) | nothing (a row pointing at the folder) | `manifest.type` or `single_molecule` | `/viewer`, `/sm-viewer` |

Already producing new formats: raw data on the server card under the prepare
flag (→ `spatialdata`), a SpatialData store or OME-Zarr image saved as it is
(→ `spatialdata` / `image`), and the image upload (→ `image`).

Row 8 writes no data, it only names old data that exists; it stays (§2, D8).

### 1.2 What Save needs that the preview does not keep

- **The dropped files.** The viewer's Save hands the modals a parsed
  `StandardizedDataset` / `SingleMoleculeDataset` (`navbar-wrapper.tsx:84-101`),
  never the `File`s. Only `localBundleStore` keeps `File`s, and only for
  processed bundles. The server path (`doServerUpload(files)` in
  `components/file-upload.tsx:454`) takes exactly those `File`s, so the preview
  must keep them.
- **Sign-in is already required everywhere.** Every save route returns 401
  signed out (`/api/datasets/initiate:83`, `initiate-zarr:101`,
  `single-molecule/initiate:81`, `requireUser` in `/api/ingest/*`). What
  differs is the experience:
  - The server card stashes the files, opens `auth-modal.tsx` and resumes
    the upload after sign-in (`file-upload.tsx:438`, `:1418`).
  - The three Save dialogs only close the sign-in popup, so the user has to
    press Upload again.

### 1.3 Duplicate detection

- `Dataset.fingerprint` is **not unique any more**. Migration
  `20261007160000_drop_fingerprint_unique` made it `@@index`. Comments in
  `lib/ingest/finalize.ts:57-63`, `:307-323` and `docs/server-side-ingestion`
  §9 still say `@unique`.
- **Server path.** The worker hashes the raw upload (SHA-256 over the files
  sorted by key, `key:size` plus bytes; `worker/entrypoint.py:483-507`)
  before it branches on the pipeline. Finalize fails a non-admin's upload
  whose hash another dataset has (`finalize.ts:67-76`).
- **Bug: the check spans every owner** (`finalize.ts:74` has no `ownerId`).
  - Re-uploading a public file someone else has fails, and the email names
    *their* dataset's id and title.
  - On the legacy combined path it then writes an overlay onto their dataset
    (`:95`) and puts their dataset in your project (`:115`).
- **The fingerprint can be lost.** Reconcile (a lost callback) finalizes
  without it, leaving the placeholder `raw_{id}` for good
  (`lib/ingest/reconcile.ts:152-156`, `:185`).
- **Browser fingerprints** are written but never checked; the dialogs say so.
  - Their check routes have no caller in the app.
  - The molecule one hashes `NaN` counts on lazily loaded datasets, because
    `getCoordinatesByGene` isn't awaited (`lib/utils/fingerprint.ts:126`).

### 1.4 What `prepare.py` cannot take yet

- **A chunked cell folder.** `recognise` exits 2 when given one. Beside a
  molecule folder, the cells are **dropped silently**.
- **`scripts/chunked_to_h5ad.py` converts one to `.h5ad`, then `prepare.py`
  takes it.** All four real folders tried validate. But the conversion loses
  data:
  - **`obs/*.bin.gz` entirely.** It reads only `*.json.gz`, and the current
    writers write `.bin.gz`, so clusters vanish (`tests/data/tiny/chunked/raw`,
    `GSE274572`: 0 of 2 and 0 of 4 columns).
  - Palettes, embeddings (`coords/umap`, `pca`, …) and the old `de/` stats.
  - **Where the cells were.** The browser saved coordinates as
    `(x − mean) / maxAbs` and kept only the scale. So `norm × scale` recovers
    distances, not the offset: an old folder converted this way would not
    line up with an image or molecules in the original frame.
    `legacy_to_points.py:76-78` has the same limit for molecules.
- **Bug: a molecule folder with no cells (alone, or alone in a folder) gives
  an OUT with only `mapping.json` and `stats.json`, and exit 0.**
  - `prepare_folder.py:772-775` links the molecules and the cleanup in
    `prepare.py:268-277` removes them.
  - The worker checks only that `mapping.json` exists (`entrypoint.py:738`,
    `:802`), so such an upload would finish COMPLETE with nothing in it.
  - An image-only folder probably takes the same path (untested).

### 1.5 The preview and the saved store, side by side

Both sides were read with their real code:
- **The preview:** the browser adapters, `StandardizedDataset` and
  `SingleMoleculeDataset.fromParquet`, run under `vite-node`.
- **The saved store:** `prepare.py` writes it; it is read the way
  `/spatialdata-viewer` reads it (`openSpatialDataStoreFromFiles` →
  `SpatialDataTableAdapter`).
- **Independently:** anndata, h5py and pandas read each source directly.

**What agrees.** On `tests/data/tiny/h5ad/tiny.h5ad` and the Xenium lung
2-FOV export the two sides agree on:
- cells: 300; 11,898
- genes, in the same order: 40; 289, the gene-expression features after the
  control-probe filter
- expression, equal to the source
- coordinates, within the preview's 2-decimal rounding
- Xenium molecules: 823,575 after the filter, with per-gene counts equal to
  the source

The control-probe patterns are the same seven regexes on both sides
(`lib/utils/gene-filters.ts`, `scripts/meyes_common.py`). Both sides use the
same column-type function (`lib/utils/column-type-detection.ts`). The CLAUDE.md
notes on "[-1, 1] coordinates" and "≤100 unique values = categorical" are out
of date.

**What differs:**

| Property | Preview | Saved store | Seen on |
|---|---|---|---|
| Xenium clusters | 10x's graph clustering from `analysis.tar.gz`, as "Cluster" | none (spatialdata-io doesn't read `analysis/`) | Xenium lung |
| Default colour column | `selectBestClusterColumn` on the real types | the alphabetically first column, even a numerical one (`SpatialDataTableAdapter.ts:221`, `:814`) | Xenium lung opens on `cell_area` |
| Colours per label | JS `.sort()`; `uns[<col>_colors]` honoured | numeric `localeCompare`; `uns` colours ignored | tiny: `cluster_2` is `#ffeb3b` vs `#2ca02c` |
| Columns offered | drops all-empty columns, coordinate columns, ids | offers them, adds `z_level`, `nucleus_count` | Xenium lung, MERSCOPE |
| MERSCOPE categories | one column from `cell_categories.csv` | every column | MERSCOPE subset |
| Embeddings | `X_*` in `obsm`; MERSCOPE UMAP from `cell_numeric_categories.csv` | every `obsm` key; MERSCOPE UMAP not read | built h5ad, MERSCOPE |
| Bool labels | `True` / `False` | `true` / `false` | built h5ad |
| Molecules' Z | kept (Xenium: 7.5–31.1) | flat unless "Keep the molecules' Z" | Xenium lung |
| MERSCOPE assigned / unassigned | split by `cell_id` | the points element has `x, y, f` only | by code |
| Xenium control feature types | filtered by name **or** feature type | spatialdata-io keeps Gene Expression only | agree on lung; would differ on 5K "Genomic Control" |
| Gene filter on an AnnData Zarr | none | applied | by code |

**Preview bugs found on the way.** They matter whatever is decided here:
- **An `.h5ad` written by current anndata (0.13, pandas 3: nullable-string
  `var` index) shows 0 genes in the preview.** Its string obs columns are
  skipped too. The store reads it.
- **Cells with a missing coordinate are dropped from the positions but not
  from the rest, so positions and cells shift out of step after the first
  one** (`H5adAdapter.ts:333`, `XeniumAdapter.ts:298-318`,
  `MerscopeAdapter.ts:377-397`). MERSCOPE subset: 2,997 positions for 3,000
  cells. The store keeps all 3,000.
- **A `var` index named in its attributes (e.g. `gene_ids`) gives 0 genes in
  the preview** (`H5adAdapter.ts:413` ignores `_index`).
- **MERSCOPE's unnamed id column is offered as a numerical colour column.**

### 1.6 What saving costs, old way and new

Both pipelines ran in the worker image (`--cpus 4`, the standard tier).

| Dataset | New: sent | Old: sent (cells; + molecules) | `prepare.py` (wall, peak) | Store |
|---|---|---|---|---|
| tiny `.h5ad` (300 × 40) | 88 KB, 1 file | 23 KB, 51 files | 3 s, 0.2 GB | 98 KB |
| Xenium lung 2-FOV (11,898 cells, 824k mol.) | 24.7 MB, 7 files | 0.9 MB; 8.2 MB | 9 s, 0.8 GB | 19.4 MB |
| Xenium breast 2-FOV (7,275 cells, 1.1M mol.) | 28.8 MB, 7 files | 1.1 MB; 10.9 MB | 10 s, 0.8 GB | 23.1 MB |
| MERSCOPE Human_Ovarian_Cancer (71k cells, 32M mol.) | 2.71 GB, 3 files | 16.9 MB; 256 MB | 1:09, 1.8 GB | 448 MB |
| Beaumont set1 `.h5ad` + molecules (124k cells, 851M mol.) | 6.53 GB, 2 files; large tier | 338 MB; ~6.4–8.6 GB (est.) | 17:41, 15.7 GB (from ONE-PIPELINE) | 12.9 GB |

**Upload size.** The new path sends more bytes, because the molecules always
go and raw tables are bigger than per-gene gzipped float32: 3–160× the old
cells-only save, 2.6–10.6× the old cells and molecules.

**What the old path could never save.** The browser cannot hold more than
**1.98 GB in one buffer**, nor more than 3.5 GB of JS heap
(`bench/results/linux-x64/limits.json`). Measured:
- h5ad: 3 and 4 GB fail to read
- MERSCOPE regions: the tab crashed from about 280 MB of counts (4.7 GB
  RSS)
- molecule CSVs: from about 8M molecules, time out
- molecule parquet: 77M molecules crash

**Fargate overhead.** From 58 runs:
- queue median 35 s (p90 43 s)
- container overhead (image pull, S3 in and out) median 44 s (p90 64 s),
  measured with the 647 MB legacy image. The `prepare.py` image is 2.06 GB;
  its pull time is not measured yet.

Uploads run at 20.3 MB/s with six parallel parts (4.5 MB/s with one).

**Budget per save, from these** (allow about 1.5× the local `prepare.py`
time on Fargate):

| Size sent | Example | Expect |
|---|---|---|
| < 100 MB | Xenium 2-FOV, small `.h5ad` | about 2 min, at most 3 |
| 0.1–1 GB | `.h5ad`, MERSCOPE cells | 3–5 min (8 on a slow link) |
| 1–5 GB | MERSCOPE Ovarian 2.7 GB | about 6 min; 15–20 on a slow link |
| 5–35 GB | Beaumont set1 6.5 GB; uterine 35 GB | 25–30 min for 6.5 GB; about 45 min for 35 GB |

### 1.7 Production — done on 2026-10-24, ahead of this plan

Most of step 5 happened while this plan was being written.

- **The worker image.** The v0.14.0 worker (`e015155e…`, pipeline `dcfdc3c`)
  is now `merfisheyes-worker:latest`; the one before it is kept as
  `:pre-v0.14.0`. Before promoting it, its **legacy** path was diffed
  against the old image on the tiny `.h5ad`, Xenium lung cells, Xenium
  transcripts and a lab table with a column mapping: byte-identical except
  the molecule manifest's timestamp and temporary id. A lab table *without*
  a mapping fails on the old image and works on the new one.
- **The prepare pipeline is on in production.** Checked from outside by
  reading the live bundle: the `preparePipelineEnabled()` guards are folded
  out and the "skipping N files" string is gone.
- **No separate production job definition is needed.** The app sends its own
  `AWS_S3_BUCKET` with every job (`lib/batch.ts:204`), which overrides the
  definition's. Production submits `merfisheyes-ingest`, dev submits
  `merfisheyes-ingest-staging`, and **both share the one queue** and its
  64-vCPU pool, so a large dev job can delay a production one.
- **Verified, no longer assumptions:** the task role reads and writes
  `merfisheyes-production`; the bucket's CORS exposes `ETag`; a
  `expire-raw-ingest` lifecycle rule was added (there was none); the
  production app user can now also track and cancel jobs and clean up
  abandoned multipart uploads (`MerfisheyesBatchTrack`).
- **It ran.** `sm_4T9Mpl2nwO`: Beaumont set1 molecules, a 6 GB table →
  851,483,449 molecules × 1,067 genes, `datasetType: "spatialdata"`,
  COMPLETE in 28 minutes, `raw/` deleted on success, a real fingerprint
  stored, opening in `/spatialdata-viewer`. It carries **no** assigned /
  unassigned split (D5a).
- **Still unverified:** the 100 GiB ephemeral disk against a 30 GB input,
  and the 64 vCPU ceiling, which allows at most four xlarge jobs at once.

### 1.8 Which saves send an email

Only three paths email at all:
- server ingest, from `finalize.ts` through `lib/ses.ts`
- the browser's Upload & Save of cells or molecules, which asks
  `/api/send-email` or `/api/send-email-single-molecule`
- the AnnData Zarr upload, the same way (`uploadZarrToS3.ts:241`)

**An as-is upload sends nothing**: an image, a SpatialData store or a
processed folder goes through `uploadZarrFolderToS3` and
`/api/datasets/[id]/complete`, neither of which mails anyone. This is what
the user hit on production on 2026-10-24 and what D6 fixes.

`/api/send-email` takes an address from the browser and mails it, with no
check that the caller owns the dataset or that the address is theirs.

## 2. Decisions

**D1. Save goes through the server.** *Decided by the user, 2026-10-24.*
- The preview stays in the browser and only reads.
- Upload & Save sends the dropped raw files (only those `prepare.py` reads,
  `preparePlatformFiles`), and the server runs `prepare.py`.
- People who want to keep data local, or have more than the browser holds,
  run `prepare.py` themselves and drop the store, which uploads as it is.

**D2. Processed data keeps uploading as it is.** *Recommended.*
- SpatialData stores and OME-Zarr images go up without processing (save
  bundle, `initiate-zarr`).
- `prepare.py` would only copy them, and would have to learn image-only
  input first (§1.4).

**D3. Old-format folders dropped for saving: warn now, convert next month.**
*Decided by the user, 2026-10-25.*
- **Now:** a chunked cell folder or a legacy molecule folder still saves as
  it is, as today, with a notice on the Save dialog that the format is
  being retired and that `scripts/prepare.py` turns it into a store.
- **Next month:** `prepare.py` learns both (step 2) and they are converted
  instead. A converted dataset says that where its cells sit is recovered
  to scale but not to offset (§1.4), so an overlay onto it may not line up.
- Nothing is refused in between.

**D4. No duplicate checking at all.** *Decided by the user, 2026-10-25.*
- **Turn the finalize check off** (`lib/ingest/finalize.ts:67-76`): no
  upload is failed for matching another dataset, no "Duplicate of …" email,
  and the overlay/project redirect that hangs off it goes with it. That
  also ends the leak of other owners' ids and titles (§1.3), so step 0 drops
  the owner-scoping fix.
- **The fingerprint is still computed and stored** by the worker — it is a
  useful record and `/api/ingest/mine` shows it. It simply blocks nothing.
- **No early client-side check** is built. Re-uploading the same file makes
  a second dataset, which is the accepted cost.
- The browser fingerprints and the two unused `check-duplicate` routes go in
  step 7 with the rest of the dead code.

**D5. What a preview and its saved store must agree on.** *Decided by the
user, 2026-10-25.*
- **Must match:**
  - the cells, in order
  - the genes, in order
  - expression
  - positions, within 0.005 after the store's placement
  - the colour columns offered, with their types
  - the default colour column
  - each label's colour
  - embedding names
  - molecules: genes and per-gene counts
  - **molecules: assigned and unassigned, drawn apart as they are today**
- **Fixed in `prepare.py` or the store reader, not accepted** (§1.5):
  - **Xenium's own clustering read out of `analysis.tar.gz`** into the
    table — the user asked for this by name
  - **MERSCOPE keeps its UMAP** (`cell_numeric_categories.csv`) — likewise
  - the default column chosen as the preview chooses it
  - the label sort and `uns` colours honoured
  - all-empty and coordinate columns dropped
  - the MERSCOPE `cell_categories` rule
  - bool labels
- **Accepted, and said in the Save dialog:**
  - molecules flat unless "Keep the molecules' Z" (default off, decided
    2026-10-23). The user confirmed this one is fine as it is.

**D5a. Assigned / unassigned in the points format.** *Decided 2026-10-25;
the spec already allows it, so this is a writer and a reader change.*
- **The spec is fine as written** (`docs/points-zarr-spec.md`): a points
  element MAY carry any further column (§4), names them in
  `points:copied_columns` (§6) and may declare `spatialdata_attrs.instance_key`
  when a column assigns points to instances. No spec change.
- **`legacy_to_points.py` already writes it**: the `<gene>_uuuuuuuuuu.bin.gz`
  files become a `u1` column `unassigned` (`scripts/legacy_to_points.py:106`,
  `:133`).
- **`transcripts_to_points.py` does not** (`:214-220` writes x, y, f and z
  only), so every store the new pipeline has written so far — including
  production's `sm_4T9Mpl2nwO` — has no split.
- **`PointsLayer` does not read it either**: "unassigned" appears nowhere in
  `lib/spatialdata/`, only in the legacy `lib/webgl/MoleculeLayer.ts`.
- **Decided:** a `u1` column `unassigned`, not the full cell id. It matches
  what the viewer draws, and it compresses to almost nothing because it is
  nearly all one value; a `cell_id` enum would cost about 4 bytes per
  molecule (≈3.4 GB on Beaumont set1). The spec leaves room to add
  `cell_id` + `instance_key` later for per-cell queries.
- **Source columns:** Xenium `cell_id` (`"UNASSIGNED"`), MERSCOPE `cell_id`
  (`-1`). The molecule column mapping already carries `cellId`
  (`lib/config/moleculeColumnMappings.ts`), so the upload dialog can ask for
  it when it is not recognised.
- **The viewer's controls already exist**: the SM viz store has
  `unassignedColor`, `unassignedLocalScale` and `unassignedShape`, which
  `MoleculeLayer` uses. `PointsLayer` has to honour the same three.
- Stores already written keep no split; re-running `prepare.py` adds it.

**D6. Waiting for a save, and an email for every save.** *Decided by the
user, 2026-10-25.*
- Save shows the upload bar, then "Processing on the server (about N
  min)", with N from §1.6 by size.
- The tab can close; the email and *Your datasets* say when it is ready.
- The preview stays open meanwhile, and "View dataset" replaces it once the
  store exists.
- **An as-is upload emails too.** An image, a SpatialData store or a
  processed folder never goes near Fargate and sends nothing today (§1.8),
  which is what the user hit on production. The email moves to the server
  step every such upload already ends in, addressed to the dataset's owner,
  one message per save however many parts it held. `/api/send-email` and
  `/api/send-email-single-molecule`, which mail any address anyone names,
  are retired with it (step 6).

**D7. Kill switch.** *Decided by the user, 2026-10-25: yes.*
- `NEXT_PUBLIC_PREVIEW_SAVE=server` turns the new Save on for one release.
- Unset means today's dialogs.
- It goes once production has run clean for two weeks.

**D8. Registering old folders by URL** (row 8). *Decided by the user,
2026-10-25: leave it alone.*
- It writes nothing, and points at data that exists.
- Converting someone's public folder is a separate question.

## 3. The contract

| Dropped | Preview (unchanged) | Upload & Save sends | Server runs | `datasetType` | Opens in |
|---|---|---|---|---|---|
| `.h5ad` (+ a molecule table beside it) | `/viewer/from-local` | the `.h5ad` (+ the table) | `prepare.py` | `spatialdata` | `/spatialdata-viewer/{id}` |
| AnnData Zarr v2 / v3 | `/viewer/from-local` | the folder | `prepare.py` | `spatialdata` | same |
| Xenium export | `/viewer/from-local` | `experiment.xenium`, `cell_feature_matrix.h5`, `cells.parquet`, `cells.zarr.zip`, both boundaries, `transcripts.parquet` (+ `morphology*` when ticked) | `prepare.py` | `spatialdata` | same |
| MERSCOPE region | `/viewer/from-local` | `cell_by_gene.csv`, `cell_metadata.csv`, `cell_categories.csv` and `cell_boundaries.parquet` when present, `detected_transcripts.*` (+ images when ticked) | `prepare.py` | `spatialdata` | same |
| Molecule table alone | `/sm-viewer/from-local` | the table (+ the column mapping when asked) | `prepare.py` → points only | `spatialdata` | same (molecule scene) |
| Chunked cell folder (D3a) | `/viewer/from-local` | the folder | `prepare.py` (chunked input, step 2) | `spatialdata` | same |
| Molecule folder alone (D3a) | `/sm-viewer/from-local` | the folder | `prepare.py` (folder alone, step 2) | `spatialdata` | same |
| SpatialData store(s), with or without separate parts | `/…/from-local` | as is (unchanged) | — | `spatialdata` | same |
| OME-Zarr image alone | `/image-viewer/from-local` | as is (unchanged) | — | `image` | `/image-viewer/{id}` |

Every row requires sign-in. The signed-out flow is the server card's: stash
the files, sign in, resume.

## 4. Steps

Each step is its own PR into develop, with its gate.

### Step 0: fix what is wrong today (`fix/no-dedup`, `fix/preview-parse`)

These are worth doing whatever is decided:
- **Duplicate checking off** (D4, `finalize.ts:67-163`): no upload is failed
  for matching another dataset, no "Duplicate of …" email, and the
  overlay/project redirect that hung off it goes too. The fingerprint is
  still computed and stored. This also ends the leak of another owner's id
  and title. The stale `@unique` comments and docs are fixed in the same
  commit.
- **Preview parse fixes:**
  - nullable-string `var` index and obs columns (anndata 0.13)
  - `var` `_index` attribute honoured
  - cells with a missing coordinate kept in step (NaN position, not a
    dropped row)
  - MERSCOPE id column not offered
- **`prepare.py` fails loudly when it writes no store** (exit 3, not an empty
  OUT). The worker fails a job whose `mapping.json` names nothing that
  exists.

**Gate:**
- Unit tests for each fix, on fixtures built from the cases in §1.5:
  `parity/make_edge_h5ad.py`, `make_named_index.py`, `make_merscope_subset.py`.
- A `tests/python` case: a molecule folder alone exits non-zero (until step
  2 makes it work).

### Step 1: one probe for "what the user sees" (`test/preview-store-parity`)

**Built 2026-10-25.** It ended up simpler than first planned: one reader in
TypeScript instead of two scripts and a Python comparator.
- **`tests/unit/helpers/seen.ts`** describes what a person sees, through one
  `describe()` for both sides:
  - `previewOf`: the worker's own parse, as a dropped file opens
  - `storeOf`: the store, opened as the SpatialData viewer opens it
  So a difference is always in the data, never in how it was read. The
  opening column is `selectBestClusterColumn` on the freshly opened dataset,
  the same call `lib/viewer/open-dataset.ts` makes for both.
- **`tests/unit/preview-store-parity.test.ts`** compares eight properties
  per source. Positions compare *whether* each cell is drawn (both
  coordinates finite) and where: within the preview's two-decimal rounding
  plus float32 precision (a CSV's 781.585 is 781.58 in the preview and
  781.5850219… in the store).
- **The stores are committed** (`scripts/testdata/make-tiny-parity-stores.py`
  → `tests/data/tiny/parity/`), so CI's unit job needs no Python.
  `tests/python/test_parity_stores.py` reruns `prepare.py` and fails, with
  the regenerate command, when they no longer match by content.
- **Sources:** the tiny `.h5ad`, the three step 0 `.h5ad`s, the tiny MERSCOPE
  region; a real Xenium export when `MERFISHEYES_XENIUM` and
  `MERFISHEYES_PREPARE_PYTHON` are set. Chunked folders join with step 2.

**Measured — what agrees and what does not.** On the tiny sources 35 checks
agree. The known differences are expected failures in the test's `KNOWN`
list, each with its cause:

| Source | Differs | Cause |
|---|---|---|
| tiny `.h5ad`, MERSCOPE | label colours | the store sorts labels numerically (`cluster_2` before `cluster_10`), the preview as text, so colours are dealt in another order |
| nullable | label colours | the store labels booleans `true`/`false`, the preview `True`/`False` |
| missing_xy | colour columns | the store offers `center_x` / `center_y` |
| MERSCOPE | colour columns | the store offers every `cell_categories` column (adds `celltype`) and the coordinates |
| Xenium lung (opt-in) | colour columns, opening column, colours | no `Cluster` in the store; it adds `z_level` and `nucleus_count`; it opens on `cell_area`; `transcript_counts` labels sort differently |

On the Xenium lung export, cells (11,898), genes (289), expression,
positions and embeddings all agree.

**Gate, met:**
- The suite runs in CI (40 checks).
- Every known difference is an `it.fails` naming step 3. A fixed one fails
  the suite until its entry is removed.

### Step 2: `prepare.py` takes every input the preview does (`feat/prepare-legacy-inputs`)

- **Chunked cell folder as an input:**
  - `chunked_to_h5ad.py` becomes `read_chunked(folder) → AnnData`
  - obs from `.bin.gz` (the layout in `lib/utils/obs-binary.ts`) as well as
    `.json.gz`
  - `coords/*` → `obsm/X_<name>`
  - palettes → `uns/<col>_colors`
  - `× spatial_scaling_factor` when normalised, with the unknown offset
    written to `stats.json` and logged
- **The parts around it:**
  - `kind_of` → `"chunked"`, counted as cells
  - a molecule folder alone → a points-only store via `legacy_to_points`
  - an image alone → a store holding the image (or refused, per D2)

**Gate:**
- `tests/python/test_specs.py` cases on `tests/data/tiny/chunked/{tiny,raw}`
  and a molecule folder alone.
- Each store passes `validate_store.py`.
- obs columns are kept 2/2.
- The step 1 probes agree on cells, genes and expression with the folder as
  the preview reads it.

### Step 3: the saved store looks like the preview (`fix/store-matches-preview`)

The D5 fixes:
- in `prepare_xenium.py`: graph clustering from `analysis.tar.gz` into the
  table
- in `prepare_merscope.py`: the `cell_categories` rule and the UMAP from
  `cell_numeric_categories.csv`
- in `SpatialDataTableAdapter.ts`:
  - the default column by real types
  - the label sort and `uns` colours
  - empty, coordinate and id columns not offered
  - bool labels

And D5a, assigned / unassigned end to end:
- `transcripts_to_points.py` reads the mapping's `cellId` column (Xenium
  `cell_id` = `"UNASSIGNED"`, MERSCOPE `cell_id` = `-1`) and writes a `u1`
  column `unassigned`, named in `points:copied_columns` — the same column
  `legacy_to_points.py` already writes, so both converters agree
- `PointsLayer` reads that column and draws the two kinds apart, honouring
  the `unassignedColor` / `unassignedLocalScale` / `unassignedShape` the SM
  viz store already has for `MoleculeLayer`
- the upload dialog asks for the cell column when the table's columns are
  not recognised, as it already asks for gene, x and y

**Gate:**
- The step 1 suite has no expected failures left, except the one D5
  accept (molecules flat unless Keep Z).
- The store and the preview are compared with the same probes, on every
  fixture.
- A `tests/python` case: a transcripts table with a cell column gives a
  points element whose `unassigned` counts equal the source's, and
  `validate_store.py` passes it.
- `tests/unit`: `PointsLayer` draws assigned and unassigned as separate
  clouds, as `molecule-layer.test.ts` already checks for the legacy layer.

### Step 4: Save sends the files (`feat/preview-save-server`)

- **The preview keeps its source.** `file-upload.tsx` keeps the `{key, File}`
  list it would upload, in a store that lives as long as the dataset.
- **One upload function.** `doServerUpload` moves out of the homepage
  component into `lib/upload/`. The homepage and every `from-local` viewer's
  Save call it.
- **The Save dialog.** It names what is sent and the expected wait (D6),
  with the images and Keep Z ticks. Signed out, it uses stash, sign in,
  resume. A chunked cell folder or legacy molecule folder saves as it does
  today, with D3's notice that the format is being retired.
- **An email for every save** (D6), sent server-side when an as-is upload
  completes, so an image, a store or a processed folder is no longer
  silent.
- **Behind `NEXT_PUBLIC_PREVIEW_SAVE=server` (D7).**

**Gate:** an `@upload` matrix, one case per row of §3, run with
`E2E_WORKER_IMAGE` and `NEXT_PUBLIC_INGEST_PIPELINE=prepare`. Each case:
drop, preview, Save, the worker runs, then:
- `datasetType` is `spatialdata` or `image`
- the dataset opens in the viewer §3 names
- `validate_store.py` exits 0
- **the store equals a local `prepare.py` run on the same files, file for
  file**
- the owner is emailed, including for an image and for a store saved as it is

The matrix extends `prepare-upload.spec.ts` and `drop-matrix-save.spec.ts`.

### Step 5: production (`docs` + console work, no code)

**Mostly done on 2026-10-24, before the rest of the plan** (§1.7): the
worker image, the pipeline flag, the task role, the bucket's CORS and
lifecycle, the app user's job permissions, and a real 851-million-molecule
run through it.

What is left for this step:
- **`NEXT_PUBLIC_PREVIEW_SAVE=server`** on production once step 4 has run on
  dev, then redeploy.
- **Measure the §1.6 budget on dev first, then production:** Xenium lung 2-FOV
  and Beaumont set1, Save from the preview, timed from click to "ready".
- **Watch the shared queue.** Dev and production submit to the same one; if
  a dev job ever delays a production save, production gets its own queue.

**Gate:**
- Both saves land within the budget.
- Both stores open.
- One admin Retry of a legacy dataset still works.

### Step 6: stop the old writers (`chore/deprecate-legacy-writers`)

Two weeks after step 5 with no failures:
- **Retire the old routes:** `/api/datasets/initiate`,
  `/api/single-molecule/initiate`, both check-duplicate routes,
  `/api/send-email` and `/api/send-email-single-molecule` answer **410 Gone**
  with a pointer to Save.
- **Keep the completion routes** (`files/…/complete`, `[id]/complete`) until
  open upload sessions have expired.
- **Close the old card:** the molecule Chunked Folder card goes; its folder
  is accepted by the Folder card (D3a).
- **Retire the kill switch.**

**Gate:** a test that every save entry point (§3) ends as `spatialdata` or
`image`, and that the retired routes answer 410.

### Step 7: delete (`chore/delete-legacy-writers`)

**Writers to delete:**
- `GeneChunkProcessor`, `SingleMoleculeProcessor`, `lib/utils/fingerprint.ts`,
  `sc-upload-files.ts`
- `uploadChunkedToS3`
- `upload-settings-modal.tsx`, `single-molecule-upload-modal.tsx`
- the `isPreChunked` branches
- `saveChunkedCells` and `saveMolecules` in `save-bundle.ts`
- the obs and DE **encoders**
- `scripts/parity/` for the old pipeline
- the legacy worker processors and their finalize and reconcile branches
  (ONE-PIPELINE step 5)

**Readers to keep:**
- `ChunkedDataAdapter`, the obs and DE **decoders**
- `StandardizedDataset.fromS3` / `fromLocalChunked`,
  `SingleMoleculeDataset.fromS3` / `fromCustomS3`, `storeLoading.ts`
- the dataset and single-molecule read routes
- the local preview
- `uploadZarrFolderToS3` with `initiate-zarr` (images and stores, D2)
- `chunked_to_h5ad.py`, `legacy_to_points.py`

**Docs to update:** CLAUDE.md, `docs/server-side-ingestion`, `app/docs`,
README.

**Gate:**
- `tsc`, the unit tests and every e2e suite pass.
- Old datasets still open: the `@quick` viewer specs on the chunked and
  molecule-folder fixtures.

## 5. Rollout and rollback

**Rolling out.** Dev first, with both flags on. Production after step 5's
gate.

**Rolling back:**
- **The Save path.** Unset `NEXT_PUBLIC_PREVIEW_SAVE` and redeploy; the old
  dialogs come back. This works until step 6.
- **Datasets in flight.** Each keeps the pipeline it was uploaded under:
  `processingParams.pipeline` is what the worker, finalize and reconcile go
  by.
- **The worker.** Re-point `:latest` at `:pre-v0.14.0`.

## 6. What shows it is right

None of these is the author's say-so:
1. **The upstream library reads it.** `spatialdata.read_zarr` opens every
   saved store, and `scripts/validate_store.py` exits 0 on it.
2. **Two ways in, one store out.** Save from the preview and a local
   `python scripts/prepare.py` on the same files give the same store, file
   for file. This was shown before for the worker against a local run, on
   the Xenium lung and breast exports.
3. **An independent reader agrees.** anndata, h5py and pandas reading the
   source agree with the saved table on cells, genes, expression and obs
   columns.
4. **Every way in is covered.** The `@upload` matrix fails if any entry point
   in §3 produces an old format, or if a retired route accepts a new upload.
5. **A budget held on real data.** Xenium lung 2-FOV and Beaumont set1 save
   from the preview within §1.6's budget on dev, then on production.
6. **What you previewed is what you saved.** The step 1 probes find no
   difference between the preview and the store on the D5 properties, apart
   from the one accepted one (molecules flat unless Keep Z).
7. **The molecules keep their assignment.** A points element written from a
   transcripts table with a cell column reports the same assigned and
   unassigned counts as the source table, and the viewer draws the two
   apart, as the legacy layer does today (D5a).
