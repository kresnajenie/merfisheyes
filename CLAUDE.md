# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

This file is a **map**, not a specification. The specifications live in
[`docs/`](docs/) and are the source of truth; where this file and a doc in
`docs/specs/` disagree, the doc wins. Read the doc before changing the thing
it describes.

---

## 1. What this is

**MERFISHEYES** — a web viewer for spatial transcriptomics, live at
[merfisheyes.com](https://www.merfisheyes.com). It draws **cells**,
**molecules**, **cell polygons** and a **microscopy image** in one Three.js
scene, reading data lazily from S3, from a public bucket by URL, or from a
folder on disk — at the scale of an 11.24-billion-molecule store.

Current version: see `package.json` (**0.13.0** at time of writing) and the top
of [`CHANGELOG.md`](CHANGELOG.md).

### Read these first

| Order | Document | Why |
|---|---|---|
| 1 | [docs/PHILOSOPHY.md](docs/PHILOSOPHY.md) | The eight rules a change is judged by. Short. Non-negotiable. |
| 2 | [docs/specs/](docs/specs/) | The format profiles, the layer contract, the loading and rendering rules, the budgets. |
| 3 | [docs/plans/](docs/plans/) | Live, multi-session work briefs with steps and gates. |

The philosophy in one line each: **one implementation per concern**; an element
is a kind and a kind is a contract; **the format is the product**; read the
view, not the dataset; **human time is the budget**; correctness is shown from
outside; start from the raw material; the largest dataset we have is the design
target.

---

## 2. Rules for agents

- **Never commit or push unless explicitly asked.** Never `git tag`, push a
  tag, create a GitHub Release, or open a `develop` → `main` PR unless the user
  asks for it by name. See [docs/RELEASING.md](docs/RELEASING.md) §"Rules for
  agents".
- **Branching**: feature branches → PR into **`develop`** (staging,
  dev.merfisheyes.com). A release is a promotion of `develop` → `main`
  (production, www.merfisheyes.com). Never work directly on `main`.
- **Plan files are living documents.** A step that changes a fact or a decision
  in `docs/plans/*.md` **edits that file in the same commit**. Same for a spec
  in `docs/specs/` when behaviour changes.
- **Measure before optimising.** Perf claims belong in
  [docs/specs/BUDGETS.md](docs/specs/BUDGETS.md) with the run that produced
  them, not in a commit message.
- **If the worker or `scripts/process_*.py` changed**, the ECR image must be
  rebuilt or Fargate keeps running the old code (§8).

---

## 3. Where the truth lives

`docs/` is current. **`scripts/README.md` and `DATA_PROCESSING_ARCHITECTURE.md`
at the repo root are older than the specs** — use them only for the legacy
per-script flag reference, and trust `docs/` on anything they contradict.

### Specifications — `docs/specs/`

| File | Covers |
|---|---|
| [store-layout.md](docs/specs/store-layout.md) | What a folder must look like to be opened: Zarr v3, consolidated metadata, `spatialdata_attrs`, CORS, **coordinate systems (§8 — the transform rule)** |
| [table-profile.md](docs/specs/table-profile.md) | A SpatialData table as cells: `obsm/spatial`, `X` CSC, obs columns |
| [shapes-profile.md](docs/specs/shapes-profile.md) | Cell polygons; **spots and bins drawn at real size** (Visium / Visium HD) |
| [image-profile.md](docs/specs/image-profile.md) | OME-NGFF 0.5 multiscale images, `c,y,x`, placement |
| [mapping.md](docs/specs/mapping.md) + [schemas/](docs/specs/schemas/) | `mapping.json`: what is drawn beside a cell dataset |
| [LAYER-CONTRACT.md](docs/specs/LAYER-CONTRACT.md) | **The `Layer` lifecycle and the element registry. Read before adding any drawable kind.** |
| [LOADING.md](docs/specs/LOADING.md) | `ElementStore`: one cache, cancellation, merged ranges, priority queue, budgets per reader |
| [RENDERING.md](docs/specs/RENDERING.md) | On-demand frames, uniforms vs attributes, picking grid, draw order |
| [BUDGETS.md](docs/specs/BUDGETS.md) | The human-time acceptance criteria and every measured number to date |
| [sddb-profile.md](docs/specs/sddb-profile.md) | The 1,884-store SpatialData-DB census against our profiles |
| [sddb-catalog.md](docs/specs/sddb-catalog.md) | How SpatialData-DB rows become Explore cards |
| [AUDIT-2026-10-17.md](docs/specs/AUDIT-2026-10-17.md) | The findings LOADING/RENDERING answer |
| [SECOND-IMPLEMENTER-2026-10-17.md](docs/specs/SECOND-IMPLEMENTER-2026-10-17.md) | What a fresh agent found missing when adding a kind from the contract alone — the contract's backlog |

### Formats and pipelines

| File | Covers |
|---|---|
| [docs/PIPELINE-SPEC.md](docs/PIPELINE-SPEC.md) | **The chunked S3 layout, byte for byte — and the rule that browser JS ≡ server Python.** Canonical for the legacy format. |
| [docs/PARITY-TESTING.md](docs/PARITY-TESTING.md) | `npm run parity` — runs both pipelines on the same inputs and diffs. 9/9 identical. |
| [docs/points-zarr-spec.md](docs/points-zarr-spec.md) | **Draft 0.1** — the multiscale, doubly-ordered Zarr points element (molecules streamed by view) |
| [docs/molecule-pyramid.md](docs/molecule-pyramid.md) | The research and measurements behind that format |
| [docs/data-requirements.md](docs/data-requirements.md) | Per-format input requirements of the **browser** adapters. ⚠ Its "normalized to [-1,1]" and "sparse X unsupported" sections are stale — see §6. |

### Systems

| File | Covers |
|---|---|
| [docs/server-side-ingestion/README.md](docs/server-side-ingestion/README.md) | **Shipped.** The "Upload & process on server" path end to end: routes, AWS Batch worker, callback, email, AWS resources, **§12 how to deploy a worker change** |
| [docs/spatialdata/README.md](docs/spatialdata/README.md) | Opening SpatialData stores, already-processed folders, images as datasets, overlays, the `mapping.json` surface |
| [docs/plans/ONE-PIPELINE.md](docs/plans/ONE-PIPELINE.md) | **Live plan: one pipeline, run by you or by the server** — the `prepare.py` contract (steps 1–3 done) |
| [docs/plans/BIL-CLUSTER.md](docs/plans/BIL-CLUSTER.md) | **Live plan: the same pipeline on the PSC cluster** — many sections in one store, the 87-dataset fleet, billion-scale molecules |
| [docs/plans/SDDB.md](docs/plans/SDDB.md) | Live plan: MERFISHEYES on SpatialData-DB (steps 1–4 done) |
| [docs/plans/UPLOAD-TUTORIAL.md](docs/plans/UPLOAD-TUTORIAL.md) | Live plan: one way to prepare and upload data (`prepare_*.py`) |
| [docs/labelled-single-molecule.md](docs/labelled-single-molecule.md) | `/lm-viewer` — molecules with N categorical labels, no expression matrix |
| [docs/sc-sm-one-space.md](docs/sc-sm-one-space.md) | SC+SM overlay; `__all__` ships, per-sample mode is future work |
| [docs/SM-STREAMING-PLAN.md](docs/SM-STREAMING-PLAN.md) | Streaming the legacy per-gene `bin.gz` files |
| [TESTING.md](TESTING.md) | The three test layers, the dataset databank, perf baselines |
| [docs/RELEASING.md](docs/RELEASING.md) | SemVer, branch model, the release checklist |

---

## 4. Commands

```bash
npm install
npm run dev                  # dev server, :3000, Turbopack
npm run build                # prisma generate + migrate deploy + next build (4 GB)
npm run build:low-memory     # same at 2 GB
npm run lint                 # eslint --fix

npx prisma generate          # after any schema change
npx prisma migrate dev       # create + apply a migration
npx prisma studio

npm test                     # unit (Vitest), tests/unit/
npm run test:e2e:quick       # @quick Playwright suite on committed tiny fixtures
npm run test:e2e             # full E2E
npm run test:upload          # @upload — needs test:infra:up (Postgres + MinIO)
npm run test:perf            # perf suite → perf/results/latest.json
npm run test:perf:compare    # vs perf/baselines/<platform>.baseline.json; non-zero on regression
npm run parity               # browser JS vs server Python pipeline diff
npm run test:data:make-tiny  # regenerate the committed fixtures
```

Node **20** (`.nvmrc`) — the Playwright runner breaks on newer majors.
Python for the `scripts/` needs `anndata numpy pandas scipy pyarrow`; the
`prepare_*.py` family additionally needs `spatialdata` (+ `spatialdata-io`,
`scanpy`, `leidenalg`, `igraph`) — see `scripts/requirements-prepare.txt`.

---

## 5. Architecture: two generations, both live

This is the single most important thing to understand. **There are two data
paths, and most confusion comes from mixing them up.**

### A. SpatialData / Zarr (current, where new work goes)

A store is read **as written** and never converted. Everything goes through one
contract and one store:

```
ElementStore (lib/storage/)         one byte cache, cancellation, merged ranges,
  └─ Backend                        priority queue · public URL | app bucket | dropped files
       │
  registry (lib/scene/registry.tsx) one ElementDefinition per kind: detect, open,
       │                            createLayer, Bind, Controls
       ├─ ImageLayer    (lib/spatialdata/ImageLayer.ts)
       ├─ PointsLayer   (lib/spatialdata/PointsLayer.ts)   tiled Zarr points
       ├─ ShapesLayer   (lib/spatialdata/ShapesLayer.ts)
       └─ MoleculeLayer (lib/webgl/MoleculeLayer.ts)       legacy per-gene files
       │
  LayerHost (lib/scene/LayerHost.ts)  one rAF loop, one view rect per frame
```

- Every layer implements `Layer` (`lib/scene/contract.ts`): `attach`, `detach`,
  `tick`, `status`, `subscribe`, `setBudget`, `abort`, `dispose`. Extend
  `BaseLayer`. **Never run your own `requestAnimationFrame` loop.**
- `tests/unit/layer-contract.test.ts` is the conformance suite — "implements the
  contract" means it passes there.
- One dataset loader: `lib/viewer/load-dataset.ts`. One overlay resolver:
  `lib/viewer/resolve-overlays.ts` (pure, tested).
- **Adding an element kind is five steps** and touches no viewer — the steps and
  what a second implementer found missing are in
  [LAYER-CONTRACT.md](docs/specs/LAYER-CONTRACT.md).

**The entry point for preparing data is [`scripts/prepare.py`](scripts/prepare.py)**
`INPUT OUT`: it recognises an `.h5ad` / AnnData Zarr, a Xenium export, a
MERSCOPE region, a Visium output, a molecule table or a folder of parts **by
content**, runs the matching `prepare_*` script, and writes exactly
`OUT/<name>.zarr` (one SpatialData store), `mapping.json` and `stats.json`,
never touching INPUT. Exit codes: **0** written and validated, **2**
unrecognised, **3** incomplete, **4** invalid. Packages are pinned in
`scripts/requirements.txt`. Shared rules live in `scripts/meyes_common.py` (the
control-probe filter, molecule column names, obs coordinate columns) —
`worker/entrypoint.py` and `lib/ingest/molecule-columns.ts` repeat the column
names and tests keep the three equal.

The server worker runs `prepare.py` for an upload made with
`NEXT_PUBLIC_INGEST_PIPELINE=prepare` (`processingParams.pipeline`, read through
`lib/ingest/pipeline.ts`; steps 2–3 of
[docs/plans/ONE-PIPELINE.md](docs/plans/ONE-PIPELINE.md)): **one dataset per
upload**, finished as `datasetType: "spatialdata"`, its compute tier chosen up
front from the upload's files (`lib/ingest/compute-tier.ts`). Without the
setting it runs the legacy chunked processors below.

### B. Legacy chunked format (still the output of both processors)

`StandardizedDataset` + adapters, and `SingleMoleculeDataset` for molecules.
This is what `process_spatial_data.py` / `process_single_molecule.py` write,
what the server ingestion worker produces without the flag above, and what
thousands of existing S3 datasets are in. **Specified byte-for-byte in
[docs/PIPELINE-SPEC.md](docs/PIPELINE-SPEC.md).**

```
datasets/{id}/
├── manifest.json              statistics, files, normalized:false
├── coords/spatial.bin.gz      [u32 n][u32 dims][f32…]  raw µm, rounded 2 dp
├── expr/index.json            genes → chunk_id, position_in_chunk
├── expr/chunk_%05d.bin.gz     one gene per chunk; sparse [idx…][val…]
├── obs/metadata.json          per column: categorical | numerical
├── obs/{col}.bin.gz           binary obs v1 (lib/utils/obs-binary.ts)
├── palettes/{col}.json        categorical only, sorted label order
└── de/{col}.bin.gz            DE stats: counts, means, pct

# single molecule
├── manifest.json.gz
└── genes/{sanitized}.bin.gz   f32 [x,y,z,…], 3 floats/molecule even in 2D
    genes/{g}_uuuuuuuuuu.bin.gz  unassigned, only when non-empty
```

**The invariant: the same input must yield the same output from the browser
pipeline and the Python pipeline.** The Python is the reference. `npm run
parity` enforces it across 9 cases; read
[docs/PARITY-TESTING.md](docs/PARITY-TESTING.md) before touching either side —
coordinate rounding, float32 narrowing, obs dictionary JSON and gene filtering
all have exact rules that are easy to break silently.

### Where they meet

`lib/adapters/ChunkedDataAdapter.ts` reads format B; `lib/spatialdata/*` reads
format A; both go through `ElementStore` where they can. A legacy molecule
folder can be converted into a points element with
`scripts/legacy_to_points.py`, after which the `mapping.json` link is not
needed.

---

## 6. Facts that are easy to get wrong

- **Coordinates are raw microns now, not normalized.** `manifest.normalized:
  false`, rounded to 2 dp in float64 (`np.round` semantics, round-half-even)
  then narrowed to float32. Older S3 datasets *are* normalized to [-1,1];
  `StandardizedDataset` reads both (`normalized` flag). **An image can only be
  drawn under cells with raw coordinates.**
- **`chunk_size` is 1** — one gene per expression chunk, both pipelines.
- **Sparse h5ad `X` (CSR/CSC) is supported** by the server pipeline and the
  worker path. The *local browser* `H5adAdapter` still reads only a dense `X` —
  keep local h5ad test fixtures dense (see TESTING.md).
- **A missing categorical label is `""`**, not `"Unknown"`.
- **Column type detection**: `lib/utils/column-type-detection.ts` ≡ Python
  `is_categorical`. Name lists, ≥80 % unique → numerical, >50 % floats →
  numerical. Roughly: ≤100 unique values is categorical in the viewer.
- **Gene filtering** (control probes, blanks, codewords) is single-sourced:
  `lib/utils/gene-filters.ts` ≡ `should_filter_gene`, pinned by
  `tests/unit/gene-filter-parity.test.ts`.
- **An element's `coordinateTransformations` is a list of *alternatives*, one
  per output coordinate system — not a chain.** Multiplying them placed Visium
  images 148× too large. The rule is in
  [store-layout.md §8](docs/specs/store-layout.md). A placement is a 2-D affine
  (`Map2`: scale, offset and, when it rotates, `skew`) applied as **one matrix
  on the layer's group**; points and shapes must still place as written.
- **A table with no `obsm/spatial` still draws** (CosMx): its cells sit at one
  point per shape of the region element it annotates — the Point
  `<region>_for_vitessce` sibling first — matched by instance id
  (`lib/spatialdata/centroids.ts`, `wkb.ts`).
- **The image layer is shared by all three viewers** and takes any channel
  count; the shader has `MAX_VISIBLE_CHANNELS` (8) slots, so an 84-channel
  CODEX image compiles and shows 8 at once.
- **Two zarrita versions are installed**: `zarrita` 0.5.1 (the h5ad-zarr cell
  path only) and `zarrita-v3` (alias for 0.7.5, everything SpatialData). Zarr v3
  string arrays need 0.7.
- **Zarr v2 SpatialData stores are not opened.** Zarr v3 only, with
  consolidated metadata at the root — a public URL cannot be listed, so that
  metadata is the only way the viewer learns what a store holds.
- The viewer shows `Dataset.title`, **not** `manifest.name`.

---

## 7. Ways data gets in

| Path | What happens | Entry points |
|---|---|---|
| **Preview in browser** | Raw H5AD/Xenium/MERSCOPE/parquet/CSV parsed in a web worker. Bounded by browser RAM. | `components/file-upload.tsx`, `lib/workers/` |
| **Already-processed folder** | A SpatialData store, an AnnData zarr table, a merfisheyes output folder, an OME-Zarr image — **opened as it is**, read piecewise from disk, any size. SpatialData parts in separate folders **join the cells' store as its elements** (`joinedStoreFileMap`), in the preview and on save alike; several folders dragged at once are walked as one. | `lib/local/classify-bundle.ts`, `open-bundle.ts`, `save-bundle.ts`, `components/bundle-upload-modal.tsx` |
| **Upload & Save** | The browser-processed or already-processed parts go to S3. Processed data dropped on the server card skips the worker and opens the same save dialog, unopened. | `app/api/datasets/*`, `app/api/single-molecule/*`, `lib/upload/` |
| **Upload & process on server** | Raw bytes to `raw/{id}/`, AWS Batch (Fargate) runs the Python, writes `datasets/{id}/`, HMAC callback, SES email. Sign-in required. `lib/ingest/finalize.ts` writes the `mapping.json` when the parts finish; under `NEXT_PUBLIC_INGEST_PIPELINE=prepare` the same folder becomes **one** dataset instead (cells and molecules always, images only when ticked; a Xenium export or MERSCOPE region sends what `prepare.py` reads — `preparePlatformFiles`). | `app/api/ingest/*`, `worker/`, `lib/ingest/` |
| **Load from S3 / by URL** | A public store, folder, or image opened directly; nothing copied. "Add to my account" registers the URL. | `lib/spatialdata/resolveUrl.ts`, `lib/datasets/register-or-claim.ts` |
| **Prepare it yourself** | `scripts/prepare_folder.py` (+ `prepare_xenium/merscope/visium.py`, `image_to_ome_zarr.py`, `transcripts_to_points.py`) build one SpatialData store, ending in `validate_store.py`. | `scripts/`, [/docs/prepare](https://merfisheyes.com/docs/prepare) |

**There is no server-side cell typing.** The user runs MapMyCells (or anything)
and uploads a one-row-per-cell CSV; the worker merges it as another cluster
column with its own palette and DE stats.

**The validator is the gate**: `python scripts/validate_store.py <store.zarr>
[--mapping]` — run in CI on the tiny fixtures and by hand on real stores. It
applies the same coordinate-system rule the viewer does.

---

## 8. Infrastructure

- **Database**: PostgreSQL + Prisma. Key models: `User`, `Dataset`,
  `CatalogDataset` + `CatalogDatasetEntry` (Explore), `Project`,
  `UploadSession`, `UploadFile`, `DatasetView`. Enums: `DatasetStatus`
  (`UPLOADING → QUEUED → PROCESSING → COMPLETE | FAILED`),
  `DatasetFaultCategory` (`USER | PLATFORM | UNKNOWN`), `FileStatus`.
- **`Dataset.datasetType` is free text.** Use `lib/datasets/kind.ts`
  (`datasetKind`, `viewerBase`, `viewerHref`) — never compare the string
  yourself. Kinds are `single_cell | single_molecule | image`; `"spatialdata"`
  is single cell that opens in its own viewer; old rows hold raw format names.
- **Auth**: NextAuth v5, Google sign-in, Prisma adapter. Admin surfaces under
  `/admin` and `/api/admin`.
- **Storage**: S3 (`merfisheyes-staging`), `raw/{id}/` (7-day lifecycle) and
  `datasets/{id}/`. Presigned URLs; multipart above ~100 MB.
- **Processing**: AWS Batch on Fargate, account 533267148861, us-west-2. Job
  definition `merfisheyes-ingest`, 4 vCPU / 16 GB. Worker image is
  `worker/Dockerfile` + `worker/entrypoint.py` + `scripts/process_*.py`
  **copied in** — ⚠ **editing those files changes nothing until the ECR image is
  rebuilt and pushed** ([§12](docs/server-side-ingestion/README.md#12-deploying-a-worker-code-change)).
- **Email**: AWS SES (`lib/ses.ts`).
- **CI** (`.github/workflows/`): `test.yml` (unit), `e2e.yml` (e2e-quick and
  e2e-upload on every PR, perf nightly).

---

## 9. Code map

```
app/
  viewer/ sm-viewer/ lm-viewer/        cells · molecules · labelled molecules
  spatialdata-viewer/ image-viewer/    each with [id], from-s3, from-local
  explore/ account/ admin/ docs/       catalog, your data, admin, user guide
  api/datasets|single-molecule|ingest|projects|admin|catalog|community

lib/
  scene/        contract.ts, BaseLayer, LayerHost, registry.tsx, useLayers — the layer system
  storage/      ElementStore, backends, fetchRange, FileMapStore, PresignedFetchStore
  spatialdata/  store, resolveUrl, coordinates, footprint, table adapter,
                image/points/shapes readers + their Layers
  viewer/       load-dataset, open-dataset, resolve-overlays
  adapters/     H5ad, H5adZarr, Xenium, Merscope, ChunkedData, ProcessedSingleMolecule
  local/        classify-bundle, open-bundle, save-bundle  (dropped folders)
  molecules/    storeLoading  (legacy per-gene files through ElementStore)
  ingest/       classify-folder, finalize, overlay, error-classification, reconcile, …
  upload/       uploadChunked|Raw|Image|Zarr ToS3
  sddb/         catalog-schema, catalog-file, catalog-row
  datasets/     kind, metadata, owner, read-genes, register-or-claim
  stores/       Zustand; create*Store factories give each split panel its own
  workers/      standardized-dataset, single-molecule, shapes, zarr-densify
  webgl/        scene-manager, shaders, point-cloud, picking-grid, MoleculeLayer
  utils/        obs-binary, gene-filters, color-palette, coordinates,
                column-type-detection, url-state-writer, test-hooks

components/     79 files + docs/ explore/ plots/ admin/ account/ ui/ auth/
worker/         Dockerfile + entrypoint.py (the Batch job)
scripts/        Python processors, prepare_*, validate_store, points_zarr,
                parity/, perf/, testdata/, bil-scripts/ (HPC SLURM)
tests/          unit/ (38 Vitest), python/, fixtures/, data/ (tiny committed)
e2e/            Playwright: tests/viewer, tests/upload, perf/
perf/           baselines/, results/, real/, sddb/ (census)
docs/           specs/, plans/, server-side-ingestion/, spatialdata/
```

### State management

Zustand. Each store has a `create*Store` factory so **split screen** can give
each panel its own instance. Notable: `visualizationStore` (single cell scene:
`mode: VisualizationMode[]` supports combined gene + celltype),
`singleMoleculeVisualizationStore`, `labelledMoleculeVisualizationStore`,
`layersStore` (what the registry published), `localBundleStore` (parts of a
dropped folder), `splitScreenStore`, `uploadStore`, `tourStore`.

Components used to branch on `pathname.startsWith("/sm-viewer")` to pick a
store. New code should take its stores from the panel context / `Bind`
component instead.

---

## 10. Testing

Three layers — full guide in [TESTING.md](TESTING.md).

- **Unit** (Vitest, `tests/unit/`) — pure helpers, readers, the **layer
  contract conformance suite**, the registry, parity of the gene filter.
- **E2E** (Playwright, `e2e/`) — `@quick` on committed tiny fixtures (no
  Docker, runs on every PR), `@upload` against Postgres + MinIO in Docker.
  `e2e/helpers/fixture-server.ts` stands in for a public bucket, answering CORS
  preflights and `Range`.
- **Perf** (`e2e/perf/`, `scripts/perf/`) — fixtures at a capped 12 MB/s plus
  `scripts/perf/real-store.ts` against real stores. Compared to
  `perf/baselines/<platform>.baseline.json`; **headless Chromium draws with
  SwiftShader**, so frame numbers are only comparable with themselves.
- **Parity** (`npm run parity`) — both processing pipelines, diffed.
- App instrumentation: `window.__merfish` (`lib/utils/test-hooks.ts`), only
  when `NEXT_PUBLIC_E2E=1` — layer hooks, frame intervals, authoritative counts.

The test dataset databank is `tests/data/manifest.json`; tiny fixtures are
committed, medium/large are generated or downloaded on demand.

---

## 11. Known gaps and active backlogs

These are recorded where they will be read; listed here so they are not
rediscovered.

- **Layer contract backlog** — no reusable `TileLayer` base (a second tiled
  kind re-derives ~200 lines); two hidden enumerations of kinds; the URL-slot
  writer has six hard-coded slots. [SECOND-IMPLEMENTER-2026-10-17.md](docs/specs/SECOND-IMPLEMENTER-2026-10-17.md)
- **Rendering not done** — cells through a lookup texture; pooled draws for the
  points element; stacked image levels. [RENDERING.md](docs/specs/RENDERING.md)
- **Loading not done** — the h5ad-zarr cell path still uses
  `PresignedFetchStore`, not `ElementStore`. [LOADING.md](docs/specs/LOADING.md)
- **Pipeline Tier-2 divergences** — two manifest key spellings, embedding name
  drift, **fingerprints differ between the browser and server paths so there is
  no cross-path dedupe**. [PIPELINE-SPEC.md](docs/PIPELINE-SPEC.md)
- **Budgets not met** — a one-event zoom still reads ~16.6 MB against a 12 MB
  target; "All" detail cannot be seamless over the internet and is labelled so.
  [BUDGETS.md](docs/specs/BUDGETS.md)
- **Not built** — a QC stage in ingestion; per-user compute tiers; labels as a
  drawable kind; SC+SM per-sample (non-`__all__`) overlay; reading a store's
  `points.parquet` directly.
- `spatialdata` 0.8 cannot read our points element — the spec has not landed
  upstream.

---

## 12. BIL / HPC pipelines

`scripts/bil-scripts/` holds the SLURM batch pipeline for processing BIL
datasets on the PSC cluster (`launch_*.sh` + `*.sbatch`, sample lists in
`*.csv`), with its own [README](scripts/bil-scripts/README.md),
[SETUP](scripts/bil-scripts/SETUP.md) and
[INPUT_REQUIREMENTS](scripts/bil-scripts/INPUT_REQUIREMENTS.md). It runs the
legacy `scripts/process_*.py` on a cluster, then syncs to S3.

**Moving it onto `prepare.py` is planned in
[docs/plans/BIL-CLUSTER.md](docs/plans/BIL-CLUSTER.md).** Read that before
changing anything here. Three facts it records that will otherwise cost you a
day:

- **Three input families**, each with its own sample list and column
  conventions: MERSCOPE multi-slice (`samples.csv`, `csv-samples.csv` —
  `sample_name,species,percentile`), H5AD (`h5ad-samples.csv`, molecules in
  `segmented_spot_table.csv` as `gene_names,x,y,z,cell_ids`) and Xenium
  (`xenium-samples.csv`).
- **Cells are combined across slices; molecules are not.**
  `combine_slices_v3.py` grid-shifts every slice into one frame, tags each cell
  with `_sample_id` and writes the exact per-slice shift to
  `sample_offsets.json`; `process_single_molecule.py` writes one folder per
  slice. Both derive `sample_id` the same way (the slice's top-level folder
  name), so joining them is a lookup, not a search.
- **Pipeline history lives on unmerged branches.** `pipeline/final-rerun`,
  `mmc-update`, `s3-update`, `testing/full-pipeline` and
  `feature/spinal-rna-spots` each carry BIL work not on `develop`, and they
  disagree about `map_my_cell.py` and `process_spatial_data.py`. Check there
  before concluding something was never written;
  `scripts/bil-scripts/find_missing_scripts.sh` searches the cluster for the
  two scripts that are in no branch at all.
