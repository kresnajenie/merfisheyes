# Plan: MERFISHEYES on SpatialData-DB

Open the 1,778 datasets of [SpatialData-DB](https://spatialdata-db-web.pages.dev/)
in MERFISHEYES, streamed from their storage, with the same loading,
rendering and format rules as everything else
([docs/PHILOSOPHY.md](../PHILOSOPHY.md), [docs/specs/](../specs/)).

This file is the brief for every session working on it: the facts
measured on 2026-10-18, the decisions taken, the steps with their gates,
and the prompt to give a session for each step. Keep it current: a step
that changes a fact or a decision edits this file in the same commit.

## 1. What the database is (measured, not read off the website)

- A public Lamin instance `scverse/spatialdata-db`: id
  `f95ff493-5f8b-4c1a-ad80-b781031e8faf`, API
  `https://aws.eu-central-1.lamin.ai/api`, lamindb 2.8.1, storage
  `s3://scverse-spatial-eu-central-1` (region eu-central-1). The settings
  endpoint `GET /api/instances/scverse/spatialdata-db/settings` answers
  without a token. The hub's REST is documented only as Swagger
  (`https://aws.us-east-1.lamin.ai/_docs`, `openapi.json`); the record
  endpoints are `POST /instances/{id}/modules/{module}/{model}/…`
  (`get-values`, `group-by`, `relatives`, `query_features`, `{id_or_uid}`),
  `GET /instances/{id}/artifacts/by-path`, `POST /storages/credentials`.
  The Python client reads anonymously: `ln.DB("scverse/spatialdata-db")`.
- The website is static (pre-rendered from the instance); it has no API.
  Each dataset page carries the artifact uid and a Vitessce config naming
  the elements it uses.
- Every artifact is `.lamindb/<uid>.zarr` in the bucket. Bytes are served
  by the storage proxy
  `https://lamin.ai/storage/s3/scverse-spatial-eu-central-1//.lamindb/<uid>.zarr/<key>`:
  CORS `*` with `Range` allowed and `Content-Range` exposed, byte ranges
  honoured (206), a missing key is 404 (JSON body), and Python's default
  user agent gets 403 (curl and browsers do not). The bucket itself has no
  CORS, so the proxy is the only path from a browser.
- Every store read: Zarr v3, spatialdata 0.8.0 (`spatialdata_attrs.version
  0.2`), consolidated metadata at the root (70 to 220 KB), curated metadata
  in the root attributes: `metadata_general` (Assay, Organism, Tissue,
  Disease, Dataset Url, Publication Date, Default Table, Preservation
  Method, H&E / IF available, …), `sddb_provenance` (`sddb_id`,
  `canonical_key`, `converted_at`), `spatialdata_io_reader`, and a
  per-technology block (`metadata_visium`, `metadata_xenium`, …).
- Technologies: Visium 912, spatial transcriptomics 552, Xenium 180,
  Visium HD 132, Atera 2. Tiers: silver 1,435, bronze 262. NonCommercial
  1,090, commercial use allowed 591, unknown 97.

### Elements, one store per technology, against our profiles

| Element | In their stores | Our reader today |
|---|---|---|
| table | `X` is CSR with 250 M-row chunks (unusable per gene); `layers/csc_counts` CSC, sharded (inner chunks 150–250 k), on the **default table only**. Visium HD 2 µm and 8 µm tables (7.5 M bins, 343 M non-zeros) have no CSC layer. `instance_key`, `region_key`, `region` set; `obsm/spatial` float64 or int64 | reads; finds any CSC layer under `layers/`; a table without one scans the matrix per gene |
| image | `c,y,x`, uint8 (H&E, 3 channels) or uint16 (IF, 1–4 channels), 1–5 levels, chunks 256² (Visium HD), 1024² (most), 4096² (Xenium: 32 MB a tile); Visium stores hold three images (`full_image`, `hires_image`, `lowres_image`) in three coordinate systems | reads; 4096² tiles break the 300 ms budget |
| shapes | one `shapes.parquet` each, `ngff:shapes` v0.3: polygons (`cell_boundaries`, `nucleus_boundaries`), circles (`cell_circles`, Visium spots: Point + `radius`), squares (Visium HD bins at 2, 8, 16 µm, up to 7.5 M) | polygons only; circles and squares are not drawn |
| points | `ngff:points` v0.2: a parquet **directory** (`points.parquet/part.N.parquet`), not our tiled points format | not read: no molecule streaming for the 180 Xenium stores |
| labels | `cell_labels`, `nucleus_labels`: uint32, `y,x`, 4096² chunks, 5 levels | not read (a `LabelsLayer` exists uncommitted in the second-implementer worktree) |

### Visium coordinate systems (the first thing to support)

A Visium store has no `global`. Its elements declare transforms into
**named** coordinate systems, one per image:

```
images/<sample>_full_image     cyx -> <sample>                   identity
images/<sample>_hires_image    cyx -> <sample>_downscaled_hires  identity
                               cyx -> <sample>                   scale [1, 6.666, 6.666]
images/<sample>_lowres_image   cyx -> <sample>_downscaled_lowres identity
                               cyx -> <sample>                   scale [1, 22.22, 22.22]
shapes/<sample>                xy  -> <sample>                   identity
                               xy  -> <sample>_downscaled_hires  scale [0.150015, 0.150015]
                               xy  -> <sample>_downscaled_lowres scale [0.0450045, 0.0450045]
```

(Visium HD is the same shape with `scale [1, 4.4477, …]` and `0.224837`;
the spinal-cord ST store uses `global` with identity everywhere; Atera and
Xenium use `global` with the points and shapes scaled by 4.70588 and the
image identity.) The table's `obsm/spatial` is in `<sample>` (full-image
pixels); `uns/spatial/<sample>/scalefactors` holds `spot_diameter_fullres`.

What our reader does today: `lib/spatialdata/image.ts` takes the **first**
transform of the shapes (else points) element and undoes it to place every
image, assuming one target system. On a Visium store that places the
full image right (identity) and the two downscaled images wrong (they are
in other systems), and the spots are not drawn at all (circles).

## 2. Decisions

- **Visium coordinate transformations are in scope now.** A store names
  no single system, so the viewer chooses one target: `global` when every
  element reaches it, else the system the full-resolution image declares
  identity into (the one the table's coordinates are in). Every element
  is composed into that target; an element that cannot reach it is not
  drawn and is reported in the layer status, not silently skipped. The
  downscaled Visium images are not drawn (they are the full image's own
  pyramid levels, worse).
- **Visium spots and Visium HD bins are drawn**: circles and squares as
  shapes kinds, instanced, coloured like polygons from the cells' colours.
  This is what "Visium support" means on screen; without it a Visium store
  shows an H&E image and dots.
- **Molecules from parquet points come afterwards.** The Xenium stores'
  `points.parquet` directories are not read in this plan; the conversion
  path (`scripts/points_zarr.py`, `legacy_to_points.py`) and a parquet
  reader by row group are a later step, after the maintainers answer
  whether they will write a tiled index.
- **Labels come afterwards** (the `LabelsLayer` from the second-implementer
  worktree, once a `TileLayer` base exists).
- The catalog is synced **from the instance, not pasted**: a job reads the
  artifacts through the hub REST (or `lamindb` in a Python job), keyed on
  uid and version, and writes our dataset records with license,
  NonCommercial, provenance and the curated metadata; the viewer shows
  license and source on every one.
- Their storage proxy is a **backend of the ElementStore** with its own
  contract (CORS, ranges, 404, no listing, user agent), tested nightly on
  one real store per technology.
- The maintainers get a **written profile**, not requests in chat: what
  passes, what we ask them to write differently.

## 3. Steps

Branches off `develop`, in order, one PR each, every PR with unit, `@quick`
and perf gates green, numbers recorded in `docs/specs/BUDGETS.md`, this
file updated. Commit only on the user's go-ahead; never push or open PRs
unless asked.

### Step 1: the census (`chore/sddb-census`, half a day, no product code)

- `scripts/sddb_census.py`: list the artifacts from the hub REST (fall
  back to the website's dataset slugs, which carry the uid), fetch each
  root `zarr.json` through the proxy (a browser user agent; 1,778 × ~150 KB),
  and tabulate per store: technology, element kinds and names, image
  dtype / chunk / levels, shapes geometry (polygon, circle, square) and row
  count where cheap, points layout, CSC layer per table, coordinate
  systems and which elements reach which, `metadata_general`, license.
  Output `perf/sddb/census.json` and a Markdown summary committed.
- Run `scripts/validate_store.py` logic over the same roots (a `--root-json`
  mode that checks what the root alone can show).
- Deliverable: `docs/specs/sddb-profile.md` (ADR): their stores against
  `table-profile.md`, `image-profile.md`, `shapes-profile.md`,
  `store-layout.md`, row by row, with counts; the list of writer-side
  requests for the maintainers (CSC on every table, chunks ≤ 1024²,
  shapes row groups, a tiled points index or our points format,
  `global` or a declared default system), filed as issues at
  `theislab/spatialdata-db-feedback`.

Gate: every one of 1,778 roots classified, none "unknown".

### Step 2: coordinate systems (`feat/sddb-coordinate-systems`)

- `lib/spatialdata/coordinates.ts`: read every element's
  `coordinateTransformations`, build the graph of named systems, choose
  the target as decided above, and give each element its composed
  transform to the target (`scale`, `translation`, `sequence`, `identity`;
  reject `affine` with a reported error until needed). Replace the
  "first transform of the shapes element" rule in `image.ts` and the
  per-element transform reads in `points.ts` and `shapes.ts` with it.
- The table's coordinates are the scene; an image in another system is
  mapped through the graph; an image that cannot reach the target is not
  opened and `openElements` reports why.
- Unit tests on the five downloaded roots (commit them under
  `tests/data/sddb/roots/`): the chosen target per store, each element's
  composed transform, the rejected ones. The conformance suite unchanged.
- `@quick`: a tiny Visium-shaped fixture (three images in three systems,
  spots as circles) added to `e2e/helpers/spatialdata.ts`, asserting the
  full image is placed and the downscaled ones are not opened.

Gate: on the cerebellum store (`fJSOO03VzCTYtD1o`) and the ovary store
(`haOzFbJ7Bc9x7NW6`) the H&E image sits under the spots exactly as in
Vitessce (same store, their viewer button).

### Step 3: circles and squares (`feat/sddb-spot-shapes`)

- `lib/spatialdata/shapes.ts` and `lib/workers/shapes.worker.ts`: read
  Point + `radius` (circles) and the square polygons as instanced
  geometry (one instanced mesh, not a polygon each), coloured from the
  cells' colours like polygons, with the same visibility and opacity.
  Visium HD at 2 µm (7.5 M bins) is bounded by the existing 60,000-per-view
  cap and the 1 mm zoom rule; the 8 and 16 µm tables are the ones with a
  CSC layer and are the default.
- A `shapesStatus` that says the geometry kind; the validator checks the
  geometry column's types.
- Unit tests on a generated circle and square parquet; `@quick` on the
  Visium-shaped fixture.

Gate: `sliderToFrame` and `hoverFrameMs` unchanged within noise on the
medium store; the ovary 16 µm table draws its bins at 1 mm.

### Step 4: the Lamin backend and the catalog (`feat/sddb-catalog`)

- `lib/storage/backends.ts`: `laminBackend(uid)` resolving
  `lamin://scverse/spatialdata-db/artifact/<uid>` through the settings
  endpoint and the storage root to the proxy URL; the quick suite's
  fixture server gains a `/lamin/` route that mimics the proxy (CORS,
  Range, 404 JSON, 403 on a non-browser user agent).
- `/spatialdata-viewer/sddb/<uid>` opens an artifact by uid through the
  one dataset loader (`kind: "spatialdata"`, a `lamin` source).
- The catalog sync: a route or job that reads the artifact records
  (uid, version, key, size, `metadata_general` fields, license,
  NonCommercial) and upserts our dataset records with provenance; the
  explore page shows them with license and source; re-run is idempotent.
- Nightly: one real store per technology opened through the proxy by the
  perf script, numbers into `docs/specs/BUDGETS.md`.

Gate: the five stores of §1 open by uid; the explore page lists the
synced catalog with the same counts as the website.

### Later, out of this plan

- Parquet points (Xenium transcripts) by row group, or conversion with
  `legacy_to_points.py` / `points_zarr.py` when the maintainers accept a
  tiled index.
- Labels (`LabelsLayer`), after a `TileLayer` base is extracted from
  `ImageLayer` (backlog in `docs/specs/LAYER-CONTRACT.md`).
- Visium HD 2 µm and 8 µm gene expression without a CSC layer.

## 4. What shows it is right

- The census: every root classified, the validator's verdict per store.
- spatialdata 0.8 (the writer) and our reader agreeing on each store's
  inventory and counts; Vitessce showing the same store the same way.
- The budgets held on one real store per technology through the proxy:
  first content ≤ 300 ms, a settled Auto view ≤ 1 s at 12 MB/s.
- The conformance and registry suites for the new kinds; the
  second-implementer backlog not growing.
- The maintainers' written answer to `docs/specs/sddb-profile.md`.

## 5. How to brief a session on a step

```
Read docs/plans/SDDB.md, docs/PHILOSOPHY.md and docs/specs/LAYER-CONTRACT.md,
LOADING.md, RENDERING.md, BUDGETS.md and the profile closest to the
work. Do step <N> of the plan: <its title>. Start from the raw material
the plan names (the downloaded roots under tests/data/sddb/, the proxy
URLs, the hub's openapi.json), not from the website. Read through the
ElementStore; never fetch, cache or loop in a layer; constants are
uniforms; nothing per frame while still. It is done when the step's gate
in the plan holds, unit + @quick + perf are green, BUDGETS.md and this
plan are updated. Commit only when told; never push or open a PR unless
told. List every place you had to touch outside the files the step names.
```

## 6. Open questions for the maintainers (ask in the profile, not in chat)

1. Will every table get a `csc_counts` layer, including Visium HD 2 µm and
   8 µm? Without it per-gene reads scan 343 M non-zeros.
2. Can image chunks be at most 1024² (the Xenium 4096² uint16 tiles are
   32 MB each)?
3. Will points get a spatial index (our tiled points format, or parquet
   row groups by tile with a manifest), so a view can be read without the
   whole element?
4. Can shapes parquet files carry row groups ordered by position?
5. Will every store declare one default coordinate system (`global`), or
   is "the system the full-resolution image is identity into" the rule we
   should rely on?
6. Is the storage proxy's 403 for non-browser user agents intended, and is
   there a rate limit a nightly check should respect?
