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
- Technologies (website): Visium 912, spatial transcriptomics 552, Xenium
  180, Visium HD 132, Atera 2. Tiers: silver 1,435, bronze 262, which
  leaves 81 stores with no tier the website shows. NonCommercial 1,090,
  commercial use allowed 591, unknown 97.
- **The census (step 1, done) supersedes these counts**:
  [`docs/specs/sddb-profile.md`](../specs/sddb-profile.md),
  [`perf/sddb/CENSUS.md`](../../perf/sddb/CENSUS.md). The instance holds
  1,884 latest SpatialData artifacts: 1,761 converted by SDDB and 123
  older stores at the bucket's top level (65 Visium, 44 Xenium, 12 earlier
  Visium HD versions, 2 tests). Licence is a Lamin label, not in the store.

### The catalog is a published file (measured 2026-10-18)

- The website is fast because it is built from a pre-computed file, not
  because it queries anything. `theislab/spatialdata-db-client` (pip
  `spatialdata-db`, `import sddb`) publishes a GitHub release
  `catalog-<date>` (nine releases 2026-10-05 → 10-18, most of them
  nightly) holding `catalog.parquet` (1.7 MB, 1,774 rows, 63 columns,
  parquet metadata `sddb_generated_at`, `sddb_n_datasets`),
  `collections.parquet` (26 rows), `gene_search.parquet` (240,151
  Ensembl entities, each a 222-byte bitset over the uid order in
  `gene_search_meta.json`, least-significant bit first, guarded by
  `order_hash`), `gene_index.parquet` (264 MB, per-dataset detection, for
  Python only) and `citations.bib`. The column contract is
  `src/sddb/_catalog_schema.py` (required: `uid`, `technology`, `assay`,
  `organism`, `validation_status`, `zarr_url`; facets and range columns
  derived from it). The engine's nightly job on their HPC writes the file;
  the website is rebuilt from it by a Cloudflare deploy hook
  (`.github/workflows/rebuild-site.yml`); the client reads it with a
  conditional GET on the ETag (`_cache.py`). Their design spec (D61) says
  consumers never query the Lamin database. The release assets carry an
  ETag and no CORS header, so a browser cannot read them and a server
  can; the website serves its own copies with `max-age=0,
  must-revalidate`, one 304 round trip per load.
- The catalog against the census: every catalog uid is in the census;
  110 census stores are not in the catalog (49 Visium HD, 48 older
  top-level stores, 9 Visium, 2 Xenium, 2 test folders), so the catalog
  is the published set and the instance holds more. `zarr_url` is always
  `s3://scverse-spatial-eu-central-1/.lamindb/<uid16>.zarr`, the proxy's
  path. Licence is an id with properties (`CC-BY-NC-SA-4.0` 1,090,
  `CC-BY-4.0` 411, `ncbi-geo-policy` 224, `unknown` 49;
  `license_noncommercial`, `license_redistributable`). Fourteen columns
  are all null, among them `thumbnail_url`, `published_at`,
  `sddb_version`, `container_tag`; `n_obs` is null in 81 rows,
  `n_features` in 112. `vitessce_url` is a complete vitessce.io link in
  1,742 rows.

### Elements, one store per technology, against our profiles

| Element | In their stores | Our reader today |
|---|---|---|
| table | `X` is CSR with 250 M-row chunks (unusable per gene); `layers/csc_counts` CSC, sharded (inner chunks 150–250 k), in 1,845 of 2,085 tables. None in any Visium HD `square_002um` (7.5 M bins, 343 M non-zeros), nor in the **default `square_008um` of 82 Visium HD stores** (45 have one). `instance_key`, `region_key`, `region` set; `obsm/spatial` float64 or int64 (654 tables) | reads; finds any CSC layer under `layers/`; a table without one scans the matrix per gene. Reads the **first table by name**, which in 90 Visium HD stores is `square_002um`, not the default |
| image | `c,y,x`, uint8 (H&E, 3 channels) or uint16 (IF, 1–4 channels), 1–5 levels, chunks 256² (Visium HD), 1024² (most, every SDDB Xenium store); 4096² (32 MB a tile) only in the 44 older top-level Xenium stores, 6000² single-chunk hires images in 13 Visium HD stores; Visium stores hold three images (`full_image`, `hires_image`, `lowres_image`) in three coordinate systems, and 77 have no full image | reads; 4096² and 6000² tiles break the 300 ms budget |
| shapes | one `shapes.parquet` each, `ngff:shapes` v0.3: polygons (`cell_boundaries`, `nucleus_boundaries`), circles (`cell_circles`, Visium spots: Point + `radius`), squares (Visium HD bins at 2, 8, 16 µm, up to 28 M); circles in every store, 313 elements over 60,000 rows (35 of them 2 µm bins filed under `visium/`) | polygons only; circles and squares are not drawn |
| points | `ngff:points` v0.2: a parquet **directory** (`points.parquet/part.N.parquet`), not our tiled points format; 193 elements, about 20.5 B molecules, up to 2.8 B in one store | not read: no molecule streaming for the Xenium and Atera stores |
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

An element's `coordinateTransformations` list holds **alternatives**,
one per output system; it is not a sequence to chain.

What our reader does today (measured 2026-10-18 on develop `978a844`,
`openSpatialDataImage` run on the live stores through the proxy):

- It opens **one** image, the first by name (`elementNames(keys,
  "images")[0]`): on Visium the `…_full_image`. The hires and lowres
  images are never opened.
- It places that image by undoing the shapes (else points) element's
  transforms, and `composeTransforms` **multiplies every entry of the
  list**: on the cerebellum store 1 × 0.150015 × 0.0450045 ≈ 0.00675. The
  full image's level-0 pixel comes out **148.1** scene units wide where
  the spots need 1, so the image is drawn 148× too large; on the Visium
  HD ovary store 197.8× (1 × 0.224837 × 0.0224837).
- Atera and Xenium are placed right only because each element declares a
  single transform: Atera's `morphology_focus` gets 0.2125 = 1 / 4.70588.
- The spots are not drawn at all (circles), nor the Visium HD bins
  (squares).

| Store | Image opened | Level-0 pixel today | Right |
|---|---|---|---|
| Visium `fJSOO03VzCTYtD1o` | `…_full_image` | 148.1 | 1 |
| Visium HD `haOzFbJ7Bc9x7NW6` | `…_full_image` | 197.8 | 1 |
| Atera `jCLxdY1tWKmbfWOB` | `morphology_focus` | 0.2125 | 0.2125 |

None of the five downloaded roots uses `affine`; Xenium stores with an
aligned H&E may (the census counts them).

## 2. Decisions

- **Visium coordinate transformations are in scope now.** A store names
  no single system, so the viewer chooses one target: `global` when every
  element reaches it, else the system the table's region (its shapes
  element) declares identity into, which is the system the table's
  coordinates are in. This rule does not depend on a full-resolution
  image existing. Each element takes **the one entry of its list whose
  output is the target** (or a path through the graph to it), never the
  product of the list. An element that cannot reach the target is not
  drawn and is reported in the layer status, not silently skipped.
- **One image per sample**: the full image when the store has one; the
  downscaled Visium images are not drawn beside it (they are the full
  image's own pyramid levels, worse). A store without a full image draws
  its hires image, mapped through its scale into the target.
- **A store with several samples** (several systems, none common) opens
  one sample at a time, chosen in the viewer; no sample is placed in
  another's system. The census found **none** among the 1,884, so this is
  not built until a store needs it.
- **The table is the default one**: `metadata_general["Default Table"]`
  when the store names it, else a table named `table`, else the first by
  name (today it is always the first by name).
- **Visium spots and Visium HD bins are drawn**: circles and squares as
  shapes kinds, instanced, coloured like polygons from the cells' colours.
  This is what "Visium support" means on screen; without it a Visium store
  shows an H&E image and dots. Spots are the data, so they do **not**
  follow the polygons' 1 mm zoom rule: circles and squares are drawn at
  any zoom while the view holds no more than the per-view cap, and above
  the cap the cell dots stand in (as polygons do today). A cheaper first
  cut, if it looks the same, is the existing cell dots sized from
  `spot_diameter_fullres`.
- **Molecules from parquet points come afterwards.** The Xenium stores'
  `points.parquet` directories are not read in this plan; the conversion
  path (`scripts/points_zarr.py`, `legacy_to_points.py`) and a parquet
  reader by row group are a later step, after the maintainers answer
  whether they will write a tiled index.
- **Labels come afterwards** (the `LabelsLayer` from the second-implementer
  worktree, once a `TileLayer` base exists).
- The catalog is synced **from their published release, not from the
  instance** (decided 2026-10-19, replacing "from the hub REST"): a job
  reads `catalog.parquet` from the GitHub release, validates it against
  their schema, and upserts our records keyed on the versioned uid with
  the release tag and `sddb_generated_at` as provenance; a re-run is a
  no-op. The Lamin REST is used for nothing the catalog carries (their
  spec D61: consumers never query the database). The 110 stores in the
  instance but not in the catalog are not listed. Licence, NonCommercial
  and the source are shown on every card and in the viewer.
- **The explore index is pre-computed and served static**, as theirs is:
  the sync writes immutable, content-hashed files (the card rows, the
  facets, one gene bitset per gene over our uid order) and a small
  `latest.json`; the explore page fetches them through Cloudflare with a
  long `max-age` and filters, facets and gene-searches in memory.
  Postgres stays for the internal, community and admin paths. Today the
  page's API answers in 2.1 s and the ISR shell is 180 KB.
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

### Step 1: the census (`chore/sddb-census`, half a day, no product code) — done

Result: [`docs/specs/sddb-profile.md`](../specs/sddb-profile.md). All 1,884
roots read and classified (9.5 minutes, 6 requests in flight, no 429).
For the later steps: 601 stores need the coordinate-system rule, 77 have
no full image, 5 use `affine` (a CytAssist image only), none has several
samples; circles in every store, squares in the Visium HD ones.

- `scripts/sddb_census.py`: list the artifacts from the hub REST (fall
  back to the website's dataset slugs, which carry the uid), fetch each
  root `zarr.json` through the proxy (a browser user agent; 1,778 × ~150 KB),
  and tabulate per store: technology, element kinds and names, image
  dtype / chunk / levels, shapes geometry (polygon, circle, square) and row
  count where cheap, points layout, CSC layer per table, coordinate
  systems and which elements reach which, `metadata_general`, license.
  Output `perf/sddb/census.json` and a Markdown summary committed.
- Count what the later steps depend on: stores with no full-resolution
  image, stores with more than one sample (systems with no common one),
  stores using `affine` or any transform other than `identity`, `scale`,
  `translation`, `sequence`, spots or bins per store (against the 60,000
  per-view cap), and the 81 stores with no silver or bronze tier.
- Send a browser `User-Agent` (Python's default gets 403) and throttle
  (a few requests in flight, back off on 429 or 5xx): the rate limit is
  unknown until §6 question 6 is answered. Cache the roots on disk so a
  re-run fetches nothing.
- Run `scripts/validate_store.py` logic over the same roots (a `--root-json`
  mode that checks what the root alone can show).
- Deliverable: `docs/specs/sddb-profile.md` (ADR): their stores against
  `table-profile.md`, `image-profile.md`, `shapes-profile.md`,
  `store-layout.md`, row by row, with counts; the list of writer-side
  requests for the maintainers (CSC on every table, chunks ≤ 1024²,
  shapes row groups, a tiled points index or our points format,
  `global` or a declared default system), filed as issues at
  `theislab/spatialdata-db-feedback`.

Gate: every one of 1,778 roots classified, none "unknown". (Held: 1,884 of
1,884, the instance's full count.)

### Step 2: coordinate systems (`feat/sddb-coordinate-systems`)

- `lib/spatialdata/coordinates.ts`: read every element's
  `coordinateTransformations` as a list of alternatives keyed by output
  system, build the graph of named systems, choose the target as decided
  above, and give each element its composed transform to the target
  (`scale`, `translation`, `sequence`, `identity`; reject `affine` with a
  reported error until the census says it is needed). Replace the
  multiply-the-whole-list rule (`composeTransforms` over the shapes
  element's list in `image.ts`) and the "first image by name" pick, and
  the per-element transform reads in `points.ts` and `shapes.ts`, with it.
  `composeTransforms` stays for a `sequence`, which is a chain.
- The table's coordinates are the scene; an image in another system is
  mapped through the graph; an image that cannot reach the target is not
  opened and `openElements` reports why.
- Unit tests on the five downloaded roots (already under
  `tests/data/sddb/roots/`): the chosen target per store, each element's
  composed transform, the rejected ones. The full image's level-0 pixel
  is **1** on the cerebellum and ovary stores (148.1 and 197.8 today) and
  stays **0.2125** on Atera. The conformance suite unchanged.
- `@quick`: a tiny Visium-shaped fixture (three images in three systems,
  spots as circles) added to `e2e/helpers/spatialdata.ts`, asserting the
  full image is placed and the downscaled ones are not opened.

Gate: on the cerebellum store (`fJSOO03VzCTYtD1o`) and the ovary store
(`haOzFbJ7Bc9x7NW6`) the H&E image sits under the spots exactly as in
Vitessce (same store, their viewer button).

Result (2026-10-18): built as above, except that `points.ts` and
`shapes.ts` read no transforms to replace — the census shows every points
and shapes element in the cells' frame, so the registry checks that and
leaves out (in `OpenedElements.problems`) one that is not. The table is the
curated default (`defaultTableName`). Unit: the Visium full image's
level-0 pixel is 1 on both stores, Atera 0.2125; `@quick`: the
`spatialdata` library's own Visium-shaped store
(`scripts/testdata/make-tiny-visium.py`), whose test fails on the old
placement (the image ran off the view) and passes now; perf unchanged
(`BUDGETS.md`). On the two real stores, opened through the proxy, the
H&E lies under the spots and the ovary's bins follow the tissue and its
holes. **Against Vitessce** (SpatialData-DB's own config, rendered from
the same stores): their config places the spots and the full image in
the same system the viewer now targets, and the spots and bins sit on the
same parts of the tissue in both. Our view is the vertical mirror of
theirs (image and cells together): the scene draws y upward for every
dataset, so an image's first row is at the bottom — not a placement error,
but every H&E reads upside down next to Vitessce; left for a decision.
Seen on the way, not placement: an H&E image showed as one blue channel — fixed in this step
(`isRgb` in `ImageLayer.ts`; channels sharing a chunk were read as if each
had its own, so only the first could load: 2,718 images in the census) —
and the ovary's default colouring is a positional obs column (stripes),
not fixed.

### Step 3: circles and squares (`feat/sddb-spot-shapes`)

- `lib/spatialdata/shapes.ts` and `lib/workers/shapes.worker.ts`: read
  Point + `radius` (circles) and the square polygons as instanced
  geometry (one instanced mesh, not a polygon each), coloured from the
  cells' colours like polygons, with the same visibility and opacity.
  Drawn at any zoom while the view holds no more than the per-view cap
  (§2), not only below 1 mm. The default table is `square_008um` in 127
  Visium HD stores, but only 45 of those have a CSC layer on it; the
  `square_016um` table has one in 92. Where the default has none, the
  viewer opens the 16 µm table and says so.
- Visium HD at 2 µm (7.5 M bins) overruns the cap even at 1 mm across
  (about 250,000 bins in view against 60,000). It must not truncate
  silently: either the cap is raised for instanced squares (one instance
  is a few bytes, not a triangulated polygon) after measuring it, or the
  status says the bins are hidden until the view is narrower.
- A `shapesStatus` that says the geometry kind; the validator checks the
  geometry column's types.
- Unit tests on a generated circle and square parquet; `@quick` on the
  Visium-shaped fixture.

Gate: `sliderToFrame` and `hoverFrameMs` unchanged within noise on the
medium store; the cerebellum spots draw at the default (whole-section)
view; the ovary 16 µm table draws its bins at 1 mm; the ovary 2 µm table
at 1 mm either draws every bin in view or says why not.

### Step 4: the catalog (`feat/sddb-catalog`) — built 2026-10-21

What was built (spec `docs/specs/sddb-catalog.md`): the vendored schema
(`lib/sddb/catalog-schema.ts`, pinned to `catalog-2026-10-20`, 3,222 rows),
the row → card mapping (`lib/sddb/catalog-row.ts`), the reader
(`lib/sddb/catalog-file.ts`, hyparquet), the sync
(`scripts/sync-sddb-catalog.ts`: newest release with the asset or `--tag`,
cached by tag, every store registered as an admin-owned `Dataset` and named
from the catalog, upsert by `sddbUid`, stale uids unpublished, the uid-set
gate; `--limit=N` for a sample on dev),
the Prisma columns (`source`, `sddbUid`, `license`, `licenseNoncommercial`,
`sourceRelease`, `sourceGeneratedAt`; migration
`20261021100000_catalog_sddb_source` backfills `source` from the flags), the
**SpatialData-DB** category and the licence chip on cards, the info modal
and the page, and the unit tests on fixtures cut from the release
(`tests/unit/sddb-catalog.test.ts`, `scripts/testdata/make-sddb-catalog-fixture.py`).
User decisions (2026-10-19): sync only when told, dev database first, then
production, no schedule; every store listed, mixed into "All".

Deferred from the list below, each its own step: the static explore index
behind a CDN, the gene-query parity test, the `/spatialdata-viewer/sddb/<uid>`
route and Lamin proxy backend, the nightly perf run, the cron. A card opens
the dataset's page (every Explore card does since #162) and the page opens
the store by its public URL; the bucket's CORS allows `www.merfisheyes.com`
only, so the stores open on production until the scverse team adds the dev
origin. Shipped in v0.13.0 (2026-10-21); production synced by hand after
the deploy. The 10-21 release carried no assets for hours after its tag, so the
sync skips releases without the file.

The original step, for the record:

Raw material: the release assets (`catalog-2026-10-16`, the tag this
plan was measured on, and `latest`), and in
`theislab/spatialdata-db-client` the files `src/sddb/_catalog_schema.py`,
`_cache.py`, `genes.py`, `docs/design/SPEC.md` §4 and §6 and the two
workflows under `.github/workflows/`; ours: `app/explore/page.tsx`,
`app/api/explore/route.ts`, `lib/explore/query.ts`,
`prisma/schema.prisma` (`CatalogDataset`, `CatalogDatasetEntry`),
`scripts/ingest-bil-projects.ts`. Not the website, not the hub REST.

- `docs/specs/sddb-catalog.md` (ADR): the release as the contract (the
  columns we read, how a row becomes a card, the licence shown, the 110
  unlisted stores), the sync's idempotence and provenance, the static
  index and its cache contract, how we learn of a new release (§6).
- Prisma: on `CatalogDataset` a `sddbUid @unique`, `source` (`bil`,
  `community`, `sddb`), `license`, `licenseNoncommercial`,
  `sourceRelease`, `sourceGeneratedAt`; one entry per row,
  `datasetType: "spatialdata"`, opened by uid. The existing `isBil` and
  `bilCode` stay as they are.
- `scripts/sync-sddb-catalog.ts`: a conditional GET on the asset (ETag
  kept beside the file), read with `hyparquet`, the required columns and
  dtypes checked against a vendored copy of their schema pinned to the
  release tag (a contract test fails when the upstream file differs),
  upsert by uid, rows whose uid left the catalog unpublished, the run's
  tag and `sddb_generated_at` recorded. Run from a cron route
  (`app/api/cron/sddb-catalog`, Vercel cron, `CRON_SECRET`) and by hand.
- The static explore index: `docs/specs/schemas/explore-index.schema.json`
  and a generator the sync runs (and the build), writing
  `explore/<hash>/cards.json`, `facets.json` and `genes.bin` (one bitset
  per gene over the index's uid order, their layout) plus
  `explore/latest.json` to the app's bucket behind Cloudflare: the hashed
  files immutable (`max-age=31536000`), `latest.json` at `max-age=60`.
  The explore page reads `latest.json`, then the files, and filters,
  facets and gene-searches in memory; `/api/explore` keeps the internal,
  community and admin tabs. Nothing else on the page changes.
- `lib/storage/backends.ts`: `laminBackend(zarrUrl)` mapping the row's
  `zarr_url` to the proxy URL; the quick suite's fixture server gains a
  `/lamin/` route that mimics the proxy (CORS, Range, 404 JSON, 403 on a
  non-browser user agent).
- `/spatialdata-viewer/sddb/<uid>` opens a catalog row through the one
  dataset loader (`kind: "spatialdata"`, a `lamin` source); the licence
  and `card_url` are in the viewer's dataset info.
- Nightly: the sync, then one real store per technology opened through
  the proxy by the perf script, numbers into `docs/specs/BUDGETS.md`.
  Check once that the proxy accepts the user agent the nightly runner
  sends (Node's `fetch`, Playwright's browser); only Python's default is
  known to get 403.
- Tests: unit on the sync against a fixture parquet in the release's
  schema (insert, re-run a no-op, a uid removed, a schema drift); a
  parity test in the Python job (`pip install spatialdata-db`): ten gene
  queries answer the same uid sets from our bitset and from their
  `GeneIndex.where_expressed`; `@quick`: the explore page lists an SDDB
  card from a fixture index and opens it in
  `/spatialdata-viewer/sddb/<uid>` on the fixture server; perf:
  `exploreToInteractive` and `exploreBytes` added to the suite.

Gate: `sddb.Catalog()` loads the file we consumed and its uid set equals
our published set (count and `order_hash`), checked on every sync and in
CI; the five stores of §1 open by uid; the explore page is interactive
within 1 s at 12 MB/s with cache hits on the index files; the licence is
on every card.

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
- Their own client agreeing with our catalog: the same uid set on every
  sync, the same gene answers in the parity test; the explore budget
  held on production through Cloudflare.

## 5. How to brief a session on a step

```
Read docs/plans/SDDB.md, docs/PHILOSOPHY.md and docs/specs/LAYER-CONTRACT.md,
LOADING.md, RENDERING.md, BUDGETS.md and the profile closest to the
work. Do step <N> of the plan: <its title>. Start from the raw material
the plan names (the downloaded roots under tests/data/sddb/, the proxy
URLs, the catalog release assets and the client's schema), not from the
website or the hub REST. Read through the
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
7. May we poll the catalog releases nightly, or will the release workflow
   (`rebuild-site.yml`, a deploy hook) ping a MERFISHEYES endpoint too?
8. Will `thumbnail_url`, `published_at` and `sddb_version` be filled
   (all null in `catalog-2026-10-16`)?
9. Are the 110 stores in the instance but not in the catalog meant to stay
   unlisted (49 Visium HD, 48 older top-level, 9 Visium, 2 Xenium)?
