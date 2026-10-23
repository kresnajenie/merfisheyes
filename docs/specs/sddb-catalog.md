# The SpatialData-DB catalog on Explore

Status: decided 2026-10-21 (`docs/plans/SDDB.md`, step 4). Measured on
release `catalog-2026-10-20` (3,222 rows).

## Decision

Explore lists every dataset of [SpatialData-DB](https://lamin.ai/scverse/spatialdata-db)
from the catalog file the project publishes, not from Lamin's database.
`theislab/spatialdata-db-client` attaches `catalog.parquet` to a GitHub
release `catalog-<date>` (near-nightly; the assets land a few hours after
the tag). Their design spec (D61) says consumers never query Lamin; the
release is the contract, and `src/sddb/_catalog_schema.py` in that
repository is its schema. We vendor the column names we read and the
required set in `lib/sddb/catalog-schema.ts`, pinned to a release tag.

## How a row becomes a card

`lib/sddb/catalog-row.ts` (pure, tested) maps one row to one
`CatalogDataset` with `source = "sddb"` and one entry:

| Card | Catalog column |
|---|---|
| title | `title` (else `sample_id`, else `uid`) |
| species / tissue / disease | `organism` / `tissue` / `disease`, capitalised; organisms mapped to the spellings Explore already filters on (human → Human) |
| platform | `technology`, mapped to the existing spellings (Merscope → MERSCOPE, VisiumHD → Visium HD, SpatialTranscriptomics → Spatial Transcriptomics) |
| description | assay, organism, tissue, disease, `n_obs` × `n_features`, tier |
| cells / genes | `n_obs` / `n_features` |
| licence | `license`, `license_noncommercial` (shown on the card, the info modal and the page) |
| links | `dataset_url` (else `card_url`) as the source link; `doi` as the publication |
| details | `card_url`, `tier`, `assay`, `product`, `chemistry_version`, `sample_id`, `study_id`, `canonical_key`, `development_stage`, `capture_area`, `default_table`, `n_transcripts`, `size_bytes`, `publication_year`, `license_url`, `collections` in `metadata` |
| entry | `datasetType: "spatialdata"`, `s3BaseUrl` = `zarr_url` rewritten from `s3://bucket/key` to the bucket's HTTPS URL (region from a bucket table; unknown bucket → the sync fails) |
| provenance | `sddbUid` (upsert key), `sourceRelease` (the tag), `sourceGeneratedAt` (`sddb_generated_at` from the file) |

Every store is also registered as an admin-owned `Dataset` row through the
same path as the viewer's "Add this dataset" banner
(`ensureDatasetForS3Url`: one read of the store's root metadata, the table's
size taken from it), named and described from the catalog, and the entry
carries its `datasetId`, so views are counted and the banner never shows.
A card opens the dataset's page (`/explore/[id]`, as every Explore card
does; `lib/explore/links.ts`), and the page opens the store in
`/spatialdata-viewer/from-s3?url=…` like any public store. Their bucket
allows CORS for `www.merfisheyes.com` only, so the stores open on
production and not on dev or localhost until the scverse team adds those
origins; the cards and pages work everywhere.

## The sync

`scripts/sync-sddb-catalog.ts`, run by hand (dev database first, then
production; no schedule):

1. Finds the newest release carrying `catalog.parquet` (or `--tag`), caches
   it under `.cache/sddb/<tag>/` (tags are immutable; `--file` reads a local
   copy).
2. Reads it with hyparquet (`lib/sddb/catalog-file.ts`): a missing required
   column or an all-null one fails the run; a missing known column warns.
3. Registers each store as an admin-owned `Dataset` (a few in flight at a
   time; a store whose root cannot be read is logged and its card synced
   without a dataset id).
4. Upserts every row by `sddbUid` and keeps its one entry in step.
5. Unpublishes rows whose uid left the release (never deletes: views and any
   later curation survive).
6. Gate: the published `sddb` uid set equals the release's uid set, else the
   run exits non-zero.

`--limit=N` syncs a sample spread over the technologies (a trial on the dev
database); a partial run unpublishes nothing and skips the gate.

A second run on the same release is a no-op. Idempotence and the gate are
what make the file, not our database, the source of truth.

## On Explore

- `source` is the column that tells the four kinds apart (`curated`, `bil`,
  `community`, `sddb`); the older flags stay.
- SpatialData-DB rows are mixed into "All" and have their own category chip.
  Filters (species, tissue, platform) include them.
- The licence is on every card that has one; rows from other sources have
  none yet.

## Deferred (recorded in the plan)

The static, content-hashed explore index behind a CDN, the gene-query parity
test against their `GeneIndex`, a `/spatialdata-viewer/sddb/<uid>` route, a
Lamin proxy backend, the nightly perf run and a cron: not needed to list
and open the datasets, and each is its own step.
