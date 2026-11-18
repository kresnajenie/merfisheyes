# Plan: one metadata schema for every dataset — SpatialData-DB's

Give every dataset MERFISHEYES serves the metadata SpatialData-DB gives its
own: one versioned schema, ontology-backed terms validated with bionty,
written **into the store** so it travels with the data, recorded in Lamin
with a verification step, and published as **one query index** of every
public dataset that the site's filters and the agents start from.

SpatialData-DB (Tim Treis, Theis lab) is the gold standard here; what
MERFISHEYES adds is the viewer. So this plan copies their schema, their
ontology sources and their release artifacts as closely as it can, and
departs only where our data does not fit (multi-section stores,
molecules at billion scale) — each departure recorded in §3.

This file is the brief for every session working on it. **Status
(2026-11-18): steps 1–4 done; step 5 done for the 63 BIL stores in S3, the
rest waits on the reprocessing now running; steps 6 and 7 built and
published to `bil-psc-zarr/.index/` (each with a gate half still open);
steps 8–9 not started. `feature/metadata` is pushed, no PR yet; the worker
and cluster images are not rebuilt.** Spec: [`docs/specs/metadata.md`](../specs/metadata.md).
Keep it current: a step that changes a fact or a decision edits this
file in the same commit. Branch: `feature/metadata`, from `develop` (the
pipeline parts touch `scripts/` and `scripts/bil-scripts/`, so they
rebuild the worker image — CLAUDE.md §8). Related plans:
[MCP.md](MCP.md) (pillar B queries this index; the remote MCP serves it),
[LAMIN.md](LAMIN.md) (the `merfisheyes/eyes-catalog` instance),
[SDDB.md](SDDB.md) (reading their catalog), [BIL-CLUSTER.md](BIL-CLUSTER.md)
(the pipeline that writes our stores).

---

## 1. What the meeting settled (2026-11-14)

1. **Copy SpatialData-DB's metadata schema as far as it goes**, including
   the outside ontologies: tissue ids (UBERON), disease (MONDO), organism
   (NCBITaxon), assay (EFO), development stage (HsapDv / MmusDv). The site
   is restructured to filter on them.
2. **Validate with bionty**: incoming metadata is checked against the
   ontologies on the way in, and a verification step (in Lamin) confirms
   it before it is public.
3. **Metadata lives in the store** — the root `zarr.json` attributes, as
   SDDB's stores carry `metadata_general` — so it stays with the dataset
   wherever the store moves (S3, another bucket, a laptop: pillar D of
   MCP.md opens it from disk with its metadata). The pipeline writes it.
4. **A query index file**: one file holding the quick facts of every
   public dataset; a filter on the site or an agent's query starts there,
   not at pages or the database. Public datasets are a facet of the site.
5. The local use case (MCP.md pillar D) is unchanged, except that a
   local store now carries its own metadata, which the viewer can show.

## 2. Facts

### SpatialData-DB, as of release `catalog-2026-11-12`

- **Three layers.** (a) Inside each store: root `zarr.json` attributes
  `metadata_general` (Assay, Organism, Tissue, Disease, Dataset Url,
  Publication Date, Default Table, Preservation Method, H&E / IF
  available, …), `sddb_provenance` (`sddb_id`, `canonical_key`,
  `converted_at`), `spatialdata_io_reader`, and a per-technology block
  (`metadata_xenium`, `metadata_visium`, …). (b) The Lamin instance
  `scverse/spatialdata-db`: every store an artifact annotated with
  ontology-backed features. (c) Release files on
  `theislab/spatialdata-db-client` (near-nightly `catalog-<date>` tags).
- **Release files** (`catalog-2026-11-12`): `catalog.parquet` (11 MB, the
  query index), `gene_index.parquet` (1 GB: dataset × gene,
  `total_counts`, `fraction_obs_detected`), `gene_search.parquet` (16 MB:
  one row per gene with a packed bitmap of the datasets that have it),
  `feature_index` / `feature_search` (the same for proteins and other
  modalities), `marker_index.parquet`, `ontology_rollup.json` (each
  ontology term with its ancestors present and name — how "brain" finds
  hippocampus), `technology.json` (vocabulary: readouts, vendors,
  technologies, variants; schema 3.0.0), `collections.parquet`,
  `citations.bib`, `release_info.json` (`catalog_schema_version: 3`,
  curation and hub commits, lamindb version, share of "other"
  tissue/disease), `release_stats.json`.
- **The schema** (`src/sddb/_catalog_schema.py`, schema 3): identity
  (`uid`, `canonical_key`, `study_id`, `collection_name`); biology
  (`organism`/`organism_id`, `tissue`/`tissue_id` with
  `tissue_system`/`tissue_organ`, `disease`/`disease_id` with
  `disease_class`/`disease_group`, `development_stage`); technology
  (`technology`, `vendor`, `variant` from `technology.json`, `assay`,
  `panel_size`); per sample (`sample_id`, `donor_id`, `timepoint`,
  `timepoint_kind`, `timepoint_days`, `treatment`, `genotype`,
  `condition`); curation (`validation_status`, `tier`, `review_status`,
  `review_reason`); licence; headline stats (`n_obs`, `n_features`,
  `total_counts`, `n_transcripts`); access (`zarr_url`, `card_url`,
  `vitessce_url`); observation unit and segmentation; publication (`doi`,
  `pmid`, verification status). Required: `uid`, `technology`, `assay`,
  `organism`, `validation_status`, `zarr_url`.
- **bionty** (Lamin): registries backed by public ontologies — Organism
  (NCBI Taxonomy, Ensembl), Gene (Ensembl, NCBI Gene), Protein (UniProt),
  CellType (Cell Ontology), CellLine, Tissue (Uberon), Disease (Mondo,
  Human Disease Ontology), Pathway (GO), Phenotype (HPO, Mammalian
  Phenotype, …), ExperimentalFactor (EFO), DevelopmentalStage (human and
  mouse), Ethnicity (HANCESTRO); `validate`, synonyms and abbreviations,
  versioned `Source`s.
- **What they pin** (read from `scverse/spatialdata-db` 2026-11-13; in
  `conventions.json` `metadata.ontologies`): NCBITaxon 2023-06-20 and
  Ensembl release-112 (organism names), UBERON 2024-08-07, MONDO
  2024-06-04, HsapDv / MmusDv 2024-05-28, EFO 3.70.0, PATO 2024-03-28, CL
  2024-05-15. bionty's defaults are newer (Ensembl release-116), so every
  lookup names its version.
- **Their stores carry labels, not ids.** `metadata_general` (General
  Metadata Version 0.3.0–0.3.2, 25 keys over 83 sampled stores) says
  `"Organism": "mouse"`, `"Tissue": "brain"`; no ontology id is in any
  store. The ids are in their Lamin features (named as the
  `metadata_general` keys, typed to bionty: `Organism`, `Tissue`,
  `Disease`, `Development Stage`, `Assay` → EFO) and in `catalog.parquet`.
  No sex column or feature exists.
- **They list multi-section data per section**: the Zhuang whole-brain
  MERFISH is one row per section (`canonical_key`
  `abc:Zhuang-ABCA-1.010`), grouped by `study_id` `abc:Zhuang-ABCA-1`;
  they tag it `merfish` (academic) and Yao et al. 2023's `vizgen_merscope`.
- **Their agent surface**: the Python client's CLI and stdio MCP server
  (`sddb-mcp`: `query`, `describe`, `genes`, `features`, `facets`,
  `cite`, `download`). Tim described a hosted, connector-style MCP at the
  site — not in their published docs yet (MCP.md §0).

### MERFISHEYES today

- **Postgres** (`Dataset`, `CatalogDataset`): free-text `species`,
  `tissue`, `platform`, `disease`, `tags`, and a JSON `metadata` blob. No
  ontology ids, ages, sex or genotype. Explore and the agent API search
  these; spellings are kept equal by `docs/specs/conventions.md`.
- **Sidecars next to a store** (BIL pipeline): `catalog.json` (BIL code,
  title, description, species, tissue, platform, links, free metadata),
  `stats.json` (counts, genes, molecules, DE columns, pipeline),
  `qc.json` — on reprocessed BIL stores only (MCP.md §10, E13). Nothing
  is written into the store's own attributes.
- **Legacy folders** carry `manifest.json`: technical, no biology.
- **Lamin** `merfisheyes/eyes-catalog` (LAMIN.md): BIL datasets
  registered by `scripts/bil-scripts/register_lamin.py`, bionty available
  but not yet used to validate. What is registered there now is not
  migrated: the reprocessed fleet is registered afresh (step 5).
- **SDDB rows on Explore** (`docs/specs/sddb-catalog.md`): a subset of
  their columns mapped to ours, the rest kept in `metadata`.
- **Production's own rows (2026-11-16)**: 92 BIL, 56 curated (47 of them
  BIL Xenium, `ace-owl-*`), 1 community. No BIL row opens a store yet (all
  legacy chunked folders); 2 curated rows do. All 92 BIL rows say
  platform "MERFISH", which does not name the instrument. 62 reprocessed
  BIL stores with sidecars exist on the workstation. Per-row counts:
  `docs/specs/metadata-inventory.md`.

## 3. Decisions

| # | Decision | Status |
|---|---|---|
| M1 | MERFISHEYES uses SpatialData-DB's catalog schema (version 3, pinned; column names, vocabularies, ontology sources and versions) for every dataset. Our additions are new columns in a `meyes_` namespace, never a redefinition of theirs | **taken** 2026-11-14 |
| M2 | Where metadata lives: **the store** (root attributes, SDDB's keys) — the copy that travels; **Lamin `eyes-catalog`** — the validated record, with verification; **the query index** — derived, what queries read; **Postgres** — the app's own database, which Lamin does not replace: accounts, ownership, projects, upload and processing state, views, private datasets, and the site's per-request reads (Lamin is a separate service the site would call on every page; it records datasets, not users) | **taken** 2026-11-14 |
| M3 | Terms are validated with bionty against the same ontology sources and versions SDDB pins. A dataset published through our publishing system is public once bionty passes it; a person's review is expected but not required — bionty is the gate when nobody reviews | **taken** 2026-11-14. **Amended 2026-11-17**: bionty lives on the Lamin side, not in the pipeline. The pipeline writes the metadata unvalidated (`not_validated`); `validate_metadata.py`, run where lamindb is, checks it before the store is registered and writes the verdict into the store and Lamin |
| M4 | One query index of every public dataset, in SDDB's catalog format plus our columns, published with each change; the site's filters and the MCP's queries start there | **taken** 2026-11-14 |
| M5 | A gene index for our datasets in SDDB's `gene_index` / `gene_search` format, so "datasets that measure / express X" is one lookup across both libraries | **taken** 2026-11-14 |
| M6 | A dataset shown in one view keeps all its sections in one store and is **one catalog row**, its sections listed (`meyes_sections`); after reprocessing, a dataset's cells and molecules are in that same store. Sections shown separately are separate stores and rows | **taken, tentatively** 2026-11-14 — the user is not certain; revisit when the reprocessed BIL stores and SDDB's per-sample columns are compared. 2026-11-16: kept for step 1 with per-section facts nested in `meyes_sections`, knowing SDDB lists the Zhuang ABCA data per section; still revisited before step 6 |
| M7 | Private datasets are never in the public index; an owner's queries read Postgres with the same columns | **taken** 2026-11-14 |
| M8 | Our stores carry the ontology ids: `metadata_general` and the per-technology block exactly as SDDB writes them (labels), and a sibling root attribute `meyes_metadata` with every id and our columns under the index's column names, plus `meyes_provenance`. SDDB's blocks are never redefined | **taken** 2026-11-16 |
| M9 | A row enters the public index only when its dataset is a store (`zarr_url` is required). Legacy chunked datasets keep the same columns in Postgres and stay on Explore until reprocessed | **taken** 2026-11-16. **Amended 2026-11-18**: every dataset live on the site is in the index; one without a checked store carries what Explore knows (`meyes_metadata_from` `explore`, `meyes_format` `merfisheyes-chunked` when it has no store, `validation_status` `not_validated`), mapped by the pipeline's own rules |
| M10 | All our data is reprocessed for this schema (step 5). Until that is done the required set is relaxed to what search and the MCP need — organism and tissue with ids, counts, title, a store and our links (16 columns); `uid`, `technology`, `vendor`, `assay` are recommended and become required at step 5's gate | **taken** 2026-11-16 |
| M11 | The tissue system / organ and disease class / group facets use our own anchor lists (`conventions.json` `metadata.groupings`), built from the terms SDDB's release uses; first anchor in list order among the term and its ancestors wins | **taken** 2026-11-16 |
| M12 | Every store gets a count-matrix hash, `meyes_matrix_hash` (`mh1:`, `metadata.md` §7.7; reference `scripts/metadata/matrix_hash.py`), whole and per section, so re-releases are found: equal hashes group rows as SpatialData-DB's `same_section_group` / `rerelease`; near matches (the hash's cell sketch, per-gene totals from the gene index) are proposals a person confirms | **taken** 2026-11-16 |
| M13 | The public index is an account of what is **live on MERFISHEYES**: a dataset is in it only when a public Explore row lists its store, so it can be opened and analysed on the site. SpatialData-DB's release is a reference: their rows enter only when Explore lists the same store. Published in the stores' own folder, `s3://merfisheyes-bil/bil-psc-zarr/.index/<date>/` and `.index/latest/` (a dot folder, as `.lamindb/`, so listings of the stores skip it) | **taken** 2026-11-18 |
| M14 | Personal indexes for private datasets are **deferred** until our data is on Lamin's R2 storage, whose infrastructure is expected to give each person their own database; the index goes there (outline in step 6, "Later: personal indexes"). Until then M7 stands: private datasets are in no index file | **taken** 2026-11-18 |
| M15 | The index is rebuilt and republished after each batch of datasets is published, and nightly to pick up SpatialData-DB's latest release (bumping the pin after re-reading their schema). Each build is kept as `.index/<date>/`; `.index/latest/` is the current one. The nightly job is not set up yet | **taken** 2026-11-18 |

## 4. Steps

Each step ends at a gate; the next does not start until it passes.

### Step 1 — The schema on paper

`docs/specs/metadata.md`: our profile of SDDB schema 3 — which of their
columns a MERFISHEYES row must fill (required, recommended, optional),
the `meyes_` additions (BIL code, MERFISHEYES ids and viewer URLs,
sections, molecule counts, DE columns, the pipeline that wrote it), the
store attribute layout (`metadata_general` and the per-technology block
as SDDB writes them, plus `meyes_provenance`), and the map from today's
fields (`species` → `organism` + `organism_id`, …). A JSON Schema beside
it; the pinned schema version and ontology source versions in
`conventions.json`. Replaces the spelling rules in `conventions.md` with
ontology ids.

**Gate:** every current BIL and community row mapped on paper, with the
count of fields each can fill today; reviewed against SDDB's
`_catalog_schema.py` column by column.

*Drafted 2026-11-16, awaiting review.* The profile covers all 93 of their
columns (each with a level; `test_metadata_schema_profiles_sddb_schema_3`
fails if one is missing or its dtype differs) plus 26 `meyes_*` columns;
16 are required while the fleet is reprocessed (M10), 20 after. All 149
production rows are counted in `metadata-inventory.md`: BIL rows fill a
median of 13 of 16 required columns today and all 16 once they open their
reprocessed store; what they lack today is the store itself.
`technology`/`vendor`/`assay` need the instrument, which "MERFISH" does
not name — recommended until reprocessing. The spelling maps in `conventions.json`
(`catalog.species`, `catalog.platforms`) stay until step 6, when Explore
reads the index; `metadata.md` §7.1 and §7.3 are what replaces them.

### Step 2 — Validation with bionty

`scripts/metadata/validate_metadata.py`, **on the Lamin side** (the
`meyes-lamin` environment; M3 as amended 2026-11-17 — bionty is not in the
pipeline): reads a store's `meyes_metadata`, checks organism, tissue,
disease, sex and development stage against the pinned sources by name,
synonym or id, standardises them ("Macaque" → `macaque`, NCBITaxon:9544),
fills `development_stage` from `timepoint_days`, and with `--write` puts the
standardised values and the verdict (`meyes_provenance.validation`) back
into the store. Exit 0 pass, 3 incomplete (organism or tissue unknown), 4
invalid, with the closest terms in the message. The hierarchy columns
(M11) are the index build's (step 6), where the full ontology graph is.

**Gate:** a wrong tissue fails with suggestions; the BIL fleet passes.

*Done 2026-11-17.* `tests/python/test_metadata_bionty.py` (5 tests, run in
`meyes-lamin`; skipped where bionty is absent). All 62 reprocessed BIL
stores' metadata (built, not written) pass: 54 mouse stores get a
development stage, macaque and marmoset have no stage ontology, and two
human ages (`234.0 days`, `NA Gestational weeks`) match no term — backfill
items. bionty's search is exact, so suggestions come from a fuzzy match over
names and synonyms. Several Ensembl entries share one NCBITaxon id (mouse
strains): the name that agrees with the id is kept. Not compared yet with
SDDB's values for their own rows.

### Step 3 — Into the store

`prepare.py` and the BIL pipeline write the metadata into the
store's root attributes (SDDB's keys), with the count-matrix hash (M12,
`scripts/metadata/matrix_hash.py`), whole and per section; `validate_store.py` checks it
against the JSON Schema; the viewer shows it (an info panel) for any
store, including one opened from disk.

The pipeline side (2026-11-17): `scripts/metadata/build_metadata.py OUT`
writes the four blocks from `catalog.json`, `qc.json`, `stats.json`, the
store and an optional `--metadata FILE` (what a person supplies, as index
columns), mapping spellings to ids only through `conventions.json`
`metadata.terms` / `technologies` / `technologyFrom`, and marks it
`not_validated`. `prepare.py` calls it after `stats.json` (new
`--metadata`; a failure is logged, the exit code is unchanged) — the worker
image now copies `scripts/metadata/` (rebuild it, CLAUDE.md §8). The BIL
pipeline has a new **metadata** stage (`metadata_store.sbatch`) after
points, qc and catalog and before validate, since the mask, MapMyCells and
points stages change the store after `prepare.py`; the cluster's
apptainer image needs `scripts/metadata/` too. `validate_store.py` checks
the blocks (WARN when missing or not validated; FAIL with
`--require-metadata`, and always when `n_obs` × `n_features` disagree with
the table). Dry run over the 62 local stores: all build, all pass the
schema; 6 get `vizgen_merscope` from `input: merscope`, 56 have no
technology yet; 42 s (243 k cells) to 235 s each, the matrix hash being
most of it. Test: `test_build_metadata_writes_the_four_blocks`.

The viewer side (2026-11-17): `readStoreMetadata` (`lib/spatialdata/metadata.ts`,
`tests/unit/store-metadata.test.ts`) reads the four blocks from the root
`zarr.json` the store is opened with, ours over SpatialData-DB's labels, and
the info button ("About this dataset") in each viewer's top-right row shows
them: by id, by URL, and from disk (`/viewer|sm-viewer|image-viewer/from-local`).
A SpatialData-DB store shows its `metadata_general` labels.

**Gate:** a store copied to another bucket, and to a laptop, opens with
its metadata shown; `validate_store.py --mapping` passes on the BIL
fleet's reprocessed stores.

### Step 4 — Lamin: the record, and verification

`register_lamin.py` annotates each artifact with the validated features
(bionty registries), sets `validation_status`, and leaves a review state
a person may clear (who reviews is decided in this step). A dataset
published through our publishing system is in the index once bionty
passes it (M3); an open review marks it, it does not hide it.

**Gate:** in `eyes-catalog`, a query by organism and tissue id returns
the BIL datasets; a dataset bionty rejects is held back.

*Built 2026-11-17, not yet run against eyes-catalog.* `register_lamin.py`
reads the store's root `zarr.json` from S3 and always records
`validation_status` (`pass` | `fail` | `not_validated` | `absent`),
`canonical_key`, `meyes_matrix_hash`, `technology` and the plain columns
(`Genotype`, `Timepoint`, `Timepoint Days`, `Study ID`, …); the typed
features (`Organism`, `Tissue`, `Disease`, `Development Stage`, `Assay`,
`meyes_sex`) only when the store's validation passed. eyes-catalog's current
bionty sources are newer than SDDB's pins (UBERON 2026-06-23, MONDO
2026-09-01, Ensembl release-116, …) and the pinned versions are not in it:
the script adds them as non-current `bt.Source`s and makes the records from
them, leaving the instance's defaults alone. `--dry-run` shows each store's
validation status and ids. Who reviews, and the review state, are still
open.

*First store, 2026-11-17:* ace-dip-van's metadata built, passed bionty,
its root `zarr.json` copied to S3 (single key, byte-identical to the local
store), registered: artifact `N8eVbF41CKh2RHE30001` (a new version — the
root file changed), every check passing. `ln.Artifact.filter(Organism=mouse,
Tissue=brain)` returns it; its features: mouse, brain, 13-week-old stage,
male, validation `pass`, `bil:ace-dip-van`, `vizgen_merscope`, the matrix
hash. The five pinned Sources were added non-current; the instance's
defaults are unchanged. **SDDB's pins disagree with each other**: 12 of the
20 EFO ids in their `technology.json` (MERSCOPE EFO:0022956, Visium HD,
CosMx, Stereo-seq, …) are not in the EFO version they pin (3.70.0), so the
typed `Assay` is left unset for those (the script says so and goes on;
`technology` holds the id) — a question for Tim. The other 61 stores wait.

### Handoff to the BIL pipeline sessions (2026-11-17)

For the sessions that reprocess the fleet (`bil/pipeline-dry-run`; the
"spatialdata pipeline" and "Restore BIL processed data from s3" chats):

- **What is in S3 now** (`s3://merfisheyes-bil/bil-psc-zarr/`): 63 stores —
  `ace-dip-vat` from raw on the cluster, 62 rebuilt from legacy output
  (`reconstruct_inputs.py` → `~/data/bil-data/rebuilt/`, `prepare.py` →
  `~/data/bil-data/processed/<code>/store/`). The local copies equal S3 file
  for file (root `zarr.json` and the four sidecars, checked 2026-11-17).
  None carries metadata yet.
- **What this branch adds to `scripts/bil-scripts/`** (identical files on
  both branches until now): a `metadata` stage (`metadata_store.sbatch`,
  after points, qc and catalog, before validate) running
  `scripts/metadata/build_metadata.py`; `register_lamin.py` reads the
  store's metadata and sets typed features. Between sync and register,
  run `scripts/metadata/validate_metadata.py --write` (needs lamindb, i.e.
  the cluster's Lamin env) on the store before it is synced — or on the
  local copy, then sync the one changed `zarr.json`.
- **The instrument is lost in the rebuild**: 56 of the 62 went through
  `prepare.py` as `input: folder`, so `technology` is unknown for them.
  If `reconstruct_inputs.py` knows the source was MERSCOPE output, record
  it (an input `metadata.json` with `"technology": "vizgen_merscope"`, passed
  as `prepare.py --metadata`), else it is a backfill item.
- **Images**: the worker image and the cluster's apptainer image must
  include `scripts/metadata/`.
- **Updating metadata on a store already in S3 needs a single-key copy, not
  a sync** (decided 2026-11-17: `s3_sync_store.sbatch` keeps
  `aws s3 sync --size-only`). The sync compares sizes only, so an in-place
  rewrite of a store's root `zarr.json` that keeps its length — a verdict
  flipped by `validate_metadata.py --write`, a corrected id, a re-release —
  is skipped with no error, and S3 keeps the old bytes. After any such edit:
  `aws s3 cp <store>/zarr.json s3://merfisheyes-bil/bil-psc-zarr/<code>/<name>.zarr/zarr.json --content-type application/json`,
  then re-register.

### Step 5 — Reprocess everything, and backfill on the way

*Started 2026-11-17 for the 63 stores already in S3*, which need metadata, not
reprocessing: `scripts/bil-scripts/metadata_backfill.sh prepare|publish
<folder>…` builds, checks with bionty and validates each store locally, then
(`publish`) uploads each changed root `zarr.json` as a single key and
registers the passing stores. The technology comes from the title for the 53
"MERSCOPE v1" stores (`metadata.md` §7.3) and from
`scripts/bil-scripts/metadata-overrides.json` for the three Xenium marmoset
stores.

*Done for the 63, 2026-11-17* (run logs in `~/data/bil-data/metadata-runs/`):
all 63 build, pass bionty and `validate_store --require-metadata`; every root
`zarr.json` in S3 equals its local store; all 63 registered in eyes-catalog
(new artifact versions; 126/126 artifacts and sheet rows pass every check),
validation `pass`. A query by organism and tissue id returns the 55 mouse
brain stores; human 3, macaque 2, marmoset 3. Technology: 60
`vizgen_merscope`, 3 `10x_xenium`. The typed `Assay` is set only on the 3
Xenium stores (SDDB's MERSCOPE EFO id is not in their pinned EFO). Left
for the backfill: disease (all `unknown`), development stage for macaque,
marmoset and two human stores, donor, DOI.

**Every dataset we serve is reprocessed** into a store carrying this
metadata (M10): the BIL fleet (92 rows, plus the 47 `ace-owl-*` Xenium
rows listed as curated), the other curated rows and community
submissions. The pipeline on `bil/pipeline-dry-run` is updated to write
the metadata (step 3) and to sync and register each store in Lamin (step
4); what eyes-catalog holds today is replaced, not migrated.

Before a dataset's run, its metadata is backfilled into its input
(`catalog.json`): an agent proposes the ontology terms, the instrument
(`merfish` or `vizgen_merscope`), ages, genotypes and conditions from
titles, descriptions, BIL and GEO records, each with the source text it
came from; a person confirms or corrects (a review sheet). Unconfirmed
fields stay null, not guessed. Clears MCP.md §10's E1–E14 where they are
metadata. A legacy chunked dataset enters the index (M9) once its store is
published.

**Gate:** every public dataset is a store with the metadata, registered;
every public row has organism and tissue ids and a technology; the share
of rows with age and disease recorded; the review list empty or each item
explained. Then the interim requirement ends: `uid`, `technology`,
`vendor`, `assay` become required (`metadata.md` §3).

### Step 6 — The query index

`scripts/metadata/build_index.py`: from Lamin (ours) and SDDB's latest
release (theirs) to `index/catalog.parquet` in SDDB's format plus the
`meyes_` columns, with `ontology_rollup.json` and `technology.json`
alongside, versioned like their releases (`release_info.json`), published
to the public bucket on each change. Rows with equal matrix hashes are
grouped as re-releases (`same_section_group`, M12). Postgres's catalog columns are
rebuilt from it (a migration adds the columns). The site's Explore and
`/api/agent/v1/catalog/search` query the index, not pages.

*Built and published 2026-11-18* (M13), at `.index/2026-11-18/` and `.index/latest/`. `scripts/metadata/build_index.py OUT`
(meyes-lamin) writes `catalog.parquet` (Snappy, as SDDB's), `ontology_rollup.json`,
`technology.json` and `release_info.json`. Rows: the SDDB rows of `catalog-2026-11-17`
whose store Explore lists (12,575, unchanged, `meyes_source` `sddb`; their other
966 are not on the site and are left out) plus ours — every eyes-catalog store whose
metadata passed bionty and that a public Explore row lists (63: 60 `bil`, 3
`curated`). 102 public Explore rows have no store yet (legacy, M9) and are
in it too since M9 was amended: what Explore knows (from each row's detail, through `build_metadata.from_catalog`), `meyes_metadata_from` `explore`. Checked with SDDB's own validator (their `_catalog_schema.py` at
the tag) and our rows against `$defs/indexRow`. The pin moved to
`catalog-2026-11-17`, whose schema added 10 columns (still version 3); our
profile now lists 103 SDDB columns. **12,740 rows — every live Explore row**
(`every_live_row_indexed` in `release_info.json`); rebuilt per M15.

**Later: personal indexes (M14 — not built; build once the data is on Lamin's R2).**
The same `catalog.parquet` format, one per person: the datasets they own or that
are shared with them through projects, from Postgres `Dataset` rows joined with
each store's root `zarr.json`. Rebuilt on every change to that person's datasets
(upload finished, metadata edited, deleted, shared). Stored in their own space —
on Lamin's R2 per-person database when that exists (to be confirmed with Lamin),
never in the public `.index/`. Read only through a signed-in route (e.g.
`GET /api/agent/v1/me/index`: the session or the MCP's OAuth token, then a
short-lived link), so the MCP runs one query over the public index and the
person's own. Until then, private search reads Postgres (M7).

**Gate:** an Explore filter (organism + tissue with descendants + age
range + technology) answered from the index, timed against BUDGETS.md;
the numbers equal Postgres's; the agent's search returns the same rows.

### Step 7 — Genes across both libraries

Our pipeline's statistics give a gene index in SDDB's format for our
stores (`stats.json` → `gene_index` rows: total counts and fraction of
cells detected); `gene_search.parquet` for both libraries together. The
per-gene totals also flag likely re-releases across both catalogs
(near-equal vectors; a person confirms, M12).

**Gate:** "mouse brain datasets that express X above 5 % of cells"
answered across SDDB and MERFISHEYES from the index files alone.

*Built and published 2026-11-18* (`scripts/metadata/build_gene_index.py`;
`metadata.md` §7.8), beside the index in `.index/<date>/` and `latest/`:
`gene_index.parquet` 63,710 rows over 151 of our 165 datasets (29,500 with
counts, from the 63 stores; 14 Explore-only datasets have no gene list),
`gene_search.parquet` 764,753 entities over all 12,740 rows (4,522 of our
genes joined SpatialData-DB's entities, 205 new), `gene_search_meta.json`.
Checks: their validator passes `gene_index`; their bitsets survive the
re-ordering exactly (cow GAD1 → the same 4 datasets); "mouse brain stores
expressing Rab3b in > 5 % of cells" → 55. The gate's cross-library half
needs their `gene_index.parquet` (1 GB) read by row group, not done.

### Step 8 — A site built to filter

Explore gets the facets the index carries: organism, tissue as a tree
(system → organ → term), disease class, technology and vendor, age range
in days, sex, genotype, condition, genes measured or expressed, tier.
Upload and edit forms choose terms from the ontologies (searched from a
pinned term list), and the worker validates them with bionty before a
dataset is public.

**Gate:** a person finds "mouse hippocampus, 1–12 months, measures Gad1"
in under the BUDGETS.md interaction budget; an upload with a free-text
tissue is turned into a term or held for review.

### Step 9 — Give it back

Our index in SDDB's format is something SpatialData-DB could ingest:
offer it, and a `merfisheyes_url` column for their catalog, when Tim
agrees.

## 5. Open questions

- M6 is taken tentatively: check it against the reprocessed BIL stores
  and SDDB's per-sample columns (`sample_id`, `timepoint`) before step 6.
- Who reviews metadata, and where (Lamin's UI, a sheet, an admin page) — review is expected, not required (M3).
- ~~How often the index is rebuilt~~ — M15 (the nightly job is still to set up).
- ~~Where the index is published~~ — `bil-psc-zarr/.index/` (M13).
- ~~Which bionty source versions SDDB pins~~ — answered 2026-11-16 (§2;
  `conventions.json` `metadata.ontologies`).
- Whether SDDB's `metadata_general` keys are meant to become a scverse
  convention — if so, follow it there.
- Nice to have from Tim, not blocking: their anchor tables (to compare
  with ours, M11) and the `cmh1:` count-matrix hash algorithm. Their
  `tier` (bronze/silver) grades their own conversions; we leave it null.
- The instrument behind each BIL "MERFISH" dataset (`merfish` or
  `vizgen_merscope`): part of step 5's backfill.

## 6. Risks

- **Two schemas drifting**: pinning SDDB's version and checking our index
  with their validator (`sddb` client) on every build.
- **Backfill guesses**: proposals carry their source text and stay null
  until confirmed.
- **Index size**: SDDB's gene index is 1 GB; queries that need it run on
  the server, the browser reads only `catalog.parquet` and
  `gene_search.parquet`-sized files.
