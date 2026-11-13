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
(2026-11-14): planned after the meeting with Tim Treis; no step started.**
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
  but not yet used to validate.
- **SDDB rows on Explore** (`docs/specs/sddb-catalog.md`): a subset of
  their columns mapped to ours, the rest kept in `metadata`.

## 3. Decisions

| # | Decision | Status |
|---|---|---|
| M1 | MERFISHEYES uses SpatialData-DB's catalog schema (version 3, pinned; column names, vocabularies, ontology sources and versions) for every dataset. Our additions are new columns in a `meyes_` namespace, never a redefinition of theirs | **taken** 2026-11-14 |
| M2 | Where metadata lives: **the store** (root attributes, SDDB's keys) — the copy that travels; **Lamin `eyes-catalog`** — the validated record, with verification; **the query index** — derived, what queries read; **Postgres** — the app's own database, which Lamin does not replace: accounts, ownership, projects, upload and processing state, views, private datasets, and the site's per-request reads (Lamin is a separate service the site would call on every page; it records datasets, not users) | **taken** 2026-11-14 |
| M3 | Terms are validated with bionty against the same ontology sources and versions SDDB pins, in the pipeline (Python) on the way in. A dataset published through our publishing system is public once bionty passes it; a person's review is expected but not required — bionty is the gate when nobody reviews | **taken** 2026-11-14 |
| M4 | One query index of every public dataset, in SDDB's catalog format plus our columns, published with each change; the site's filters and the MCP's queries start there | **taken** 2026-11-14 |
| M5 | A gene index for our datasets in SDDB's `gene_index` / `gene_search` format, so "datasets that measure / express X" is one lookup across both libraries | **taken** 2026-11-14 |
| M6 | A dataset shown in one view keeps all its sections in one store and is **one catalog row**, its sections listed (`meyes_sections`); after reprocessing, a dataset's cells and molecules are in that same store. Sections shown separately are separate stores and rows | **taken, tentatively** 2026-11-14 — the user is not certain; revisit when the reprocessed BIL stores and SDDB's per-sample columns are compared |
| M7 | Private datasets are never in the public index; an owner's queries read Postgres with the same columns | **taken** 2026-11-14 |

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

### Step 2 — Validation with bionty

`scripts/metadata/validate_metadata.py` (in `meyes-prepare`'s
environment): reads a dataset's metadata, validates each ontology field
with bionty against the pinned sources, standardises synonyms ("mouse"
→ *Mus musculus*, NCBITaxon:10090), fills the hierarchy columns from the
ontology (tissue system and organ, disease class and group — SDDB's
`ontology_rollup.json` logic), and writes a report. `prepare.py` calls it;
an invalid term is exit code **4** (invalid), with the closest terms in
the message.

**Gate:** the tiny fixtures pass; a wrong tissue fails with suggestions;
the result equals SDDB's for a handful of their own rows.

### Step 3 — Into the store

`prepare.py` and the BIL pipeline write the validated metadata into the
store's root attributes (SDDB's keys); `validate_store.py` checks it
against the JSON Schema; the viewer shows it (an info panel) for any
store, including one opened from disk.

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

### Step 5 — Backfill what exists

The ~87-dataset BIL fleet, community submissions and older public
datasets: an agent proposes the ontology terms, ages, genotypes and
conditions from titles, descriptions, BIL and GEO records, each with the
source text it came from; a person confirms or corrects (in Lamin or a
review sheet). Unconfirmed fields stay null, not guessed. Clears MCP.md
§10's E1–E14 where they are metadata.

**Gate:** every public row has organism, tissue and assay ids; the share
of rows with age and disease recorded; the review list empty or each
item explained.

### Step 6 — The query index

`scripts/metadata/build_index.py`: from Lamin (ours) and SDDB's latest
release (theirs) to `index/catalog.parquet` in SDDB's format plus the
`meyes_` columns, with `ontology_rollup.json` and `technology.json`
alongside, versioned like their releases (`release_info.json`), published
to the public bucket on each change. Postgres's catalog columns are
rebuilt from it (a migration adds the columns). The site's Explore and
`/api/agent/v1/catalog/search` query the index, not pages.

**Gate:** an Explore filter (organism + tissue with descendants + age
range + technology) answered from the index, timed against BUDGETS.md;
the numbers equal Postgres's; the agent's search returns the same rows.

### Step 7 — Genes across both libraries

Our pipeline's statistics give a gene index in SDDB's format for our
stores (`stats.json` → `gene_index` rows: total counts and fraction of
cells detected); `gene_search.parquet` for both libraries together.

**Gate:** "mouse brain datasets that express X above 5 % of cells"
answered across SDDB and MERFISHEYES from the index files alone.

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
- How often the index is rebuilt (on each change, or a nightly release
  like SDDB's).
- Which bionty source versions SDDB pins, so ours match term for term.
- Whether SDDB's `metadata_general` keys are meant to become a scverse
  convention — if so, follow it there.

## 6. Risks

- **Two schemas drifting**: pinning SDDB's version and checking our index
  with their validator (`sddb` client) on every build.
- **Backfill guesses**: proposals carry their source text and stay null
  until confirmed.
- **Index size**: SDDB's gene index is 1 GB; queries that need it run on
  the server, the browser reads only `catalog.parquet` and
  `gene_search.parquet`-sized files.
