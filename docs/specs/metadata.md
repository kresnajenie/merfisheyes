# Dataset metadata: our profile of SpatialData-DB catalog schema 3

Status: **draft, 2026-11-16** — step 1 of [`docs/plans/METADATA.md`](../plans/METADATA.md).
Measured against SpatialData-DB release `catalog-2026-11-12` (12,575 rows,
`catalog_schema_version: 3`, `technology.json` 3.0.0), the root `zarr.json` of
83 of their stores (4 per technology), their Lamin instance
`scverse/spatialdata-db`, and every public curated, BIL and community row of
MERFISHEYES (149; [`metadata-inventory.md`](metadata-inventory.md)).

Machine-readable copies: [`schemas/metadata.schema.json`](schemas/metadata.schema.json)
(the store attributes and one index row, every column with its level) and the
`metadata` key of [`conventions.json`](conventions.json) (what is pinned).
`tests/python/test_specs.py::test_metadata_schema_profiles_sddb_schema_3`
checks the two against each other and against SpatialData-DB's required set.

SpatialData-DB is the reference. A column of theirs means what it means there;
ours are named `meyes_*` (M1). Every place we depart is listed in §9.

---

## 1. What is pinned

| What | Pinned to | Where |
|---|---|---|
| Catalog schema | version **3** (`theislab/spatialdata-db-client` `src/sddb/_catalog_schema.py`, identical at tag `catalog-2026-11-12` and `main` on 2026-11-13) | `conventions.json` `metadata.sddb` |
| Technology vocabulary | `technology.json` **3.0.0** of that release | same |
| In-store block | `metadata_general` **General Metadata Version 0.3.2** (their stores carry 0.3.0–0.3.2) | same |
| Per-technology block | `metadata_merscope` 0.1.0 (MERFISH and MERSCOPE), `metadata_xenium` 0.2.0 | same, `technologyBlocks` |
| Ontologies | the bionty `Source`s SpatialData-DB marks `currently_used` (read 2026-11-13) | `conventions.json` `metadata.ontologies` |
| Licences | their five ids we may meet plus two of ours (§9 X4) | `conventions.json` `metadata.licenses` |
| Facet groupings | our own anchor lists for tissue system / organ and disease class / group, built from the terms their release uses (§7.6) | `conventions.json` `metadata.groupings` |

The ontology sources, from `scverse/spatialdata-db`:

| Columns | bionty registry | Source | Version |
|---|---|---|---|
| `organism_id` | Organism | NCBITaxon | 2023-06-20 |
| `organism` (the name) | Organism | Ensembl (vertebrates) | release-112 |
| `tissue`, `tissue_id` | Tissue | UBERON | 2024-08-07 |
| `disease`, `disease_id` | Disease | MONDO | 2024-06-04 |
| `development_stage` | DevelopmentalStage | HsapDv (human), MmusDv (mouse) | 2024-05-28 |
| `assay` (EFO id via `technology.json`) | ExperimentalFactor | EFO | 3.70.0 |
| `meyes_sex_id`, `normal` | Phenotype | PATO | 2024-03-28 |
| cell types (later) | CellType | CL | 2024-05-15 |
| genes (step 7) | Gene | Ensembl | release-112 |

bionty's own defaults are newer (the installed one reads Ensembl release-116),
so every lookup names its source and version explicitly. Organism names come
from Ensembl's list, lower-cased, as SpatialData-DB writes them: Macaca
mulatta is `macaque` there and `rhesus macaque` in NCBITaxon; Callithrix
jacchus is `white-tufted-ear marmoset` in both.

## 2. Where each piece lives (M2)

| Layer | Holds | Written by | Read by |
|---|---|---|---|
| **The store**, root `zarr.json` attributes | `metadata_general` and `metadata_<technology>` as SpatialData-DB writes them (labels); `meyes_metadata` (the ids and our columns); `meyes_provenance` (§4) | the pipeline (`prepare.py`, the BIL pipeline), step 3 | the viewer's info panel (any store, from S3 or disk), `validate_store.py` |
| **Lamin** `merfisheyes/eyes-catalog` | one artifact per store, annotated with ontology-typed features under SpatialData-DB's feature names (§5); `validation_status`; the review state | the pipeline's sync, once it writes this metadata (steps 4–5); what is registered today is replaced then | the index build; people reviewing |
| **The query index** `catalog.parquet` | one row per public store: §3's columns, SpatialData-DB's first then `meyes_*` | `build_index.py`, step 6 | Explore's filters, `/api/agent/v1/catalog/search`, the MCP |
| **Postgres** | the same columns on `CatalogDataset` (public) and `Dataset` (private, M7), plus everything only the app needs: owner, projects, publishing and upload state, views | rebuilt from the index (public), the owner's edits (private) | the site's per-request reads |

Columns that are facts about the data travel with the store (`meyes_metadata`).
Columns about where and when it is published — `uid`, `zarr_url`, `card_url`,
`size_bytes`, `published_at`, `validation_status`, `meyes_id`, `meyes_source`,
`meyes_page_url`, `meyes_viewer_url` — are in the index and Lamin only: a copy
of the store on a laptop is the same data, not the same publication.

**A row needs a store** (decided 2026-11-16). `zarr_url` is one of
SpatialData-DB's required columns; a legacy chunked dataset has no store, so it
enters the index when it is reprocessed. Until then Explore shows it from
Postgres with whatever columns it has. **Every dataset we serve will be
reprocessed** by the pipeline that writes this metadata (METADATA.md step 5):
the stores, their metadata and their Lamin records all come from that run.

## 3. The index row

Levels, for a MERFISHEYES row (the schema's `x-meyes-level`):

- **required** — the row is not published without it. **Interim set (M10),
  16 columns**, while the fleet is reprocessed: `organism`, `organism_id`,
  `tissue`, `tissue_id`, `validation_status`, `zarr_url`, `canonical_key`,
  `n_obs`, `n_features`, `published_at`, `title`, and five `meyes_*`
  (`meyes_id`, `meyes_source`, `meyes_kind`, `meyes_page_url`,
  `meyes_viewer_url`) — what a search by organism and tissue and a link to
  the viewer need. **After reprocessing** four more join it (the schema's
  `x-meyes-required-after-reprocessing`): `uid`, `technology`, `vendor`,
  `assay`, which makes SpatialData-DB's six required columns ours again.
  Until then those four are recommended.
- **recommended** — filled whenever the source says it; a null is a backfill item.
- **optional** — filled when known; no one chases it.
- **derived** — computed by the index build from other columns; never stored by hand.
- **n/a** — SpatialData-DB's own bookkeeping, or a technology we do not serve; null.

`unknown` is a value, as in SpatialData-DB: `tissue` and `disease` are
`unknown` with a null id when the source does not say. It does not count as
filled in the inventory.

### 3.1 SpatialData-DB's columns (93, in their order)

**Identity**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `uid` | recommended → required after reprocessing | the store's artifact uid in eyes-catalog |
| `canonical_key` | required | `bil:<code>` for a BIL dataset (also the curated `ace-owl-*` rows); else `meyes:<meyes_id>` |
| `study_id` | recommended | `doi:<doi>` when the dataset is one of a paper's (SpatialData-DB's `doi:` prefix); else as `canonical_key`. Their client derives `source` from the prefix: ours read `bil` or `meyes` |
| `collection_name`, `collections` | optional | slugs grouping datasets of one atlas or paper (the Allen whole-brain BIL codes) — backfill |

**Biology**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `organism`, `organism_id` | required | §7.1 |
| `tissue`, `tissue_id` | required | §7.2 |
| `disease`, `disease_id` | recommended | MONDO; `normal` = PATO:0000461; `unknown` + null. No row says today (one curated row says "Breast cancer") |
| `development_stage` | recommended | §7.4 |
| `tissue_system(_id)`, `tissue_organ(_id)`, `disease_class(_id)`, `disease_group(_id)` | derived | from the pinned UBERON / MONDO graphs and our anchor lists (§7.6) |

**Technology**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `technology` | recommended → required after reprocessing | a `technology.json` id, by instrument (§7.3) |
| `vendor` | recommended → required after reprocessing | that technology's `vendor_id` |
| `variant` | optional | one of its `variants` |
| `assay` | recommended → required after reprocessing | its first `legacy_assay` (`MERFISH`, `Vizgen MERSCOPE`, `10x Xenium`) |
| `panel_size` | recommended | genes measured (`stats.json` `genes`) |
| `product` | recommended | the vendor's product name as SpatialData-DB writes it (`MERFISH`, `MERSCOPE`, `In Situ Gene Expression`) |
| `modality`, `modalities` | recommended | `RNA`, `["RNA"]` for every dataset we serve today |
| `chemistry_version`, `staining_method`, `preservation_method`, `magnification` | optional | when the source says |

**Curation**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `validation_status` | required | `pass` once bionty passes every ontology column (M3). Their release holds only `pass` rows; so does ours |
| `tier` | optional | null. SpatialData-DB grades its own conversions `bronze` (8,346 rows) or `silver` (4,214) by an unpublished rule; it says nothing about our data |
| `review_status` | derived | `under_review` while a person's review is open; the row stays listed (M3) |
| `review_reason` | optional | the reviewer's text |
| `pathological` | derived | false for `normal`, true for a MONDO disease, null for `unknown` |
| `low_obs` | derived | `n_obs < 100` |

**Licence**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `license` | recommended | an id of `conventions.json` `metadata.licenses` (BIL's are SPDX ids already) |
| `license_url`, `license_noncommercial`, `license_redistributable`, `license_unknown` | derived | from that table |

**Counts and the store**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `n_obs` | required | rows of the default table (`stats.json` `cells`) |
| `n_features` | required | its columns (`stats.json` `genes`, blanks and controls filtered) |
| `total_counts` | recommended | sum of X (`qc.json` `total_counts`) |
| `n_transcripts` | recommended | molecules in the points elements, all sections (`stats.json` `molecules`) |
| `non_integer_counts` | optional | from X; the pipeline knows |
| `default_table` | recommended | `table` |
| `obs_unit`, `observation_unit`, `observation_unit_source` | recommended | `cells`, `cell`, `technology_map` (our tables do not record `uns["spatialdata_io"]`) |
| `unit_size_um` | optional | null for cells |
| `segmentation_source` | recommended | `vendor` when the cells came from the instrument's output (`stats.json` `input` `merscope` / `xenium`), `depositor` when from the depositor's table |
| `segmentation_method` | optional | |
| `hne_image`, `if_image` | recommended | from the store's image elements (DAPI, PolyT, morphology → `if_image`) |
| `ftu_annotation`, `mpp`, `image_width`, `image_height` | optional | from the store |
| `sdata_repr` | optional | `print(sdata)`'s text, by our formatter (§9 X9) |

**Access**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `zarr_url` | required | `s3://<bucket>/<key>.zarr`, as they write it; the viewer gets the HTTPS form from `meyes_viewer_url` |
| `card_url` | recommended | the artifact's page in eyes-catalog |
| `size_bytes` | recommended | the store's bytes from the S3 listing (`register_lamin.py` checks it) |
| `published_at` | required | when the row first entered the index |
| `slug` | derived | their pattern: `slug(assay)-slug(tissue)-lower(uid)` |
| `title` | required | the dataset's title |
| `vitessce_url` | n/a | null |

**Publication and sample**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `publication_date` | recommended | the date the data was published at its source: BIL's submission date (`bildate`; we keep only its year today) |
| `doi` | recommended | bare DOI of the paper |
| `pmid`, `disease_details`, `sample_title`, `donor_id`, `treatment`, `condition` | optional | when stated; `donor_id` per section in `meyes_sections` when sections differ |
| `publication_year` | recommended | the paper's year — not BIL's `publicationYear`, which is the submission's |
| `publication_verify_status` | recommended | `provisional` when a DOI is recorded and unchecked; `verified` once checked; `no_associated_publication` only when a person says so |
| `sample_id` | recommended | the BIL code; else the source's id |
| `dataset_url` | recommended | the source's page (BIL dataset page; `externalLink`) |
| `timepoint`, `timepoint_kind`, `timepoint_days` | recommended | §7.4 |
| `genotype` | recommended | free text, as SpatialData-DB (`C57BL/6J`; BIL's `Strain_Name='…'` unwrapped) |

**Re-releases**

| Column | Level | A MERFISHEYES row |
|---|---|---|
| `same_section_group` | derived | rows holding the same matrix (§7.7): the newest member's `uid` (its `canonical_key` until uids exist); null outside a group |
| `same_section_relation` | derived | `rerelease` within a group, as SpatialData-DB writes it |

**SpatialData-DB's own bookkeeping** (n/a, null): `sddb_id`, `database_version`,
`capture_area` (Visium HD), `same_section_release` (Space Ranger release),
`same_section_panel_of` (10x).

Their release also carries seven columns not in the schema (`n_cells`,
`median_transcripts_per_cell`, `median_genes_per_cell`, `decoded_transcripts`,
`imaged_area`, `onboard_version`, `panel`: Atera, 4–19 rows). We do not write
them.

### 3.2 Our columns (`meyes_*`, 26)

| Column | Level | Meaning |
|---|---|---|
| `meyes_id` | required | the `CatalogDataset` id |
| `meyes_source` | required | Explore's category: `curated`, `bil`, `community` |
| `meyes_kind` | required | `single_cell`, `single_molecule`, `image` (`lib/datasets/kind.ts`) |
| `meyes_page_url` | required | `https://www.merfisheyes.com/explore/<meyes_id>` |
| `meyes_viewer_url` | required | opens the store in the viewer |
| `meyes_description` | recommended | the card's text (SpatialData-DB rows have none) |
| `meyes_sex`, `meyes_sex_id` | recommended | PATO: `male` 0000384, `female` 0000383, `mixed sex` 0001338 |
| `meyes_development_stage_id` | recommended | the HsapDv / MmusDv id of `development_stage` |
| `meyes_n_sections` | recommended | sections in the store |
| `meyes_sections` | recommended | one struct per section: `id` (its `_sample_id`), `sample_id`, `donor_id`, `timepoint`, `timepoint_days`, `meyes_sex`, `n_obs`, `n_transcripts`, `meyes_matrix_hash` (M6) |
| `meyes_matrix_hash` | recommended | `mh1:` + SHA-256 of the default table's counts, the same however the matrix is stored or ordered (§7.7); equal hashes = the same data |
| `meyes_dims` | recommended | 2 or 3 |
| `meyes_de_columns` | recommended | obs columns with DE statistics (`stats.json` `de_stats`) |
| `meyes_thumbnail_url` | recommended | the card's image |
| `meyes_pipeline` | recommended | git short sha that wrote the store; null for a store we did not write |
| `meyes_pipeline_input` | optional | what `prepare.py` recognised the input as |
| `meyes_bil_code` | optional | |
| `meyes_investigator`, `meyes_co_investigators`, `meyes_institution`, `meyes_funding` | optional | BIL's contributor and funder records |
| `meyes_tags` | optional | free-text tags |
| `meyes_n_molecule_genes` | optional | genes in the points elements (may include blanks) |
| `meyes_cell_type_method`, `meyes_cell_type_levels` | optional | e.g. `MapMyCells hierarchical`; the taxonomy levels present as obs columns |

QC numbers (`pct_artifact`, `mask_percentile`, bootstrapping) stay in `qc.json`
and the Lamin record; the index carries what a filter needs.

## 4. In the store

The root `zarr.json` `attributes` of a store we write, beside `spatialdata_attrs`
(the reprocessed `ace-dip-van`; the test's example):

```json
{
  "metadata_general": {
    "Product": "MERSCOPE", "Assay": "Vizgen MERSCOPE",
    "Biomaterial Type": "Specimen from Organism",
    "Organism": "mouse", "Tissue": "brain", "Modality": "RNA", "Disease": "unknown",
    "Dataset Url": "https://api.brainimagelibrary.org/web/view?bildid=ace-dip-van",
    "Sample ID at Source": "ace-dip-van", "Default Table": "table",
    "License": "CC BY-NC 4.0", "Development Stage": "13-week-old stage",
    "HnE Image Available": false, "IF Image Available": false,
    "FTU Annotation Available": false, "Number of Transcripts": 472686422,
    "General Metadata Version": "0.3.2"
  },
  "metadata_merscope": { "Merscope Metadata Version": "0.1.0" },
  "meyes_metadata": {
    "version": "1",
    "canonical_key": "bil:ace-dip-van", "study_id": "bil:ace-dip-van",
    "organism": "mouse", "organism_id": "NCBITaxon:10090",
    "tissue": "brain", "tissue_id": "UBERON:0000955",
    "disease": "unknown", "disease_id": null,
    "development_stage": "13-week-old stage", "meyes_development_stage_id": "MmusDv:0000159",
    "technology": "vizgen_merscope", "vendor": "vizgen", "assay": "Vizgen MERSCOPE",
    "timepoint": "91.0 Days", "timepoint_kind": "age", "timepoint_days": 91.0,
    "genotype": "C57BL/6J", "meyes_sex": "male", "meyes_sex_id": "PATO:0000384",
    "license": "CC-BY-NC-4.0", "n_obs": 243224, "n_features": 300,
    "meyes_matrix_hash": "mh1:d0ff37d6072dd0d8ef11bae317b5fa6e573f59364d40c61962dd183432dfb025",
    "…": "…"
  },
  "meyes_provenance": {
    "canonical_key": "bil:ace-dip-van", "pipeline": "93ab1c0", "input": "merscope",
    "written_at": "…",
    "sources": { "organism": "bil-api", "technology": "store", "…": "…" },
    "validation": { "status": "pass", "at": "…", "bionty": "…",
                    "sources": { "organism": "ncbitaxon 2023-06-20", "tissue": "uberon 2024-08-07" } }
  }
}
```

- **`metadata_general`** is written as SpatialData-DB writes it: their keys,
  Title Case, **labels not ids** (their stores carry no ontology ids at all —
  the ids live in their Lamin and their catalog). Of their 25 keys we write
  every one we know, and not `Database Version` or `Count Matrix Hash`
  (theirs; §9 X9 — our hash is `meyes_matrix_hash`). `License` is the licence's label (`CC BY 4.0`), as in their
  stores, not its id.
- **`metadata_<technology>`**: the block their converter writes for the
  technology (`conventions.json` `technologyBlocks`). For MERFISH and
  MERSCOPE it holds only its version.
- **`meyes_metadata`** (decided 2026-11-16, §9 X1): every index column that is
  a fact about the data, under the index's column names — so the ids travel
  with the store and a store opened from disk can be filtered by them. Its
  `version` is this profile's.
- **`meyes_provenance`** mirrors their `sddb_provenance`: what wrote the store,
  where each value came from (`bil-api`, `bil-csv`, `owner`, `curator`, `store`,
  `agent` confirmed by a person, `sddb`), what bionty checked it against, and
  `matrix_hash`: how `meyes_matrix_hash` was computed (version, the matrix and
  gene names hashed, cells, non-zeros) and its 64-entry sketch (§7.7).
- We never write `sddb_provenance` or `sddb`: those say SpatialData-DB
  converted the store.

`metadata.schema.json` is the schema of these four; `validate_store.py` will
check it (step 3).

## 5. In Lamin

What is in eyes-catalog today is not the target and is not migrated: once
the pipeline writes this metadata, its sync (the scripts on
`bil/pipeline-dry-run`, updated for it) registers every reprocessed store
afresh. The records then use SpatialData-DB's feature names, which are the
`metadata_general` keys, typed to bionty registries where theirs are:

| Feature | Type (theirs) |
|---|---|
| `Organism` | `cat[bionty.Organism]` |
| `Tissue` | `cat[bionty.Tissue]` |
| `Disease` | `cat[bionty.Disease]` |
| `Development Stage` | `cat[bionty.DevelopmentalStage]` |
| `Assay` | `cat[bionty.ExperimentalFactor]` (the technology's `efo_id`: MERFISH EFO:0008992, MERSCOPE EFO:0022956, Xenium EFO:0022615) |
| `Product`, `Modality`, `License`, `Biomaterial Type`, `Preservation Method`, … | `cat[ULabel]` |
| `Genotype`, `Donor ID`, `Timepoint`, `Timepoint Kind`, `Condition`, `Treatment`, `Sample ID at Source`, `Study ID` | `str` |
| `Timepoint Days` | `float` |

Our own features keep `meyes_*` names. They have no sex feature; ours is
`meyes_sex` (`cat[bionty.Phenotype]`). These replace eyes-catalog's
free-text, snake_case features (`species`, `tissue`, `platform`), kept free
text by LAMIN.md's decision of 2026-11-09, which M1 and M3 supersede.

## 6. From today's fields

**Postgres** (`CatalogDataset`, and `Dataset` for private rows):

| Field | Column(s) |
|---|---|
| `title` | `title` |
| `description` | `meyes_description` |
| `species` | `organism`, `organism_id` (§7.1) |
| `tissue` | `tissue`, `tissue_id` (§7.2) |
| `disease` | `disease`, `disease_id` |
| `platform` | `technology`, `vendor`, `assay`, `product`, `modality(ies)` (§7.3) |
| `institute` | `meyes_institution` |
| `tags` | `meyes_tags` |
| `genes` | the gene index (step 7), not a column |
| `numCells`, `numGenes` | `n_obs`, `n_features` — the store's values win |
| `thumbnailUrl` | `meyes_thumbnail_url` |
| `bilCode` | `meyes_bil_code`, `canonical_key`, `sample_id` |
| `externalLink` | `dataset_url` |
| `publicationLink` | `doi` when it holds one (doi.org, bioRxiv URLs) |
| `license`, `licenseNoncommercial` | `license` (+ derived flags) |
| `source` | `meyes_source` |
| `entries[]` | `meyes_kind`, `meyes_viewer_url`, `zarr_url` (the store entry) |
| `metadata.age` | `timepoint`, `timepoint_kind`, `timepoint_days`, `development_stage` (§7.4) |
| `metadata.sex` | `meyes_sex`, `meyes_sex_id` |
| `metadata.genotype` | `genotype` |
| `metadata.license` | `license` |
| `metadata.publicationYear` (BIL) | the year of `publication_date` (the submission), not `publication_year` |
| `metadata.investigator`, `coInvestigators`, `institution`, `funding` | `meyes_investigator`, `meyes_co_investigators`, `meyes_institution`, `meyes_funding` |
| `metadata.technique` | evidence for `technology` only |
| `isPublished`, `isFeatured`, `isInternal`, `isCommunity`, `reviewStatus`, `reviewNote`, `viewCount`, `sortOrder`, ownership | Postgres only |
| `sddbUid`, `sourceRelease` (SpatialData-DB rows) | theirs: their rows enter the index as their rows (step 6), not remapped |

**Beside the store** (BIL pipeline):

| File | Field | Column(s) |
|---|---|---|
| `catalog.json` | as the Postgres fields above (it is one `bil-completed.json` item) | |
| | `source` | `meyes_provenance.sources` |
| `stats.json` | `cells`, `genes`, `molecules`, `molecule_genes`, `dims`, `de_stats` | `n_obs`, `n_features` and `panel_size`, `n_transcripts`, `meyes_n_molecule_genes`, `meyes_dims`, `meyes_de_columns` |
| | `images` | `hne_image`, `if_image` |
| | `input`, `pipeline`, `packages` | `meyes_pipeline_input` (and `technology`, `segmentation_source` when it names an instrument), `meyes_pipeline`, `meyes_provenance.packages` |
| `qc.json` | `n_samples` | `meyes_n_sections` |
| | `total_counts` | `total_counts` |
| | `mmc_methods`, `mmc_levels` | `meyes_cell_type_method`, `meyes_cell_type_levels` |
| | the rest | Lamin and `qc.json` only |
| `mapping.json` | — | not metadata (what the viewer draws beside the cells) |

**Not kept today, and available**: the BIL API's full submission date (`bildate`),
and per section, the sections' names and cell counts (the store's `_sample_id`).

## 7. Rules

### 7.1 Organism

Our spellings → SpatialData-DB's name (Ensembl release-112, lower case) and id:

| Ours | `organism` | `organism_id` |
|---|---|---|
| Mouse, Mice | mouse | NCBITaxon:10090 |
| Human | human | NCBITaxon:9606 |
| Rhesus Macaque, Macaque | macaque | NCBITaxon:9544 |
| Marmoset | white-tufted-ear marmoset | NCBITaxon:9483 |
| Rat | rat | NCBITaxon:10116 |
| Drosophila melanogaster, Schmidtea mediterranea | bionty's NCBITaxon name | (step 2) |

The table replaces `conventions.json` `catalog.species` when Explore reads the
index (step 6); until then both exist.

### 7.2 Tissue

The most specific UBERON term the source states. Today: Brain → `brain`
(UBERON:0000955), Brainstem → `brainstem` (UBERON:0002298), Breast → `breast`
(UBERON:0000310), Skin → `skin of body` (UBERON:0002097, SpatialData-DB's
choice), Whole embryo → `embryo` (UBERON:0000922). "Whole body", "Whole
organism" and "Whole embryo, heart" (two tissues; a row has one) go to a person.
"Brain" is right but coarse for most BIL datasets; a region from the title
("prefrontal cortex", "hypothalamus") is a backfill proposal, not a rule.

### 7.3 Technology

By the **instrument**, as SpatialData-DB does: it lists the Zhuang lab's
whole-brain data as `merfish` (academic) and Yao et al. 2023's as
`vizgen_merscope`. So:

- `stats.json` `input` `merscope` → `vizgen_merscope`; `xenium` → `10x_xenium`.
- Platform spellings that name one technology: Xenium → `10x_xenium`,
  MERSCOPE → `vizgen_merscope`, Stereo-seq → `stomics_stereo_seq`.
- **"MERFISH" alone does not decide** (all 92 BIL rows; the community row) nor
  does "MERFISH+": a person picks `merfish` or `vizgen_merscope` from the BIL
  record, the title ("MERSCOPE v1: 609889") or the paper. Step 5.

`vendor`, `assay` and the EFO id follow from `technology.json`.

### 7.4 Age and development stage

From the source's text (`28.0 Days`, `50.0 Years`, `20.0 Gestational weeks`):

- `timepoint` = the text as given.
- Postnatal: `timepoint_kind` `age`; `timepoint_days` = days since birth
  (months × 30.4375, years × 365.25, as SpatialData-DB's `4 months` = 121.75).
- Prenatal: `timepoint_kind` `stage`; `timepoint_days` = days
  post-fertilization. Embryonic days are already that (E16.5 → 16.5, as
  theirs). **Gestational weeks count from the last menstrual period**, two
  weeks earlier: GW20 → 126 days, development stage `18th week
  post-fertilization stage` (HsapDv). §9 X7.
- `development_stage`: the most specific HsapDv / MmusDv term for the age
  (mouse 28 days → `4-week-old stage` MmusDv:0000149; 91 days →
  `13-week-old stage` MmusDv:0000159; human 50 years → `50-year-old stage`
  HsapDv:0000144). No such ontology for macaque or marmoset: null, the age is
  in `timepoint`.

### 7.5 The rest

- **Sex**: `meyes_sex` from BIL's `Male` / `Female`; a store whose sections
  differ is `mixed sex`, each section's in `meyes_sections`.
- **Genotype**: BIL's `Strain_Name='C57Bl/6J'` → `C57BL/6J` (unwrapped, the
  strain's usual spelling); free text otherwise.
- **Licence**: the id is the index's; the label goes in the store (§4).
- **Sections** (M6, kept 2026-11-16): one store, one row; its sections in
  `meyes_sections`, in the order of the store's `_sample_id` categories.
  `sample_id`, `donor_id`, `timepoint` on the row hold the value when every
  section shares it, else null with the values per section.

### 7.6 Tissue and disease groupings

The four facet columns (`tissue_system`, `tissue_organ`, `disease_class`,
`disease_group`) let "brain" match hippocampus by equality. SpatialData-DB
computes them from anchor tables it does not publish; ours are our own
convention, `conventions.json` `metadata.groupings` (decided 2026-11-16, §9
X11):

- **The anchors** are ordered `[id, label]` lists: 12 tissue systems, 40
  organs, 8 disease classes, 39 disease groups — the terms SpatialData-DB's
  `catalog-2026-11-12` actually uses, with their labels (`neoplasm, other`,
  `cancer, other site` are theirs). So where both catalogs have a term they
  group it alike, and a filter means the same over both.
- **The rule**: take the term and all its ancestors in the pinned UBERON /
  MONDO graph (is_a and part_of); the group is the **first anchor in list
  order** that is among them. The lists run specific before general, so
  order decides: brain cancer is in `cancer`, not `neurological`; retina is
  organ `retina`, not `eye`, and system `sensory system`, not `nervous
  system`. No graph distances are needed.
- **No organ anchor** → the organ is the system's term (as theirs:
  `digestive system` appears as an organ). **No anchor at all** → `other
  tissues` / `other conditions`, id null. `unknown` stays `unknown`;
  `normal` has no class or group.
- Ours today: every BIL row is brain → organ `brain`, system `nervous system`.
- A query for any other term ("hippocampus and below") uses
  `ontology_rollup.json` (each term → its ancestors present), which the index
  build writes for our terms in SpatialData-DB's format.

Lists grow by adding entries; reordering changes groupings, so a change is a
spec change. If SpatialData-DB publishes their tables, we compare and adopt.

### 7.7 The count-matrix hash and re-releases

Data gets re-released: reprocessed, re-uploaded, published again under
another accession, one section on its own. The hash says when two stores hold
the same matrix whatever their bytes (decided 2026-11-16, M12). Reference:
[`scripts/metadata/matrix_hash.py`](../../scripts/metadata/matrix_hash.py);
`test_matrix_hash_names_the_same_data_the_same` pins it.

**mh1**, over the default table's counts (`X`, or a counts layer when `X` is
normalised — recorded in `meyes_provenance.matrix_hash.matrix`):

1. Gene names stripped and lower-cased; they must be unique.
2. Each value rounded half-even to 6 decimals; integers written without a
   fraction (`3`), others as Python's shortest repr (`0.25`); zeros left out.
3. Each cell: its non-zero `gene\tvalue\n` lines sorted by gene; the cell's
   digest is the first 16 bytes of their SHA-256 (a cell with no counts hashes
   the empty string).
4. The digests sorted, after the bytes `mh1\n`, SHA-256 → `mh1:<hex>`.

So storage (CSR, CSC, dense, chunks, compression), dtype width, cell order,
gene order and genes with no counts do not change it; any count does. Cell
names, annotations and coordinates are not in it. Each section gets its own
(`meyes_sections[].meyes_matrix_hash`), so a section re-released alone is
found too. Measured on `ace-dip-van` (243,224 cells, 39.9 M non-zeros, 3
sections): 38 s and 1.7 GB in plain Python — a few minutes for a 4 M-cell
store; the pipeline computes it while the matrix is in memory.

**Near matches** — a re-release with some cells filtered, or renormalised:

- **The sketch**: the 64 smallest cell digests, kept in `meyes_provenance`.
  Two sketches estimate the share of cells two stores have in common
  (`sketch_overlap`), for free from the digests above.
- **Per-gene totals**: the gene index (step 7) holds each dataset's total
  counts and fraction of cells detected per gene — for SpatialData-DB's
  datasets too, from their published `gene_index.parquet`. Near-equal vectors
  on the shared genes flag a likely re-release across both catalogs without
  opening a store.

**How a match is used** (the index build, step 6): rows with equal
`meyes_matrix_hash` form a `same_section_group` with relation `rerelease`;
a near match is a proposal a person confirms before it groups. Search can
then show one row per group. SpatialData-DB's `Count Matrix Hash` (`cmh1:`)
uses an unpublished algorithm; if they publish it we compute it as well and
match their rows exactly.

## 8. What today's rows can fill

From [`metadata-inventory.md`](metadata-inventory.md): the interim 16
required and 45 recommended columns counted; "filled" is a value present or
made by a fixed rule from one.

| Rows | Count | Required filled (median) | Recommended filled (median) | Once they open their reprocessed store |
|---|---|---|---|---|
| BIL | 92 | 13 of 16 | 20 of 45 | 59 have one: 16 of 16, 32 recommended |
| curated | 56 (47 of them BIL Xenium, `ace-owl-*`) | 13 | 20 | 3 have one: 16, 35 |
| community | 1 | 12 | 11 | — |

- **Every row has organism and tissue** that a rule maps (except five curated
  tissues for a person). Title, description and dataset URL are there.
- **What stops a row today is the store**: `zarr_url`, `n_obs` and
  `meyes_viewer_url` come with it. The 59 BIL rows with a reprocessed store
  fill all 16 required columns; the rest wait for reprocessing.
- **Technology** (now recommended) is derived for the 6 stores made from
  MERSCOPE output and is a person's call for the 56 made from a folder of
  parts (`input: folder`) and for every row without a store: "MERFISH" does
  not name the instrument.
- **Age, sex, genotype, licence** are on 84–87 of the 92 BIL rows (genotype on
  68): the recommended set's best-covered part.
- **Missing almost everywhere**: disease (one curated row has it), DOI (three
  rows; two more link a paper without one), donor, condition.

## 9. Departures from SpatialData-DB

| # | Departure | Why |
|---|---|---|
| X1 | Ontology ids are also in the store, in `meyes_metadata` beside `metadata_general` (decided 2026-11-16) | M2: the store is the copy that travels; theirs carry labels only |
| X2 | One row per multi-section store, sections in `meyes_sections` (M6, kept 2026-11-16) | Their Zhuang ABCA rows are one per section, grouped by `study_id`; ours are viewed as one store |
| X3 | A row needs a store; legacy chunked datasets stay in Postgres until reprocessed (decided 2026-11-16) | `zarr_url` is required |
| X4 | Licence ids `CC-BY-NC-4.0` and `CC-BY-SA-4.0` added to their six | BIL uses them (21 and 4 rows) |
| X5 | Sex: `meyes_sex`, `meyes_sex_id` (PATO) | They have no sex column or feature; BIL records it for 84 of 92 |
| X6 | `meyes_development_stage_id` | Their catalog publishes the stage's name only (their Lamin feature is typed) |
| X7 | Prenatal `timepoint_days` is days post-fertilization, `timepoint_kind` `stage`; gestational weeks shifted by two | Their rows disagree (`E16.5` stage 16.5; `9 weeks of gestation` age 63; `GW_18` age null) |
| X8 | `vitessce_url`, `sddb_id`, `database_version`, `same_section_release`, `same_section_panel_of`, `capture_area` null; `same_section_group` / `relation` filled from our matrix hash, not from 10x library records; `card_url` is our Lamin page | Theirs, or not our data |
| X9 | `metadata_general` without `Database Version` and `Count Matrix Hash`; `SpatialData Repr` from our formatter | Their database's; their hash's algorithm (`cmh1:`) is unpublished, so ours is `meyes_matrix_hash` (`mh1:`, §7.7) in `meyes_metadata`; spatialdata 0.8 cannot read our points element |
| X10 | `meyes_description` | Their rows have no description column |
| X11 | Our own grouping anchors (§7.6) and `tier` null | Their anchor tables and tier rule are unpublished; our anchors are the terms their release uses, so groupings agree where both have a term |
| X12 | Interim required set (M10): `uid`, `technology`, `vendor`, `assay` recommended until the fleet is reprocessed; ids required throughout | What search and the MCP need now is organism, tissue and a store; the instrument is a person's call for most rows today. Their validator still passes the combined index (it refuses only an all-null required column) |

## 10. Open questions

1. **The instrument of each BIL dataset** (§7.3) — part of reprocessing
   (step 5); needed before `technology` is required again.
2. **Version drift.** Their stores carry `General Metadata Version` 0.3.0–0.3.2
   and their release moves near-nightly; we follow by bumping the pins in
   `conventions.json` after re-reading `_catalog_schema.py`, as
   `lib/sddb/catalog-schema.ts` does.
3. Whether their `metadata_general` keys become an scverse convention (plan §5).
4. Nice to have from Tim, not blocking: their anchor tables (to compare with
   §7.6) and the `cmh1:` hash algorithm (to spot one matrix in both
   catalogs).
