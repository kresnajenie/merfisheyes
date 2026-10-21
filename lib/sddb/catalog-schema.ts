/**
 * The SpatialData-DB catalog contract, vendored from
 * `theislab/spatialdata-db-client` `src/sddb/_catalog_schema.py` at the
 * release this sync was last checked against (`PINNED_RELEASE`). Their SPEC
 * D61: consumers read the published release, never Lamin's database.
 *
 * The sync fails when a required column is missing and warns when a known
 * column is absent, so a change upstream shows up as a failed run, not as
 * empty cards. Bump `PINNED_RELEASE` after re-checking the file.
 */

export const CATALOG_REPO = "theislab/spatialdata-db-client";
export const CATALOG_ASSET = "catalog.parquet";
export const PINNED_RELEASE = "catalog-2026-10-20";

/** `REQUIRED_COLUMNS` upstream: a release without one of these is unusable. */
export const REQUIRED_COLUMNS = [
  "uid",
  "technology",
  "assay",
  "organism",
  "validation_status",
  "zarr_url",
] as const;

/** The columns a card is built from (a subset of their `CATALOG_COLUMNS`). */
export const READ_COLUMNS = [
  ...REQUIRED_COLUMNS,
  "canonical_key",
  "study_id",
  "tissue",
  "disease",
  "development_stage",
  "tier",
  "license",
  "license_url",
  "license_noncommercial",
  "license_unknown",
  "n_obs",
  "n_features",
  "card_url",
  "size_bytes",
  "title",
  "product",
  "chemistry_version",
  "n_transcripts",
  "capture_area",
  "doi",
  "publication_year",
  "sample_id",
  "default_table",
  "dataset_url",
  "thumbnail_url",
  "collections",
] as const;

export type CatalogColumn = (typeof READ_COLUMNS)[number];

/** One row of `catalog.parquet`, as hyparquet hands it back. */
export type CatalogRow = Record<CatalogColumn, unknown>;
