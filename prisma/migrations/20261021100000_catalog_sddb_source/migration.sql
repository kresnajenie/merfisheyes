-- Where a catalog row came from, and what a SpatialData-DB row carries.
--
-- Explore has grown three kinds of curated row (admin-authored, BIL import,
-- community submission) told apart by flags, and now a fourth: rows synced
-- from a SpatialData-DB catalog release (scripts/sync-sddb-catalog.ts,
-- docs/specs/sddb-catalog.md). `source` is the one column to filter on; the
-- flags stay as they are. The sddb_* columns are the sync's upsert key, the
-- licence shown on every card, and the release a row was last written from.
ALTER TABLE "catalog_datasets"
  ADD COLUMN "source" VARCHAR(20) NOT NULL DEFAULT 'curated',
  ADD COLUMN "sddb_uid" VARCHAR(40),
  ADD COLUMN "license" VARCHAR(100),
  ADD COLUMN "license_noncommercial" BOOLEAN,
  ADD COLUMN "source_release" VARCHAR(100),
  ADD COLUMN "source_generated_at" VARCHAR(40);

CREATE UNIQUE INDEX "catalog_datasets_sddb_uid_key" ON "catalog_datasets"("sddb_uid");
CREATE INDEX "catalog_datasets_source_idx" ON "catalog_datasets"("source");

-- Backfill from the flags that identified the existing kinds.
UPDATE "catalog_datasets" SET "source" = 'bil' WHERE "is_bil";
UPDATE "catalog_datasets" SET "source" = 'community' WHERE "is_community";
