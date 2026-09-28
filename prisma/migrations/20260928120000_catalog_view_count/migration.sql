-- Denormalized view count on catalog rows, so Explore can order by popularity
-- without joining through catalog_dataset_entries -> datasets on every listing.
ALTER TABLE "catalog_datasets"
  ADD COLUMN "view_count" INTEGER NOT NULL DEFAULT 0;

CREATE INDEX "catalog_datasets_view_count_idx"
  ON "catalog_datasets"("view_count");

-- Backfill from the existing per-dataset counts. An entry resolves to a
-- dataset either directly (dataset_id, app-uploaded rows) or by matching
-- s3_base_url (curated rows, which is the large majority). EXISTS rather than
-- a join so a dataset reachable through several of a row's entries is still
-- counted once -- matching how the live increment behaves.
UPDATE "catalog_datasets" c
   SET "view_count" = COALESCE((
         SELECT SUM(d."view_count")
           FROM "datasets" d
          WHERE EXISTS (
                SELECT 1
                  FROM "catalog_dataset_entries" e
                 WHERE e."catalog_id" = c."id"
                   AND (e."dataset_id" = d."id" OR e."s3_base_url" = d."s3_base_url")
          )
       ), 0);
