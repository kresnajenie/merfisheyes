-- Gene symbols per uploaded dataset, so a user's own uploads are searchable
-- by gene (and community submissions carry their genes into the catalog).
ALTER TABLE "datasets" ADD COLUMN "genes" TEXT[] NOT NULL DEFAULT '{}';

-- Same substring-search index pattern as catalog_datasets (migration
-- 20260803120000, which created the IMMUTABLE catalog_genes_text() helper —
-- generic text[] -> text, reused here).
CREATE INDEX IF NOT EXISTS "datasets_genes_trgm_idx"
  ON "datasets"
  USING gin (catalog_genes_text("genes") gin_trgm_ops);
