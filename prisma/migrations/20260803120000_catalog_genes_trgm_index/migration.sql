-- Speed up explore gene search (substring ILIKE over the genes array).
--
-- The search matches `EXISTS (unnest(genes) g WHERE g ILIKE '%term%')`, which
-- a plain GIN(text[]) index can't accelerate. A pg_trgm GIN index over the
-- genes joined into one text string does. `array_to_string` is only STABLE, so
-- it can't be used directly in an index expression — wrap it in an IMMUTABLE
-- function (safe: text[] elements are already text and the separator is fixed).

CREATE EXTENSION IF NOT EXISTS pg_trgm;

CREATE OR REPLACE FUNCTION catalog_genes_text(text[])
  RETURNS text
  LANGUAGE sql
  IMMUTABLE
  PARALLEL SAFE
  AS $$ SELECT array_to_string($1, ' ') $$;

CREATE INDEX IF NOT EXISTS "catalog_datasets_genes_trgm_idx"
  ON "catalog_datasets"
  USING gin (catalog_genes_text("genes") gin_trgm_ops);
