-- Add genes array, metadata JSONB, and bilCode to catalog_datasets
ALTER TABLE "catalog_datasets" ADD COLUMN "genes" TEXT[] DEFAULT ARRAY[]::TEXT[];
ALTER TABLE "catalog_datasets" ADD COLUMN "bil_code" VARCHAR(50);
ALTER TABLE "catalog_datasets" ADD COLUMN "metadata" JSONB DEFAULT '{}';

-- bilCode unique constraint
ALTER TABLE "catalog_datasets" ADD CONSTRAINT "catalog_datasets_bil_code_key" UNIQUE ("bil_code");

-- B-tree index on bil_code for lookups
CREATE INDEX "catalog_datasets_bil_code_idx" ON "catalog_datasets"("bil_code");

-- GIN index on genes array for containment queries (@> ARRAY['gene'])
CREATE INDEX "catalog_datasets_genes_gin_idx" ON "catalog_datasets" USING GIN ("genes");

-- GIN index on metadata JSONB for flexible key-value queries
CREATE INDEX "catalog_datasets_metadata_gin_idx" ON "catalog_datasets" USING GIN ("metadata");
