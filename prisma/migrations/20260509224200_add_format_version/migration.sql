-- DropIndex
DROP INDEX "catalog_datasets_genes_gin_idx";

-- DropIndex
DROP INDEX "catalog_datasets_metadata_gin_idx";

-- AlterTable
ALTER TABLE "datasets" ADD COLUMN     "format_version" VARCHAR(20) NOT NULL DEFAULT 'chunked';
