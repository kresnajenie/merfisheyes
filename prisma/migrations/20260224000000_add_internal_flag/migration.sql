-- AlterTable
ALTER TABLE "catalog_datasets" ADD COLUMN "is_internal" BOOLEAN NOT NULL DEFAULT false;

-- CreateIndex
CREATE INDEX "catalog_datasets_is_internal_idx" ON "catalog_datasets"("is_internal");
