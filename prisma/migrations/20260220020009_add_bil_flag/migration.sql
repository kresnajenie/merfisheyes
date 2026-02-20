-- AlterTable
ALTER TABLE "catalog_datasets" ADD COLUMN     "is_bil" BOOLEAN NOT NULL DEFAULT false;

-- CreateIndex
CREATE INDEX "catalog_datasets_is_bil_idx" ON "catalog_datasets"("is_bil");
