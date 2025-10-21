-- AlterTable
ALTER TABLE "datasets" ADD COLUMN     "dataset_type" VARCHAR(50),
ADD COLUMN     "manifest_json" JSONB;

-- CreateIndex
CREATE INDEX "datasets_dataset_type_idx" ON "datasets"("dataset_type");
