-- CreateTable
CREATE TABLE "catalog_dataset_entries" (
    "id" TEXT NOT NULL,
    "catalog_id" TEXT NOT NULL,
    "label" VARCHAR(200) NOT NULL,
    "dataset_type" VARCHAR(50) NOT NULL,
    "s3_base_url" VARCHAR(1000),
    "dataset_id" VARCHAR(50),
    "sort_order" INTEGER NOT NULL DEFAULT 0,

    CONSTRAINT "catalog_dataset_entries_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE INDEX "catalog_dataset_entries_catalog_id_idx" ON "catalog_dataset_entries"("catalog_id");

-- AddForeignKey
ALTER TABLE "catalog_dataset_entries" ADD CONSTRAINT "catalog_dataset_entries_catalog_id_fkey" FOREIGN KEY ("catalog_id") REFERENCES "catalog_datasets"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- Migrate existing data: create entries from existing catalog_datasets rows
INSERT INTO "catalog_dataset_entries" ("id", "catalog_id", "label", "dataset_type", "s3_base_url", "dataset_id", "sort_order")
SELECT
    gen_random_uuid()::text,
    "id",
    CASE
        WHEN "dataset_type" = 'single_cell' THEN 'Single Cell'
        WHEN "dataset_type" = 'single_molecule' THEN 'Single Molecule'
        ELSE 'Dataset'
    END,
    "dataset_type",
    "s3_base_url",
    "dataset_id",
    0
FROM "catalog_datasets"
WHERE "dataset_type" IS NOT NULL
  AND ("s3_base_url" IS NOT NULL OR "dataset_id" IS NOT NULL);

-- Drop legacy columns
ALTER TABLE "catalog_datasets" DROP COLUMN "dataset_type";
ALTER TABLE "catalog_datasets" DROP COLUMN "s3_base_url";
ALTER TABLE "catalog_datasets" DROP COLUMN "dataset_id";
