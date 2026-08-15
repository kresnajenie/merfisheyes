-- CreateEnum
CREATE TYPE "CatalogReviewStatus" AS ENUM ('PENDING', 'APPROVED', 'REJECTED');

-- AlterTable: owner-editable catalog-style metadata on uploaded datasets
ALTER TABLE "datasets" ADD COLUMN "description" TEXT;
ALTER TABLE "datasets" ADD COLUMN "species" VARCHAR(200);
ALTER TABLE "datasets" ADD COLUMN "disease" VARCHAR(200);
ALTER TABLE "datasets" ADD COLUMN "tissue" VARCHAR(200);
ALTER TABLE "datasets" ADD COLUMN "platform" VARCHAR(200);
ALTER TABLE "datasets" ADD COLUMN "institute" VARCHAR(500);
ALTER TABLE "datasets" ADD COLUMN "tags" TEXT[] NOT NULL DEFAULT ARRAY[]::TEXT[];
ALTER TABLE "datasets" ADD COLUMN "external_link" VARCHAR(1000);
ALTER TABLE "datasets" ADD COLUMN "publication_link" VARCHAR(1000);

-- AlterTable: community moderation fields on catalog rows
ALTER TABLE "catalog_datasets" ADD COLUMN "is_community" BOOLEAN NOT NULL DEFAULT false;
ALTER TABLE "catalog_datasets" ADD COLUMN "review_status" "CatalogReviewStatus" NOT NULL DEFAULT 'APPROVED';
ALTER TABLE "catalog_datasets" ADD COLUMN "review_note" TEXT;
ALTER TABLE "catalog_datasets" ADD COLUMN "submitted_by" TEXT;
ALTER TABLE "catalog_datasets" ADD COLUMN "source_dataset_id" VARCHAR(50);
ALTER TABLE "catalog_datasets" ADD COLUMN "source_project_id" TEXT;

-- CreateIndex
CREATE INDEX "catalog_datasets_is_community_idx" ON "catalog_datasets"("is_community");
CREATE INDEX "catalog_datasets_review_status_idx" ON "catalog_datasets"("review_status");
CREATE INDEX "catalog_datasets_submitted_by_idx" ON "catalog_datasets"("submitted_by");
CREATE INDEX "catalog_datasets_source_dataset_id_idx" ON "catalog_datasets"("source_dataset_id");
CREATE INDEX "catalog_datasets_source_project_id_idx" ON "catalog_datasets"("source_project_id");

-- CreateTable: user-owned project folders
CREATE TABLE "projects" (
    "id" TEXT NOT NULL,
    "title" VARCHAR(500) NOT NULL,
    "description" TEXT,
    "species" VARCHAR(200),
    "disease" VARCHAR(200),
    "tissue" VARCHAR(200),
    "platform" VARCHAR(200),
    "institute" VARCHAR(500),
    "tags" TEXT[] NOT NULL DEFAULT ARRAY[]::TEXT[],
    "thumbnail_url" VARCHAR(1000),
    "external_link" VARCHAR(1000),
    "publication_link" VARCHAR(1000),
    "owner_id" TEXT NOT NULL,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,
    "updated_at" TIMESTAMP(3) NOT NULL,

    CONSTRAINT "projects_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE INDEX "projects_owner_id_idx" ON "projects"("owner_id");

-- CreateTable: project ↔ dataset membership (many-to-many)
CREATE TABLE "project_datasets" (
    "project_id" TEXT NOT NULL,
    "dataset_id" VARCHAR(50) NOT NULL,
    "sort_order" INTEGER NOT NULL DEFAULT 0,
    "added_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "project_datasets_pkey" PRIMARY KEY ("project_id","dataset_id")
);

-- CreateIndex
CREATE INDEX "project_datasets_dataset_id_idx" ON "project_datasets"("dataset_id");

-- AddForeignKey
ALTER TABLE "projects" ADD CONSTRAINT "projects_owner_id_fkey" FOREIGN KEY ("owner_id") REFERENCES "users"("id") ON DELETE CASCADE ON UPDATE CASCADE;
ALTER TABLE "project_datasets" ADD CONSTRAINT "project_datasets_project_id_fkey" FOREIGN KEY ("project_id") REFERENCES "projects"("id") ON DELETE CASCADE ON UPDATE CASCADE;
ALTER TABLE "project_datasets" ADD CONSTRAINT "project_datasets_dataset_id_fkey" FOREIGN KEY ("dataset_id") REFERENCES "datasets"("id") ON DELETE CASCADE ON UPDATE CASCADE;
