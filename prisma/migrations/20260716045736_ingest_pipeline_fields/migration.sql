-- AlterEnum
ALTER TYPE "DatasetStatus" ADD VALUE 'QUEUED';

-- AlterTable
ALTER TABLE "datasets" ADD COLUMN     "batch_job_id" TEXT,
ADD COLUMN     "ingest_source" TEXT,
ADD COLUMN     "owner_id" TEXT,
ADD COLUMN     "processing_params" JSONB,
ADD COLUMN     "processing_progress" JSONB;

-- AlterTable
ALTER TABLE "users" ADD COLUMN     "compute_tier" TEXT NOT NULL DEFAULT 'standard';

-- CreateIndex
CREATE INDEX "datasets_owner_id_idx" ON "datasets"("owner_id");

-- AddForeignKey
ALTER TABLE "datasets" ADD CONSTRAINT "datasets_owner_id_fkey" FOREIGN KEY ("owner_id") REFERENCES "users"("id") ON DELETE SET NULL ON UPDATE CASCADE;
