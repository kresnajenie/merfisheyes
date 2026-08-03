-- AlterTable
ALTER TABLE "datasets" ADD COLUMN     "s3_base_url" VARCHAR(1000),
ADD COLUMN     "view_count" INTEGER NOT NULL DEFAULT 0,
ADD COLUMN     "viewer_config" JSONB;

-- CreateTable
CREATE TABLE "dataset_views" (
    "id" TEXT NOT NULL,
    "dataset_id" VARCHAR(50) NOT NULL,
    "session_key" TEXT NOT NULL,
    "day" DATE NOT NULL,
    "user_id" TEXT,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "dataset_views_pkey" PRIMARY KEY ("id")
);

-- CreateTable
CREATE TABLE "dataset_views_daily" (
    "dataset_id" VARCHAR(50) NOT NULL,
    "day" DATE NOT NULL,
    "count" INTEGER NOT NULL DEFAULT 0,

    CONSTRAINT "dataset_views_daily_pkey" PRIMARY KEY ("dataset_id","day")
);

-- CreateIndex
CREATE INDEX "dataset_views_dataset_id_day_idx" ON "dataset_views"("dataset_id", "day");

-- CreateIndex
CREATE UNIQUE INDEX "dataset_views_dataset_id_session_key_day_key" ON "dataset_views"("dataset_id", "session_key", "day");

-- CreateIndex
CREATE INDEX "dataset_views_daily_day_idx" ON "dataset_views_daily"("day");

-- CreateIndex
CREATE UNIQUE INDEX "datasets_s3_base_url_key" ON "datasets"("s3_base_url");

-- CreateIndex
CREATE INDEX "datasets_s3_base_url_idx" ON "datasets"("s3_base_url");

-- AddForeignKey
ALTER TABLE "dataset_views" ADD CONSTRAINT "dataset_views_dataset_id_fkey" FOREIGN KEY ("dataset_id") REFERENCES "datasets"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "dataset_views_daily" ADD CONSTRAINT "dataset_views_daily_dataset_id_fkey" FOREIGN KEY ("dataset_id") REFERENCES "datasets"("id") ON DELETE CASCADE ON UPDATE CASCADE;

