-- CreateEnum
CREATE TYPE "DatasetStatus" AS ENUM ('UPLOADING', 'PROCESSING', 'COMPLETE', 'FAILED');

-- CreateEnum
CREATE TYPE "FileStatus" AS ENUM ('PENDING', 'UPLOADING', 'COMPLETE', 'FAILED');

-- CreateTable
CREATE TABLE "datasets" (
    "id" VARCHAR(50) NOT NULL,
    "fingerprint" VARCHAR(64) NOT NULL,
    "title" VARCHAR(500),
    "num_cells" INTEGER NOT NULL,
    "num_genes" INTEGER NOT NULL,
    "status" "DatasetStatus" NOT NULL DEFAULT 'UPLOADING',
    "manifest_url" VARCHAR(1000),
    "error_message" TEXT,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,
    "completed_at" TIMESTAMP(3),

    CONSTRAINT "datasets_pkey" PRIMARY KEY ("id")
);

-- CreateTable
CREATE TABLE "upload_sessions" (
    "id" VARCHAR(50) NOT NULL,
    "dataset_id" VARCHAR(50) NOT NULL,
    "total_files" INTEGER NOT NULL,
    "completed_files" INTEGER NOT NULL DEFAULT 0,
    "expires_at" TIMESTAMP(3) NOT NULL,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "upload_sessions_pkey" PRIMARY KEY ("id")
);

-- CreateTable
CREATE TABLE "upload_files" (
    "id" UUID NOT NULL,
    "upload_session_id" VARCHAR(50) NOT NULL,
    "file_key" VARCHAR(500) NOT NULL,
    "file_size" BIGINT,
    "status" "FileStatus" NOT NULL DEFAULT 'PENDING',
    "uploaded_at" TIMESTAMP(3),
    "error_message" TEXT,

    CONSTRAINT "upload_files_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE UNIQUE INDEX "datasets_fingerprint_key" ON "datasets"("fingerprint");

-- CreateIndex
CREATE INDEX "datasets_fingerprint_idx" ON "datasets"("fingerprint");

-- CreateIndex
CREATE INDEX "datasets_status_idx" ON "datasets"("status");

-- CreateIndex
CREATE INDEX "datasets_created_at_idx" ON "datasets"("created_at" DESC);

-- CreateIndex
CREATE INDEX "upload_sessions_dataset_id_idx" ON "upload_sessions"("dataset_id");

-- CreateIndex
CREATE INDEX "upload_sessions_expires_at_idx" ON "upload_sessions"("expires_at");

-- CreateIndex
CREATE INDEX "upload_files_upload_session_id_idx" ON "upload_files"("upload_session_id");

-- CreateIndex
CREATE INDEX "upload_files_status_idx" ON "upload_files"("status");

-- CreateIndex
CREATE UNIQUE INDEX "upload_files_upload_session_id_file_key_key" ON "upload_files"("upload_session_id", "file_key");

-- AddForeignKey
ALTER TABLE "upload_sessions" ADD CONSTRAINT "upload_sessions_dataset_id_fkey" FOREIGN KEY ("dataset_id") REFERENCES "datasets"("id") ON DELETE CASCADE ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "upload_files" ADD CONSTRAINT "upload_files_upload_session_id_fkey" FOREIGN KEY ("upload_session_id") REFERENCES "upload_sessions"("id") ON DELETE CASCADE ON UPDATE CASCADE;
