-- AlterTable
ALTER TABLE "users" ADD COLUMN     "permissions" TEXT[] DEFAULT ARRAY[]::TEXT[];

-- CreateTable
CREATE TABLE "technology_preset_versions" (
    "hash" VARCHAR(16) NOT NULL,
    "technology" VARCHAR(200) NOT NULL,
    "config" JSONB NOT NULL,
    "note" VARCHAR(1000),
    "tested_on" VARCHAR(1000),
    "app_version" VARCHAR(50) NOT NULL,
    "created_by_id" TEXT,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "technology_preset_versions_pkey" PRIMARY KEY ("hash")
);

-- CreateTable
CREATE TABLE "technology_preset_publishes" (
    "id" TEXT NOT NULL,
    "technology" VARCHAR(200) NOT NULL,
    "version_hash" VARCHAR(16),
    "published_by_id" TEXT,
    "created_at" TIMESTAMP(3) NOT NULL DEFAULT CURRENT_TIMESTAMP,

    CONSTRAINT "technology_preset_publishes_pkey" PRIMARY KEY ("id")
);

-- CreateIndex
CREATE INDEX "technology_preset_versions_technology_created_at_idx" ON "technology_preset_versions"("technology", "created_at");

-- CreateIndex
CREATE INDEX "technology_preset_publishes_technology_created_at_idx" ON "technology_preset_publishes"("technology", "created_at");

-- AddForeignKey
ALTER TABLE "technology_preset_versions" ADD CONSTRAINT "technology_preset_versions_created_by_id_fkey" FOREIGN KEY ("created_by_id") REFERENCES "users"("id") ON DELETE SET NULL ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "technology_preset_publishes" ADD CONSTRAINT "technology_preset_publishes_version_hash_fkey" FOREIGN KEY ("version_hash") REFERENCES "technology_preset_versions"("hash") ON DELETE RESTRICT ON UPDATE CASCADE;

-- AddForeignKey
ALTER TABLE "technology_preset_publishes" ADD CONSTRAINT "technology_preset_publishes_published_by_id_fkey" FOREIGN KEY ("published_by_id") REFERENCES "users"("id") ON DELETE SET NULL ON UPDATE CASCADE;

