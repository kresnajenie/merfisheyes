-- Owner-set thumbnail for a dataset (screenshot uploaded to the app bucket,
-- served via /api/datasets/[id]/thumbnail).
ALTER TABLE "datasets" ADD COLUMN "thumbnail_url" VARCHAR(1000);
