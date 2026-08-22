-- Compute tier the current/last processing job ran on, and how many jobs have
-- been submitted — used by admin retry + automatic OOM escalation.
ALTER TABLE "datasets" ADD COLUMN "compute_tier" VARCHAR(20);
ALTER TABLE "datasets" ADD COLUMN "processing_attempts" INTEGER NOT NULL DEFAULT 0;
