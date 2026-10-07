-- Drop the unique constraint on datasets.fingerprint.
--
-- A failed or abandoned upload leaves its Dataset row in place. The duplicate
-- pre-check only looks at COMPLETE/PROCESSING rows, so it reported "safe to
-- upload", and then the insert hit this constraint and returned P2002
-- ("A dataset with this fingerprint already exists"). Re-uploading a file that
-- had previously failed was therefore impossible without manual cleanup.
--
-- Fingerprints are still computed and stored. The non-unique index below keeps
-- fingerprint lookups fast.
DROP INDEX IF EXISTS "datasets_fingerprint_key";
CREATE INDEX IF NOT EXISTS "datasets_fingerprint_idx" ON "datasets"("fingerprint");
