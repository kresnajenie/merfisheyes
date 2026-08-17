-- Flexible display-metadata bag on projects and datasets (mirrors
-- catalog_datasets.metadata): investigator, institution, age, authors,
-- citation, funding, sex, genotype, technique, publication year, license, etc.
ALTER TABLE "projects" ADD COLUMN "metadata" JSONB DEFAULT '{}';
ALTER TABLE "datasets" ADD COLUMN "metadata" JSONB DEFAULT '{}';
