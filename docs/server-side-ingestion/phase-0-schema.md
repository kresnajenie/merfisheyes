# Phase 0 — Schema & scaffolding

← [Back to design doc](README.md)

## Goal

Extend the Prisma schema so the new ingestion flow can track ownership, the
`QUEUED` state, the Batch job, processing params, and live progress — without
touching any existing behaviour. This is pure app-side, no AWS work, and unblocks
every later phase.

## Depends on

Nothing. This is the foundation.

## Tasks

- [ ] In [`prisma/schema.prisma`](../../prisma/schema.prisma):
  - [ ] `Dataset`: add `ownerId String? @map("owner_id")` + relation
        `owner User? @relation(fields:[ownerId], references:[id], onDelete: SetNull)`.
  - [ ] `Dataset`: add `batchJobId String? @map("batch_job_id")`,
        `processingParams Json? @map("processing_params")`,
        `processingProgress Json? @map("processing_progress")`,
        `ingestSource String? @map("ingest_source")`.
  - [ ] `Dataset`: add `@@index([ownerId])`.
  - [ ] `User`: add `computeTier String @default("standard") @map("compute_tier")`
        and inverse relation `datasets Dataset[]`.
  - [ ] `DatasetStatus` enum: add `QUEUED` between `UPLOADING` and `PROCESSING`.
- [ ] `npx prisma migrate dev --name ingest_pipeline_fields` (creates + applies the migration).
- [ ] `npx prisma generate` (refresh the client).
- [ ] Sanity-check that the migration is picked up by the existing build-time
      migration step ([`scripts/safe-migrate.js`](../../scripts/safe-migrate.js)).

## Files touched

- `prisma/schema.prisma`
- `prisma/migrations/**` (generated)

## Design reference

See [README §6 "Database changes"](README.md) for the exact field list and the
fingerprint decision (worker computes the canonical hash post-processing).

## Verification

- [ ] `npx prisma migrate dev` applies cleanly on the dev DB.
- [ ] `npx prisma studio` shows the new columns; existing rows are untouched
      (all new columns nullable / defaulted).
- [ ] Existing upload/viewer flows still work (nothing reads the new columns yet).
- [ ] `npm run build` succeeds (client regenerated, types compile).

## Notes / risks

- All new columns are **nullable or defaulted** → the migration is additive and
  safe on existing prod data.
- `QUEUED` is inserted into the enum; Postgres enum additions are safe, but keep it
  in logical order for readability.
- No data backfill needed — existing datasets simply have `ownerId = null`.
