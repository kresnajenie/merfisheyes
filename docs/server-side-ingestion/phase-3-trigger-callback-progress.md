# Phase 3 — Trigger + callback + live progress

← [Back to design doc](README.md)

## Goal

Close the loop: `/complete` submits an AWS Batch job, the Fargate worker processes
and reports back via a secured callback, the browser sees **live progress** via
Supabase Realtime, and the user gets an email with the viewer link. This is the phase
with the **real AWS provisioning**.

## Depends on

[Phase 1](phase-1-raw-upload.md) (raw in S3, `QUEUED`) + [Phase 2](phase-2-worker-container.md) (working image).

## Tasks

### App

- [ ] `/api/ingest/[id]/complete` → after `QUEUED`, call Batch `SubmitJob`
      (`@aws-sdk/client-batch`), sizing vCPU/RAM from the uploader's `computeTier`
      (README §5.6); store `batchJobId`.
- [ ] `POST /api/ingest/[id]/callback` — shared-secret/HMAC auth. Accepts
      `status` (`PROCESSING`/`COMPLETE`/`FAILED`), `progress` (stage + %), stats,
      manifest key. Writes `Dataset.status` / `processingProgress`, and on `COMPLETE`
      triggers the SES email (reuse `app/api/send-email*`).
- [ ] `GET /api/ingest/[id]/status` — polling fallback for non-Realtime clients.
- [ ] Client: subscribe to the dataset row via **Supabase Realtime**; render the
      GitHub-Actions-style staged progress (README §5.5).

### Worker (extend Phase 2 entrypoint)

- [ ] `POST callback status=PROCESSING` at start.
- [ ] Map processor stdout stages → periodic `POST callback progress {stage, %}`.
- [ ] `POST callback status=COMPLETE {stats, fingerprint, manifestKey}` on success;
      `status=FAILED {error}` on failure.

### AWS provisioning (**you — console or IaC**)

- [ ] **ECR** repo; CI builds & pushes the Phase 2 image.
- [ ] **AWS Batch on Fargate**: compute environment + job queue + job definition
      (default `computeTier` vCPU/RAM; ephemeral storage sized for large data).
- [ ] **IAM**: app creds get `batch:SubmitJob` + `batch:DescribeJobs`; the job task
      role gets read/write on the S3 bucket + `logs:*`. (SES stays in the app.)
- [ ] **Networking**: VPC/subnets/SG for Fargate with outbound HTTPS (S3 + callback);
      NAT or S3 VPC endpoint.
- [ ] **S3 lifecycle rule** on `raw/` (expire after N days; abort orphaned multipart).
- [ ] **Supabase**: enable Realtime on the `datasets` table (or a progress table).
- [ ] **Secrets**: `CALLBACK_SECRET` shared app↔worker; AWS creds available to Vercel.

## Files touched

- `app/api/ingest/[id]/complete/route.ts`, `.../callback/route.ts`, `.../status/route.ts`
- `worker/entrypoint.py` (callback POSTs)
- Client progress component + Supabase Realtime subscription
- `package.json` (`@aws-sdk/client-batch`)

## Verification

- [ ] Drop raw → upload → `SubmitJob` fires → Batch runs the container.
- [ ] Live progress advances in the UI through the stages (Realtime).
- [ ] On success: `datasets/{id}/` populated, `raw/{id}/` deleted, status `COMPLETE`,
      email arrives, viewer link works.
- [ ] On induced failure: status `FAILED`, `errorMessage` set, UI shows the error.
- [ ] Callback rejects unauthenticated/incorrect-secret requests.

## Notes / risks

- Vercel functions are short-lived — `SubmitJob` returns immediately (fire-and-forget);
  all long work is in Batch, all status via callback. Never block a Vercel route on processing.
- Callback is the single writer of DB status + email → keeps DB/SES creds out of the worker.
- Watch the callback URL: it must be reachable from Fargate (public HTTPS to the app).
