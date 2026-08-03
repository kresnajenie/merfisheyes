# Phase 1 — Raw upload (app-side, no worker yet)

← [Back to design doc](README.md)

## Goal

Let a logged-in user drop raw data and upload the **raw bytes** directly to S3
(`raw/{datasetId}/`), with the full upload UX (overlay, MB/MB, ETA, close guards,
cancel) and the mandatory column-confirm step. At the end of this phase the dataset
reaches `QUEUED` — **no processing yet**. Everything here is app-side; no AWS
provisioning beyond the S3 bucket you already have.

## Depends on

[Phase 0](phase-0-schema.md) (needs `QUEUED`, `ownerId`, `processingParams`).

## Tasks

### API routes (new `app/api/ingest/*`, mirroring `app/api/datasets/*`)

- [ ] `POST /api/ingest/initiate` — auth-gate (NextAuth session); create `Dataset`
      (`status=UPLOADING`, `ownerId` from session, `ingestSource="server"`,
      `processingParams` from the request), `UploadSession`, `UploadFile[]`; return
      presigned upload URLs (single-PUT or multipart create).
- [ ] `POST /api/ingest/[id]/files/[key]/complete` — mark a raw `UploadFile`
      `COMPLETE`, increment `completedFiles` (reuse the single-molecule pattern).
- [ ] `POST /api/ingest/[id]/complete` — verify all files uploaded → set `QUEUED`.
      **(SubmitJob is Phase 3 — this phase stops at `QUEUED`.)**
- [ ] `POST /api/ingest/[id]/abort` — `AbortMultipartUpload`, delete uploaded
      `raw/{id}/` objects, mark dataset aborted/deleted.

### S3 helpers ([`lib/s3.ts`](../../lib/s3.ts))

- [ ] Multipart helpers for the **raw path only, files > 100 MB**:
      `CreateMultipartUpload` → presigned `UploadPart` URLs → `CompleteMultipartUpload`
      → `AbortMultipartUpload`. Small files keep single-PUT (unchanged).
- [ ] Part size 16–64 MB; return part URLs in batches.

### Client — column extraction & confirm (README §7)

- [ ] CSV: read first `~64 KB` via `File.slice` → parse header + sample rows.
- [ ] Parquet: read footer via `hyparquet` (tail bytes) → column names.
- [ ] Column-mapping UI: pre-fill Xenium/MERSCOPE defaults, **require explicit
      "Continue"** even when defaults match; write result into `processingParams`.

### Client — upload overlay (README §10)

- [ ] Full-screen blocking overlay + big blue bar; copy "please don't close this tab".
- [ ] Per-PUT progress via `XMLHttpRequest` (`fetch` has no upload progress);
      multipart aggregate = Σ(done parts) + in-flight part.
- [ ] Aggregate `MB/MB` + rolling-average ETA (~5 s window; "Calculating…" until warm).
- [ ] Close guard: `beforeunload` (native generic) + custom in-app-nav modal.
- [ ] Cancel button → calls `/abort`, aborts in-flight XHRs (`AbortController`),
      removes guards.
- [ ] Completion state: "✅ Uploaded — we'll email you"; user free to leave.

## Files touched

- `app/api/ingest/**` (new)
- `lib/s3.ts` (multipart helpers)
- New drop-zone component (complementary to [`components/file-upload.tsx`](../../components/file-upload.tsx))
- New upload-overlay + column-confirm components

## Verification

- [ ] Raw files land at `raw/{id}/…` in S3 (single-PUT small, multipart >100 MB).
- [ ] `Dataset` reaches `QUEUED` with the chosen `processingParams` and `ownerId`.
- [ ] Byte counter + ETA are accurate on a large multi-file upload.
- [ ] Cancel aborts the upload and leaves **no** orphaned S3 objects or multipart uploads.
- [ ] Tab-close prompts (native) and in-app-nav modal (custom) both fire.
- [ ] Column-confirm blocks submit until the user clicks Continue.

## Notes / risks

- The existing browser drag-drop (local parse) is **untouched**; this is a separate,
  complementary path.
- Add an S3 lifecycle rule to abort orphaned incomplete multipart uploads (belt-and-braces
  with `/abort`).
- Fingerprint/dedup is **not** computed here — the worker does it post-processing (README §6).
