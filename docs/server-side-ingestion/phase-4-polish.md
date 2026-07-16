# Phase 4 — Polish

← [Back to design doc](README.md)

## Goal

Harden the happy path from Phase 3 and lay the groundwork for the user-facing
surfaces (dashboard, logs, tiers). Nothing here is required for a first working
end-to-end flow, but it's what makes it production-comfortable.

## Depends on

[Phase 3](phase-3-trigger-callback-progress.md) (full loop working).

## Tasks

- [ ] **Failure surfacing + retries**: Batch job retry policy; on terminal failure,
      `status=FAILED` + human-readable `errorMessage` shown in the UI; optional
      failure email.
- [ ] **CloudWatch "view log" expander**: surface the worker's stdout logs in the
      dashboard as an expandable panel (GitHub-Actions "view raw logs" feel).
- [ ] **User dashboard groundwork**: page listing `Dataset where ownerId = me`
      with status/progress; entry points to `/viewer/{id}` / `/sm-viewer/{id}`.
      (Login is currently admin-oriented — extend to normal logged-in users.)
- [ ] **`computeTier` admin controls**: view/set a user's tier; map tier → vCPU/RAM
      (and optionally per-tier Batch job queues with priority).
- [ ] **Abandoned-upload cleanup**: reconcile `UPLOADING` datasets past
      `UploadSession.expiresAt` (cron or lazy sweep) with the S3 lifecycle rule.
- [ ] **Dedup UX**: when the worker's canonical fingerprint matches an existing
      dataset, decide behaviour (link to existing vs. keep both) and surface it.

## Files touched

- New dashboard page(s) under `app/**`
- Admin tier controls (extend existing admin area)
- Cleanup route/cron
- Callback handler (failure email, fingerprint-collision handling)

## Verification

- [ ] A failed job shows a clear error + optional email; retry works.
- [ ] Dashboard lists the current user's datasets with live status.
- [ ] Changing a user's `computeTier` changes the submitted job's resources.
- [ ] Expired `UPLOADING` datasets + their partial S3 objects get cleaned up.

## Notes / risks

- Keep the dashboard read-only first; editing/deleting datasets is a separate scope.
- Tier controls are admin-only until you decide on self-serve limits.
