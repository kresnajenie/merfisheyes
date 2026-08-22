# Handoff — worker memory optimization (`feat/worker-memory-opt`)

Temporary handoff note so this work can be continued from another machine
(e.g. the Ubuntu server). Delete before merging to `develop`.

Start a Claude Code session in this repo and say:
> "Read docs/HANDOFF-worker-memory-opt.md and continue from 'Next steps'."

## Where we are (2026-08-23)

Branch `feat/worker-memory-opt` (cut from `develop` after PR #106), commit
`93c998a`, pushed to origin. No PR yet.

### What changed (outputs are byte-identical to the previous scripts)

- `scripts/process_spatial_data.py`
  - h5ad is read element-wise with h5py (`X`, `obs`, `obsm`, `var` only) —
    `anndata.read_h5ad` also materialised `raw.X`, every `layers` matrix,
    `obsp` graphs and `uns`. `X` is read block-wise straight into float32
    (dense and CSR/CSC).
  - Step 4 streams genes chunk by chunk (extract → write → drop) instead of
    holding every gene's sparse copy (`all_gene_sparse`, as large as the
    matrix) until the end. DE stats accumulate during the same pass
    (`DeAccumulator`); step 4.5 only writes them out.
  - Dense matrices are served from column-major blocks (≤512 MB, ≤1/8 of the
    matrix) — on the 3M-cell retro run step 4 was 51 of 71 minutes because
    each gene was one strided pass over the whole matrix.
  - uint32 cell indices, obs arrays freed before the peak, bounded in-flight
    chunks when writing in parallel, peak RSS logged per stage.
- `scripts/process_single_molecule.py`
  - Parquet: columns straight from the Arrow table (no table + DataFrame
    double copy). CSV: `pyarrow.csv` with a dictionary-encoded gene column
    (pandas `dtype={gene: "category"}` during parsing DOUBLES the peak — do
    not use it); pandas object → category is the no-pyarrow fallback.
  - Coordinates rounded per column into one float32 N×3 array; gene index
    via stable argsort over category codes (no groupby DataFrames); gene
    files written one at a time (bounded in-flight when parallel).
- `worker/entrypoint.py`
  - `describe_oom()`: if the processor is SIGKILLed / the cgroup records an
    `oom_kill`, the failure reason starts with "Out of memory: …" so the app
    classifies it (OOM / PLATFORM) and the auto tier escalation from PR #106
    triggers. Proven in a 400 MB container (child killed, counter 0→1).

### Measured

| Case | Old | New |
|---|---|---|
| `MERFISH_Data.h5ad` (432K cells × 300 genes, file has raw + 2 layers), **same image, Linux cgroup peak** | **13.0 GB** | **2.13 GB** |
| same, native macOS RSS (understates the old code — memory compression) | 3.34 GB | 1.83 GB |
| SM 10M molecules, parquet (macOS RSS) | 1.65 GB | 1.13 GB |
| SM 10M molecules, CSV (macOS RSS) | 1.60 GB | 1.67 GB (faster) |

Equivalence: 26 h5ad fixtures (`tests/fixtures/h5ad/valid` + `edge_cases`),
SM parquet/CSV with and without `--cell-id-col`, serial and `--workers 2`,
plus the 5.5 GB real h5ad — all decompressed output bytes identical, natively
and inside the `linux/amd64` image (only `manifest.json` "name" differs when
the input filename differs).

### AWS state

- ECR `533267148861.dkr.ecr.us-west-2.amazonaws.com/merfisheyes-worker:memopt`
  pushed 2026-08-23, digest `sha256:5282fdc1a4f5b3d85cff42701cf0874c4f5620436518633b384eda7ad07d4b4c`.
- Batch job definition `merfisheyes-ingest-memopt` rev 1 = clone of
  `merfisheyes-ingest` rev 3 (Fargate, 4 vCPU / 16 GB default, same roles,
  100 GiB ephemeral, retry 2) pointing at that image.
- `:latest` and the default `merfisheyes-ingest` definition are UNCHANGED
  (prod and dev.merfisheyes.com still run the old image).
- The app submits to the scratch definition when run with
  `AWS_BATCH_JOB_DEFINITION=merfisheyes-ingest-memopt`.

### Background you need

- retroB30 (`retroB30_all_BRBB_Axon_Virus_AnnotBRBBTransfer.h5ad`, 2.68 GB,
  3,009,652 cells × 992 genes, dense X) OOM'd on 16 GB (at load) and 32 GB
  (entering step 4) with the old scripts, and completed on 64 GB in 71 min.
  Its raw S3 copy was deleted after that success — a fresh upload is needed.
- Worker callbacks go to `NEXT_PUBLIC_APP_URL` (dev.merfisheyes.com), so the
  deployed app records progress/complete/failed. dev has PR #106 (admin
  retry + auto OOM escalation), so a 16 GB OOM would be auto-retried at 32 GB
  on the OLD image via the default job definition — visible in the job list.
- dev.merfisheyes.com (Vercel Preview) is missing `AWS_SES_ACCESS_KEY_ID` /
  `AWS_SES_SECRET_ACCESS_KEY`, so its emails fail; add them in Vercel →
  Preview and redeploy.
- AWS CLI profile for everything below: `merfisheyes-admin` (us-west-2).

## Next steps

1. **Retro test on 16 GB with the new image.** Run the app with
   `AWS_BATCH_JOB_DEFINITION=merfisheyes-ingest-memopt npm run dev`, sign in
   as the owner (standard tier = 16 GB), upload retroB30 via "Upload &
   process on server". Watch:
   ```bash
   aws batch list-jobs --job-queue merfisheyes-ingest --job-status RUNNING --profile merfisheyes-admin --region us-west-2
   aws batch describe-jobs --jobs <jobId> --profile merfisheyes-admin --region us-west-2 --query 'jobs[0].{status:status,reason:statusReason,container:container.reason,log:container.logStreamName}'
   aws logs get-log-events --log-group-name /aws/batch/job --log-stream-name <logStream> --no-start-from-head --limit 40 --profile merfisheyes-admin --region us-west-2
   ```
   The log prints `peak RSS` after load, after chunk writing, and at the end.
   Expected: ≤ ~13 GB peak, step 4 well under the old 51 min. Success =
   dataset COMPLETE with `computeTier = standard`, `processingAttempts = 1`.
2. **Promote** once it passes:
   ```bash
   aws ecr get-login-password --region us-west-2 --profile merfisheyes-admin | docker login --username AWS --password-stdin 533267148861.dkr.ecr.us-west-2.amazonaws.com
   docker pull 533267148861.dkr.ecr.us-west-2.amazonaws.com/merfisheyes-worker:memopt
   docker tag 533267148861.dkr.ecr.us-west-2.amazonaws.com/merfisheyes-worker:memopt 533267148861.dkr.ecr.us-west-2.amazonaws.com/merfisheyes-worker:latest
   docker push 533267148861.dkr.ecr.us-west-2.amazonaws.com/merfisheyes-worker:latest
   ```
   (Fargate pulls `:latest` fresh per job; no job-definition change needed.)
   Optionally deregister `merfisheyes-ingest-memopt`.
3. **PR** `feat/worker-memory-opt` → `develop` (delete this file first), then
   the `develop` → `main` release per docs/RELEASING.md.
4. If 16 GB still OOMs: the log's last `peak RSS` line says where; the
   remaining big consumers are the float32 matrix itself (~12 GB here) and
   obs/obsm — next lever would be reading X in row blocks from disk per
   chunk (slower) or a per-dataset tier recommendation from n_cells × n_genes.
