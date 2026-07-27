# Later — Optional pipeline stages (MapMyCells, QC)

← [Back to design doc](README.md)

## Goal

Extend the single-stage worker into the multi-stage pipeline described in
[README §5.8](README.md), adding **MapMyCells cell-type annotation** and later
**QC** — without touching the core v1 flow. These are **not** part of the initial
build; this doc keeps the design ready so they slot in cleanly.

## Depends on

[Phase 3](phase-3-trigger-callback-progress.md) (working single-stage pipeline + callback + progress).

## Design principle

Stages are declared in `Dataset.processingParams.stages` and executed as **chained
AWS Batch jobs** with `dependsOn` (the cloud equivalent of the BIL SLURM `afterok`
chain). Each stage threads its output to the next via the processor's existing
flags, and posts its own progress step. When no optional stage is enabled, it stays
a single chunk job.

```
annotate ──(mapping_output.csv)──▶ chunk --mmc-csv
qc       ──(mask CSV)───────────▶ chunk --mask
```

## Stage: `annotate` (MapMyCells)

- [ ] **Opt-in per upload**: checkbox + species select (mouse/human) in the
      column-confirm step; off by default. Write into `processingParams.stages.annotate`.
- [ ] Constraints to enforce in the UI: **single-cell only**, mouse/human taxonomy,
      needs **raw counts**; **Xenium's native layout needs an adapter** (map_my_cell
      expects h5ad or `cell_metadata`+`cell_by_gene`).
- [ ] Reference taxonomy files (mouse/human `.h5` + markers `.json` + `gene.csv`)
      hosted on **S3, mounted via EFS** to the annotate job; set
      `MERFISHEYES_REFERENCE_DIR` (the script already reads it).
- [ ] Chain: annotate job → `mapping_output.csv` → chunk job `--mmc-csv`.
- [ ] Sized to the **default `computeTier`** (we do **not** replicate BIL's 512 GB —
      that was spare-HPC headroom). Use `--flatten` for a fast/light default if needed.
- [ ] Add `@aws-sdk/client-batch` `dependsOn` wiring in `/complete`.

## Stage: `qc`

- [ ] Compute per-cell QC metrics (n_genes, total_counts, pct_mito, doublet score) →
      emit as **numerical obs columns** that flow into the chunked output.
- [ ] Emit a **QC report artifact** (PNG/HTML) to `datasets/{id}/qc/` for the dashboard.
- [ ] Optional filter mask via the existing `--mask` path (the BIL artifact-mask
      mechanism from `combine_slices` is the template).

## Verification (when built)

- [ ] With annotate on: viewer shows MapMyCells cluster columns (`class_name`, etc.).
- [ ] With annotate off: behaves exactly like the v1 single-stage flow.
- [ ] Progress UI shows each stage as a distinct step.
- [ ] QC columns appear as numerical clusters; QC report is viewable.

## Notes / risks

- MapMyCells can be slow; `--flatten` (bootstrap=1) is the fast default for interactive
  uploads. Full hierarchical mapping is much heavier.
- Reference-data hosting + EFS mount is a one-time infra setup, separate from the v1 bucket.
- Keep annotate/qc strictly additive — the single-stage chunk path must remain the
  default and untouched.
