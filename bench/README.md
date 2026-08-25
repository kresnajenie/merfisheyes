# Benchmark databank + methodology

Drop real-world datasets here, run the discovery script, and the benchmark
harness can measure both upload paths against them. Nothing in
`bench/datasets/` is committed — only the folder structure is.

## Where things go

```
bench/datasets/
├── h5ad/              <- one .h5ad file per dataset
├── xenium/            <- one folder per export (cells.csv + cell_feature_matrix/)
├── merscope/          <- one folder per export (cell_metadata.csv + cell_by_gene.csv)
└── single-molecule/   <- one .parquet or .csv per dataset (transcripts / detected_transcripts)
```

One dataset = one top-level entry in a format folder. Name it whatever you
like; the name becomes the benchmark id (`h5ad__retroB30`).

**Big files: symlink instead of copying.** Everything here follows symlinks, so
a 30 GB export never needs duplicating:

```bash
ln -s /Volumes/data/retroB30.h5ad bench/datasets/h5ad/retroB30.h5ad
ln -s /Volumes/data/xenium_run7   bench/datasets/xenium/xenium-run7
```

A Xenium or MERSCOPE export that also contains a transcripts file can be
benchmarked on both paths — put the export in `xenium/` and symlink the
transcripts file into `single-molecule/`.

### The shared databank on this server

The real-world databank already lives at
`/data/merfisheyes-test-data/test-data-sizes/` with exactly this layout — point
the tools at it with `--root` instead of copying anything:

```bash
python3 scripts/bench/discover.py --root /data/merfisheyes-test-data/test-data-sizes
```

### Molecule tables inside single-cell exports

A Xenium export ships `transcripts.parquet` and a MERSCOPE export ships
`detected_transcripts.csv` — the same download is both a single-cell **and** a
single-molecule dataset. Collect them into `single-molecule/` as symlinks:

```bash
python3 scripts/bench/link-transcripts.py --root <databank> [--dry-run]
```

Idempotent, copies nothing, and reports what it skipped (broken links, and
`.csv.gz` — the processor accepts only `.parquet` and `.csv`, so gunzip those
first: `gunzip -k transcripts.csv.gz`).

## Inventory what you dropped

```bash
python3 scripts/bench/discover.py [--root <databank>] [--out <file.json>]
```

Writes `bench/datasets.json` and prints a table. It reads **metadata only** (no
dataset is loaded), so it stays fast on a full drop folder, and it uses the
pipeline's own `detect_input_format`, so a dataset in the wrong folder is
reported instead of being silently benchmarked as the wrong thing.

It records what's needed to interpret a result later:

- **file size** — the headline axis for the claim
- **shape** — cells × genes, or molecule count: what actually drives memory
- **notes** — structures known to blow the peak up, e.g.

  ```
  h5ad__retroB30      2.49 GB   3,009,652 cells   992 genes
    · dense f32 matrix ≈ 11.1 GB
    · obsm['X_raw'] is full-width (3,009,652×992 ≈ 11.1 GB dense f32)
  ```

  That second line is why a 2.5 GB file needed 13.3 GB of RAM — without it the
  benchmark number looks inexplicable.

### Folder size ≠ uploaded size (important for the claim)

A Xenium export is mostly morphology images, zarr and boundary files that the
pipeline never opens. Measured on the real databank:

| export | on disk | actually read |
|---|---:|---:|
| Whole_mouse_pup_65gb | 64.2 GB | **1.51 GB** |
| Human_Renal_Carcinoma_50gb | 67.5 GB | **0.28 GB** |
| Human_Lymph_Node | 42.0 GB | **0.81 GB** |
| Human_Melanoma_10gb | 13.6 GB | **0.05 GB** |

So "we handle a 65 GB Xenium export" and "we upload 65 GB" are very different
statements — the browser filters to the needed files before uploading, and the
`used` column records exactly that. Quote whichever you mean, deliberately.
MERSCOPE exports are the opposite: `cell_by_gene.csv` + `cell_metadata.csv` are
nearly the whole folder, so on-disk ≈ uploaded.

## What gets timed (and what doesn't)

The headline number is **processing only**. Byte transfer is network-bound and
says more about the tester's link than about the product, so it is measured and
reported separately, never folded into the headline.

| Phase | Browser path | Server path | In headline? |
|---|---|---|---|
| Read + parse input | h5wasm / hyparquet / PapaParse | `process_*.py` STEP 1 | **yes** |
| Coordinates, obs, palettes, DE, chunking | `GeneChunkProcessor` | STEP 2–5 | **yes** |
| Raw bytes → S3 | — | client multipart upload | no — reported as `client_upload_ms` + MB/s |
| Job submit → container RUNNING | — | Fargate provisioning | no — reported as `queue_ms` |
| S3 → container scratch | — | worker download | no — reported as `fetch_ms` |
| Chunked output → S3 | browser upload of blobs | worker upload | no — reported as `publish_ms` |

So `process_ms` means the same thing on both sides: the CPU work of turning the
input into the chunked layout the viewer reads. That is the number the claim
uses, and it is the number a faster network cannot flatter.

Every result also records **peak memory** (browser: renderer RSS + JS heap;
server: the worker's own `peak RSS` line) and an **outcome**: `ok`,
`oom` (the browser tab died / the container was OOM-killed), or `timeout`.
For the browser path the failure *is* the measurement — it's where the ceiling
is.

## Machine fingerprint

Results are namespaced per machine via `PERF_ENV` (as the existing perf suite
does), and each run records CPU model, core count, total RAM, OS, and browser
version, so "on a Mac" and "on a high-end PC" are comparable rather than
anecdotal. Results land in `bench/results/<env>/`.

## Status

Built: the drop folder and `scripts/bench/discover.py`.
Next: the local/browser runner (ceiling detection), the server runner (reusing
`scripts/perf/ingest-upload.mjs` + the CloudWatch stage parser), and the
aggregator that turns all machines' results into the comparison table.
