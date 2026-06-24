# BIL Pipeline Test Run — Handoff

End-to-end smoke test of the MERFISH Eyes pipeline on BIL. Produces both
single-cell and single-molecule output for 7 samples. No S3 upload.

## What this does

For each sample in `samples-ivan.csv`, two pipelines run in parallel:

- **Single cell**: `combine_slices → map_my_cell → process_spatial`
  → chunked binary in `meyes_output/`.
- **Single molecule**: discovers `detected_transcripts.csv` per slice
  → per-gene binary files in `sm_output/`.

Output lands at **`/bil/data/meyes/test-ivan/{sample_name}/`**.

## Prerequisites you need to have working on BIL

Before running anything, make sure:

1. **Python env** is set up per `README.md` ("Environment Setup"). The
   sbatch scripts activate it via `module load gcc/11.2.0 python/3.10.2`
   and `source ~/merfisheyes_env/bin/activate` (or your equivalent path).
2. **MapMyCells reference data** for mouse exists at
   `/bil/data/meyes/mapmycells-reference/mouse/`:
   - `precomputed_stats_ABC_revision_230821.h5`
   - `mouse_markers_230821.json`
   - `gene.csv`

   Download from <https://knowledge.brain-map.org/data/LVDBJAW34Y7YOLTLWKGM/summary>.
3. **Log directory** that the sbatch scripts write to exists and is
   writable by your user. The scripts hardcode
   `/bil/users/ijenie/meyes_process_logs/` — change those paths and the
   `--mail-user=` lines to your user before running. (`README.md` →
   "Environment Setup" → "Log Directory".)

If those aren't ready, the SLURM jobs will fail at the appropriate step
with a clear error.

## CSV format

`samples-ivan.csv` is two columns:

```
sample_name,input_path
ace-dip-tug,/bil/data/63/05/63053886eb9f9771/660403
...
```

- Header line is auto-skipped.
- All samples in this CSV are mouse (passed via the command line below).

The CSV mixes two input shapes — the pipeline handles both, but the
generated `_sample_id` differs:

- **Multi-slice** (`ace-dip-*`, `ace-dug-fin`): `input_path` is a parent
  directory whose top-level children are individual slices. Each child
  dir name becomes one `_sample_id`.
- **Single-slice** (`ace-low-*`): `input_path` is a single MERSCOPE
  output folder with `cell_metadata.csv` / `cell_by_gene.csv` /
  `detected_transcripts.csv` at the top. Produces one `_sample_id`.

## How to run

From `scripts/bil-scripts/` on BIL:

```bash
# Single-cell pipeline
./launch_pipeline.sh samples-ivan.csv /bil/data/meyes/test-ivan mouse

# Single-molecule pipeline (run in a second terminal, they're parallel)
./launch_sm_pipeline.sh samples-ivan.csv /bil/data/meyes/test-ivan
```

Both commands submit jobs and return immediately. Monitor with:

```bash
squeue -u $USER
```

Per-sample logs go to your log directory (see prerequisite #3).

## Output structure

After all jobs complete, each sample has:

```
/bil/data/meyes/test-ivan/{sample_name}/
├── combined_output/      # from combine_slices (SC)
├── mmc_output/           # cell type annotations (SC)
├── meyes_output/         # chunked binary for the viewer (SC)
│   ├── manifest.json
│   ├── coords/
│   ├── expr/
│   ├── obs/
│   └── palettes/
└── sm_output/            # per-gene binary files (SM)
    ├── mapping.json
    └── {sample_id}/
        ├── manifest.json.gz
        └── genes/*.bin.gz
```

## Verifying it worked

For any sample:

```bash
ls /bil/data/meyes/test-ivan/{sample_name}/meyes_output/
# Expected: manifest.json  coords/  expr/  obs/  palettes/

ls /bil/data/meyes/test-ivan/{sample_name}/sm_output/
# Expected: mapping.json  {sample_id}/  ...
```

If both directories have the expected files, the test ran successfully.

## Viewing without S3

S3 sync is off for this test, so the viewer can't open these via
`/viewer/from-s3?url=...`. To check the output visually, you can
drag-drop the local folders into the MERFISH Eyes web viewer:

- Single cell: drop `/bil/data/meyes/test-ivan/{sample}/meyes_output/`
  onto the **Folder** card on <https://merfisheyes.com> (single-cell
  mode). It will be detected as a chunked dataset and rendered locally,
  no upload.
- Single molecule: drop
  `/bil/data/meyes/test-ivan/{sample}/sm_output/{sample_id}/` onto the
  **Chunked Folder** card in single-molecule mode.

S3 upload can be added later with `--sync s3://bucket/prefix` once the
bucket is configured (see `README.md` → "S3 Sync Configuration").

## If something fails

- Check the per-sample log files in your log directory.
- See `README.md` → "Troubleshooting" for the common failure modes
  (MapMyCells OOM, BFS not finding files, etc.).
