# BIL SLURM Pipeline Documentation

This directory contains SLURM batch scripts and pipeline launchers for processing spatial transcriptomics data on the BIL HPC cluster. The processed output is used by the [MERFISH Eyes](https://merfisheyes.com) web viewer.

## Table of Contents

- [Quick Start](#quick-start)
- [Pipelines](#pipelines)
  - [Pipeline A: MERSCOPE Multi-Slice](#pipeline-a-merscope-multi-slice)
  - [Pipeline B: Single Molecule](#pipeline-b-single-molecule)
  - [Pipeline C: H5AD (Single Cell + Single Molecule)](#pipeline-c-h5ad-single-cell--single-molecule)
- [Input Formats](#input-formats) (see also [INPUT_REQUIREMENTS.md](INPUT_REQUIREMENTS.md) for detailed column specs)
- [Output Structure](#output-structure)
- [Environment Setup](#environment-setup)
- [S3 Sync Configuration](#s3-sync-configuration)
- [Verifying Output](#verifying-output)
- [Resource Requirements](#resource-requirements)
- [Troubleshooting](#troubleshooting)

---

## Quick Start

```bash
# 1. Set up environment (see Environment Setup section below)
# 2. Create a samples.csv file
# 3. Run the appropriate pipeline

# For MERSCOPE multi-slice data (end-to-end: combine → annotate → process):
./launch_pipeline.sh samples.csv                               # from CSV, default species=mouse
./launch_pipeline.sh samples.csv /custom/output/base human     # custom output base + species
./launch_pipeline.sh /path/to/input /path/to/output            # single sample (absolute paths)
# CSV mode output:    {output_base}/{sample_name}/ (default output_base: /bil/data/meyes)
# Single mode output: exactly /path/to/output/

# For single molecule data:
./launch_sm_pipeline.sh samples.csv                               # from CSV
./launch_sm_pipeline.sh /path/to/input /path/to/output            # single sample

# For H5AD data (single cell + single molecule):
./launch_h5ad_pipeline.sh h5ad-samples.csv

# Monitor jobs:
squeue -u $USER
```

S3 sync is **disabled by default**. Pass `--sync <s3_prefix>` to enable S3 upload after processing.

---

## End-to-End Examples

### Case 1: MERSCOPE Multi-Slice (Single Cell + Single Molecule)

**Scenario:** You have one parent directory with multiple slices of MERSCOPE data. Each slice contains `cell_metadata.csv`, `cell_by_gene.csv`, and `detected_transcripts.csv`. You want both single cell and single molecule output.

**Input structure:**

```
/bil/data/raw/my-experiment/
  slice_1/
    .../cell_metadata.csv
    .../cell_by_gene.csv
    .../detected_transcripts.csv
  slice_2/
    .../cell_metadata.csv
    .../cell_by_gene.csv
    .../detected_transcripts.csv
  ...
```

**Run both pipelines (they run in parallel):**

```bash
# Single cell pipeline: combine → annotate → process
./launch_pipeline.sh /bil/data/raw/my-experiment /bil/data/meyes/my-experiment

# Single molecule pipeline: discover transcripts → per-gene files
./launch_sm_pipeline.sh /bil/data/raw/my-experiment /bil/data/meyes/my-experiment
```

Both pipelines use the same input and output paths. They run independently in parallel:

- **SC pipeline** combines slices into one dataset, annotates cell types, and creates chunked binary output in `meyes_output/`
- **SM pipeline** discovers `detected_transcripts.csv` in each slice, creates per-gene binary files in `sm_output/`
- The `mapping.json` is copied to `meyes_output/` so the SC viewer can right-click to open SM data

**Output structure:**

```
/bil/data/meyes/my-experiment/
  combined_output/      # from combine_slices (SC)
  mmc_output/           # from map_my_cell (SC)
  meyes_output/         # from process_spatial (SC)
    manifest.json
    coords/ expr/ obs/ palettes/
    mapping.json        # copied from sm_output (SM)
  sm_output/            # from process_single_molecule (SM)
    mapping.json
    slice_1/
      manifest.json.gz
      genes/*.bin.gz
    slice_2/
      ...
```

**With S3 sync:**

```bash
./launch_pipeline.sh /bil/data/raw/my-experiment /bil/data/meyes/my-experiment --sync s3://my-bucket/prefix
./launch_sm_pipeline.sh /bil/data/raw/my-experiment /bil/data/meyes/my-experiment --sync s3://my-bucket/prefix
```

**To view in MERFISH Eyes:**

- Single cell: `https://merfisheyes.com/viewer/from-s3?url=https://my-bucket.s3.us-west-2.amazonaws.com/prefix/my-experiment/meyes_output`
- Single molecule: `https://merfisheyes.com/sm-viewer/from-s3?url=https://my-bucket.s3.us-west-2.amazonaws.com/prefix/my-experiment/sm_output`
- Right-clicking a cell in the SC viewer opens the corresponding SM data in a split panel

---

### Case 2: Single Slice MERSCOPE (Single Cell + Single Molecule)

**Scenario:** You have a single MERSCOPE slice (not multiple slices to combine). The directory contains `cell_metadata.csv`, `cell_by_gene.csv`, and `detected_transcripts.csv`.

**Input structure:**

```
/bil/data/raw/my-single-slice/
  cell_metadata.csv
  cell_by_gene.csv
  detected_transcripts.csv
```

**Run both pipelines:**

```bash
# Single cell pipeline (combine_slices still works on a single slice — passes through)
./launch_pipeline.sh /bil/data/raw/my-single-slice /bil/data/meyes/my-single-slice

# Single molecule pipeline
./launch_sm_pipeline.sh /bil/data/raw/my-single-slice /bil/data/meyes/my-single-slice
```

**Key difference from Case 1:** Since there's only one slice, the SM pipeline generates a `mapping.json` with `linkColumn: "__all__"` instead of `_sample_id`. This means every cell in the SC viewer links to the same SM dataset (no per-slice splitting). The `--mapping-url` flag is used internally for single-file mode.

**Output structure:**

```
/bil/data/meyes/my-single-slice/
  combined_output/      # from combine_slices (passes through single slice)
  mmc_output/           # from map_my_cell
  meyes_output/         # from process_spatial
    mapping.json        # linkColumn: "__all__" → one SM dataset for all cells
  sm_output/            # from process_single_molecule
    mapping.json        # linkColumn: "__all__"
    manifest.json.gz
    genes/*.bin.gz
```

**To view:** Same as Case 1 — right-clicking any cell opens the SM viewer.

---

### Case 3: H5AD (Single Cell + Single Molecule)

**Scenario:** You have pre-converted H5AD files with a spot table. Typically from a single sample that has already been processed into H5AD format.

**Input structure:**

```
/bil/data/raw/my-h5ad-sample/
  cell_by_gene.h5ad              # Single cell expression matrix
  segmented_spot_table.csv       # Molecule coordinates
```

**Run the H5AD pipeline (handles both SC and SM in one command):**

```bash
./launch_h5ad_pipeline.sh /bil/data/raw/my-h5ad-sample /bil/data/meyes/my-h5ad-sample
```

Or with a CSV for multiple samples:

```bash
./launch_h5ad_pipeline.sh h5ad-samples.csv
```

**What happens:**

1. `map_my_cell` annotates cell types from the H5AD file
2. `process_h5ad_sm` converts the spot table CSV to per-gene binary files (runs in parallel with step 1)
3. `process_h5ad_sc` converts the H5AD + cell type annotations to chunked binary (after step 1)
4. `mapping.json` is copied to `meyes_output/` (after steps 2+3)
5. S3 sync (only with `--sync`)

**Note:** The H5AD SM processing uses custom column mappings (`--x-col x --y-col y --z-col z --gene-col gene_names --cell-id-col cell_ids`) since the spot table format differs from standard MERSCOPE `detected_transcripts.csv`.

**Output structure:**

```
/bil/data/meyes/my-h5ad-sample/
  mmc_output/           # from map_my_cell (H5AD mode)
  meyes_output/         # from process_h5ad_sc
    mapping.json        # linkColumn: "__all__"
  sm_output/            # from process_h5ad_sm
    mapping.json
    manifest.json.gz
    genes/*.bin.gz
```

**With S3 sync:**

```bash
./launch_h5ad_pipeline.sh h5ad-samples.csv --sync s3://my-bucket/prefix
```

---

## Pipelines

### Pipeline A: MERSCOPE Multi-Slice

**Purpose:** End-to-end pipeline for MERSCOPE multi-slice data. Combines raw slices, annotates cell types, converts to chunked binary format for the web viewer, and optionally uploads to S3.

**Input:** A parent directory containing sample sub-folders with raw MERSCOPE output.

**How input discovery works (BFS):** The `combine_slices` step uses breadth-first search (BFS) from each top-level child directory to find directories containing both `cell_metadata.csv` and `cell_by_gene.csv`. File matching is fuzzy — it matches keywords like "metadata" and "cell"/"gene" in filenames, so variations like `cellpose_metadata.csv` or `cell_by_gene_exons.csv` are also detected. Each top-level child directory becomes a sample slice that gets combined.

```
input_path/
  slice_1/                         # top-level child → BFS searches here
    .../cell_metadata.csv          # found via fuzzy match
    .../cell_by_gene.csv           # found via fuzzy match
  slice_2/
    .../cellpose_metadata.csv      # also matched (fuzzy)
    .../cell_by_gene_exons.csv     # also matched (fuzzy)
  ...
```

**Command:**

```bash
# CSV mode
./launch_pipeline.sh samples.csv [output_base] [species]

# Single sample mode (absolute paths)
./launch_pipeline.sh /path/to/input /path/to/output [species]

# With S3 sync
./launch_pipeline.sh samples.csv --sync s3://bucket/prefix
./launch_pipeline.sh /path/to/input /path/to/output --sync s3://bucket/prefix --region us-east-1

# Examples
./launch_pipeline.sh                                                          # uses samples.csv
./launch_pipeline.sh my_samples.csv                                           # custom CSV
./launch_pipeline.sh my_samples.csv /custom/output/base human                 # custom output + species
./launch_pipeline.sh /bil/data/raw/sample1 /bil/data/meyes/sample1            # single sample
./launch_pipeline.sh samples.csv --sync s3://merfisheyes-bil/bil-psc-data2    # + S3 sync
```

**Arguments:**

- `samples.csv` — CSV file with `sample_name,input_path` columns (default: `samples.csv` in this directory)
- `/path/to/input` + `/path/to/output` — single sample mode with absolute paths. Sample name is derived from the basename of the output path
- `output_base` — (CSV mode only) parent directory for output (default: `/bil/data/meyes`). Output goes to `{output_base}/{sample_name}/`
- `species` — `mouse` (default) or `human`
- `--sync <prefix>` — enable S3 upload. Accepts `s3://` or `https://` format
- `--region <region>` — AWS region (default: `us-west-2`). Only needed with `s3://` format

**MapMyCells Reference Data:**

The `map_my_cell` step requires Allen Brain Cell Atlas taxonomy reference files. Currently stored at `/bil/data/meyes/mapmycells-reference`.

The reference directory must contain these files:

```
mapmycells-reference/
  mouse/
    precomputed_stats_ABC_revision_230821.h5    # Precomputed statistics
    mouse_markers_230821.json                    # Marker genes
    gene.csv                                     # Gene mapping
  human/
    precomputed_stats.siletti.training.h5        # Precomputed statistics
    query_markers.n10.20240221800.json            # Marker genes
    gene.csv                                     # Gene mapping
```

To set up on a new system:

1. Download the taxonomy files from Allen Brain Map:
   - Mouse whole brain taxonomy: https://knowledge.brain-map.org/data/LVDBJAW34Y7YOLTLWKGM/summary
   - Human whole brain taxonomy: https://knowledge.brain-map.org/data/Y4E2MJPILJNA6BMIP5W/summary
2. Download the precomputed stats, marker genes, and gene mapping files
3. Place them in `mouse/` and `human/` subdirectories as shown above
4. Update `--reference_dir` in `map_my_cell.sbatch` if using a different path than `/bil/data/meyes/mapmycells-reference`
5. You can also set the `MERFISHEYES_REFERENCE_DIR` environment variable instead

**Job Chain:**

```
Step 1: combine_slices   → Combine multi-slice data into one dataset (2 CPUs, 64G, ~30min-2h)
Step 2: map_my_cell      → Cell type annotation (64 CPUs, 512G, ~1-4h)
Step 3: process_spatial  → Chunked binary output (32 CPUs, 256G, ~1-4h)
Step 4: s3_sync          → Upload to S3 (1 CPU, 4G, varies) [only with --sync <prefix>]
```

Each step depends on the previous. Steps 1-3 run for each sample; step 4 requires `--sync <prefix>`.

**Output:** `{output_base}/meyes_output/` containing:

```
meyes_output/
  manifest.json          # Dataset metadata
  coords/spatial.bin.gz  # Spatial coordinates
  expr/                  # Gene expression chunks
  obs/                   # Cell metadata columns
  palettes/              # Color palettes
```

---

### Pipeline B: Single Molecule

**Purpose:** Converts detected transcript coordinates into per-gene binary files for the single molecule viewer.

**Input:** A parent directory containing sample sub-folders. The same input directory used for the single cell pipeline.

**How input discovery works (BFS):** The `process_single_molecule` step uses breadth-first search (BFS) from each top-level child directory to find `detected_transcripts.csv` files. File matching is fuzzy — it matches keywords like "detected" and "transcript" in filenames. Each top-level child directory becomes a `_sample_id` in the output, matching the `_sample_id` column generated by `combine_slices`.

```
input_path/
  slice_1/                              # top-level child → _sample_id = "slice_1"
    .../detected_transcripts.csv        # found via fuzzy match
  slice_2/                              # top-level child → _sample_id = "slice_2"
    .../detected_transcripts.csv
  ...
```

The `_sample_id` values match between the single cell and single molecule pipelines because both derive them from the same top-level directory names. This is what enables right-click linking between the SC and SM viewers.

**Command:**

```bash
# CSV mode
./launch_sm_pipeline.sh samples.csv [output_base]

# Single sample mode (absolute paths)
./launch_sm_pipeline.sh /path/to/input /path/to/output

# With S3 sync
./launch_sm_pipeline.sh samples.csv --sync s3://bucket/prefix
./launch_sm_pipeline.sh /path/to/input /path/to/output --sync https://bucket.s3.region.amazonaws.com/prefix
```

**Arguments:**

- `samples.csv` — CSV file with `sample_name,input_path` columns
- `/path/to/input` + `/path/to/output` — single sample mode with absolute paths
- `output_base` — (CSV mode only) parent directory for output (default: `/bil/data/meyes`)
- `--sync <prefix>` — enable S3 upload. Accepts `s3://` or `https://` format. When enabled:
  - `mapping.json` contains full S3 URLs for right-click SM viewer linking
  - `mapping.json` is copied to `meyes_output/` for SC viewer integration
  - `sm_output/` is uploaded to S3
- `--region <region>` — AWS region (default: `us-west-2`). Only needed with `s3://` format.

**Without `--sync`:** `mapping.json` uses relative paths. You can add S3 URLs later with `update_mapping_prefix.py`.

**Job Chain:**

```
Step 1: process_single_molecule → Per-gene binary files (64 CPUs, 512G, ~1-8h)
Step 2: copy mapping.json       → Links SM data to SC viewer [only with --sync]
Step 3: s3_sync_sm              → Upload to S3 [only with --sync]
```

**Output:** `{output_path}/sm_output/` containing:

```
sm_output/
  mapping.json              # Maps _sample_id → S3 URLs (or relative paths without --sync)
  {sample_id}/
    manifest.json.gz        # Per-sample metadata
    genes/
      GENE1.bin.gz          # Gzipped Float32Array coordinates
      GENE2.bin.gz
      GENE1_uuuuuuuuuu.bin.gz  # Unassigned molecules (if present)
      ...
```

**Parallelism:** The processing script uses two levels:

- `--sample-workers 4` — process 4 samples concurrently
- `--workers 8` — write gene files in parallel within each sample
- Default: 4 × 8 = 32 cores

---

### Pipeline C: H5AD (Single Cell + Single Molecule)

**Purpose:** Processes pre-converted H5AD files alongside spot tables. Produces both single cell and single molecule output in one pipeline.

**Input:** A directory containing:

- `cell_by_gene.h5ad` — Single cell expression matrix (required)
- `segmented_spot_table.csv` — Molecule coordinates (required)

The launch script validates both files exist before submitting.

**Command:**

```bash
./launch_h5ad_pipeline.sh h5ad-samples.csv [--sync]
```

**Job Chain:**

```
Step 1: map_my_cell (H5AD)       → Cell type annotation (64 CPUs, 512G)
Step 2: process_h5ad_sm (CSV)    → SM per-gene files (64 CPUs, 512G) [parallel with step 1]
Step 3: process_h5ad_sc (H5AD)   → Chunked SC output (32 CPUs, 256G) [after step 1]
Step 4: copy mapping.json        → Links SM to SC (1 CPU, 1G) [after steps 2+3]
Step 5: s3_sync                  → Upload both outputs (8 CPUs, 4G) [after step 4]
```

Steps 1 and 2 run in parallel. Step 3 waits for step 1. Step 4 waits for both 2 and 3.

**Output:** Both `meyes_output/` and `sm_output/` (see above for structure).

---

## Input Formats

### samples.csv

Standard CSV format used by all pipelines:

```csv
sample_name,input_path
ace-dip-use,/bil/data/path/to/merfish_output
ace-dud-vex,/bil/data/path/to/other/merfish_output
```

**Rules:**

- First line can be a header (will be skipped if it starts with a letter)
- Lines starting with `#` are comments (skipped)
- Empty lines are skipped
- Whitespace around values is trimmed
- `sample_name` is used for the output directory name under the output base (default `/bil/data/meyes/`)
- `input_path` points to the raw data directory

### h5ad-samples.csv

Same format as `samples.csv`, but `input_path` must point to a directory containing:

- `cell_by_gene.h5ad`
- `segmented_spot_table.csv`

```csv
sample_name,input_path
ace-low-cub,/bil/data/50/51/505118bf1edfb3e6/702265/1296400199/merfish_output/1296400199_SISSeg_mousedevmodel_20241010
```

---

## Output Structure

After processing, data is organized under `{output_base}/{sample_name}/` (default output base: `/bil/data/meyes`):

```
{output_base}/{sample_name}/
  combined_output/           # From combine_slices (MERSCOPE pipeline only)
    cell_metadata.csv
    cell_by_gene.csv
    check_spatial.png        # Sanity check plots
    check_expression.png
    artifact_mask_p*.csv     # Artifact filter masks (p5 through p80)
  mmc_output/                # From map_my_cell
    mapping_output.csv       # Cell type annotations
  meyes_output/              # From process_spatial (→ synced to S3)
    manifest.json
    coords/
    expr/
    obs/
    palettes/
    mapping.json             # Copied from sm_output (only when --sync is used)
  sm_output/                 # From process_single_molecule (→ synced to S3)
    mapping.json
    {sample_id}/
      manifest.json.gz
      genes/
        *.bin.gz
```

---

## Environment Setup

### Option 1: Python Virtual Environment

```bash
# Load modules
module load gcc/11.2.0
module load python/3.10.2

# Create and activate venv
python3 -m venv ~/merfisheyes_env
source ~/merfisheyes_env/bin/activate

# Install dependencies
pip install --upgrade pip
pip install numpy pandas scipy matplotlib
pip install h5py --only-binary=:all:
pip install anndata

# Verify
python -c "import numpy, pandas, scipy, matplotlib, anndata; print('All imports OK')"
```

### Option 2: Singularity Container

If using Singularity, ensure the container has:

- Python 3.10+
- numpy, pandas, scipy, matplotlib
- h5py, anndata
- AWS CLI (for S3 sync)

Modify the sbatch scripts to run inside your container:

```bash
# Example: replace the python command in sbatch scripts
singularity exec /path/to/container.sif python process_spatial_data.py ...
```

You'll also need to update the `module load` and `source activate` lines in the sbatch scripts.

### Log Directory

All job logs are written to `/bil/users/ijenie/meyes_process_logs/`. Create it if needed:

```bash
mkdir -p /bil/users/ijenie/meyes_process_logs
```

**Important:** If running under a different user, update the log directory path in all sbatch scripts and the `--mail-user` email addresses.

---

## S3 Sync Configuration

S3 sync is **disabled by default**. To enable, pass `--sync <prefix>` with your S3 bucket URL:

```bash
# Using s3:// format (region defaults to us-west-2)
./launch_pipeline.sh samples.csv --sync s3://my-bucket/my-prefix
./launch_sm_pipeline.sh samples.csv --sync s3://my-bucket/my-prefix

# Using https:// format (region auto-detected from URL)
./launch_pipeline.sh samples.csv --sync https://my-bucket.s3.us-west-2.amazonaws.com/my-prefix

# Custom region (only needed with s3:// format)
./launch_sm_pipeline.sh samples.csv --sync s3://my-bucket/prefix --region us-east-1
```

The launch scripts automatically convert between `s3://` and `https://` formats:

- `s3://bucket/prefix` → used for `aws s3 sync` upload
- `https://bucket.s3.region.amazonaws.com/prefix` → used for `mapping.json` browser URLs (SM pipeline)

### S3 Upload Destinations

When `--sync s3://bucket/prefix` is provided:

- **Single cell:** uploads to `s3://bucket/prefix/{sample_name}/`
- **Single molecule:** uploads to `s3://bucket/prefix/{sample_name}/sm_output/`

### AWS Credentials

AWS credentials must be configured before using `--sync`. Run `aws configure` or set environment variables:

```bash
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
export AWS_DEFAULT_REGION=us-west-2
```

### Updating mapping.json After Processing

If you processed without `--sync` (relative paths in mapping.json), you can add S3 URLs later:

```bash
python scripts/update_mapping_prefix.py /path/to/sm_output/mapping.json \
  --new-prefix https://bucket.s3.us-west-2.amazonaws.com/prefix/sample_name/sm_output
```

---

## Verifying Output

### Single Cell Output

After `process_spatial` completes, verify:

```bash
ls {output_base}/{sample_name}/meyes_output/
# Expected: manifest.json  coords/  expr/  obs/  palettes/

# Check manifest
cat {output_base}/{sample_name}/meyes_output/manifest.json | python -m json.tool | head -20
# Should show: normalized: false, statistics with total_cells and total_genes
```

### Single Molecule Output

After `process_single_molecule` completes, verify:

```bash
ls {output_base}/{sample_name}/sm_output/
# Expected: mapping.json  {sample_id_1}/  {sample_id_2}/ ...

# Check mapping
cat {output_base}/{sample_name}/sm_output/mapping.json | python -m json.tool
# Should show: linkColumn and links with sample_id → S3 URL mappings

# Check a sample directory
ls {output_base}/{sample_name}/sm_output/{sample_id}/
# Expected: manifest.json.gz  genes/

ls {output_base}/{sample_name}/sm_output/{sample_id}/genes/ | head -5
# Expected: GENE1.bin.gz  GENE2.bin.gz ...
```

### MapMyCells Output

```bash
ls {output_base}/{sample_name}/mmc_output/
# Expected: mapping_output.csv

head -2 {output_base}/{sample_name}/mmc_output/mapping_output.csv
# Should show cell type annotation columns
```

---

## Resource Requirements

| Job                     | CPUs | Memory | Typical Duration | Max Time |
| ----------------------- | ---- | ------ | ---------------- | -------- |
| combine_slices          | 2    | 64G    | 30min–4h         | 2 days   |
| map_my_cell             | 64   | 512G   | 1–2h             | 2 days   |
| process_spatial         | 32   | 256G   | 30min            | 2 days   |
| process_single_molecule | 64   | 512G   | 30min            | 2 days   |
| s3_sync_sample          | 1    | 4G     | 10min–2h         | 2 days   |
| s3_sync_sm              | 8    | 4G     | 10min–4h         | 2 days   |

## Troubleshooting

### Job failed — where are the logs?

```bash
ls /bil/users/ijenie/meyes_process_logs/
# Logs are named: {job_type}_{sample_name}_{job_id}.log
```

### MapMyCells hangs or OOM

- Check that `MKL_NUM_THREADS=1` is set (prevents thread oversubscription)
- Default uses 50 parallel processes × 1 thread = 50 cores on 64 CPUs
- Reference data must exist at `/bil/data/meyes/mapmycells-reference`

### Process spatial runs out of memory

- Default: 256G. For very large datasets (>1M cells), increase to 512G
- Edit `process_spatial.sbatch` line: `#SBATCH --mem=512G`

### S3 sync is slow

- SM sync uses 8 parallel uploads by default
- AWS concurrent requests set to 50
- For very large datasets, increase `--cpus-per-task` in sync sbatch

### combine_slices fails to find files

- Script uses fuzzy file matching (keywords: "cell_metadata", "cell_by_gene")
- Check that CSV files exist in sample subdirectories
- Run with single sample for debugging: `sbatch combine_slices.sbatch /path/to/input /path/to/output`

### "pyarrow not found" error in process_single_molecule

- pyarrow is only needed for `.parquet` input files
- For `.csv` input (like `detected_transcripts.csv`), pyarrow is not required
- To install: `pip install pyarrow` (may fail on BIL due to old cmake — use CSV instead)

---

## Script Reference

### Launch Scripts (run these)

| Script                    | Purpose                                                     | Input                 |
| ------------------------- | ----------------------------------------------------------- | --------------------- |
| `launch_pipeline.sh`      | MERSCOPE SC pipeline (map_my_cell → process_spatial → sync) | samples.csv           |
| `launch_sm_pipeline.sh`   | Single molecule pipeline (process → copy mapping → sync)    | samples.csv           |
| `launch_h5ad_pipeline.sh` | H5AD SC + SM pipeline (parallel processing → sync)          | h5ad-samples.csv      |
| `launch_combine_mmc.sh`   | Combine multi-slice + artifact mask + map_my_cell           | samples.csv           |
| `launch_process_sync.sh`  | Process + sync only (for already-combined data)             | samples.csv           |
| `launch_sync.sh`          | Re-sync to S3 only                                          | CSV with sample names |

### SBATCH Jobs (called by launch scripts)

| Script                           | Purpose                           | Resources    |
| -------------------------------- | --------------------------------- | ------------ |
| `combine_slices.sbatch`          | Combine multi-slice MERSCOPE data | 2 CPU, 64G   |
| `map_my_cell.sbatch`             | Cell type annotation              | 64 CPU, 512G |
| `process_spatial.sbatch`         | SC chunked binary output          | 32 CPU, 256G |
| `process_single_molecule.sbatch` | SM per-gene binary output         | 64 CPU, 512G |
| `process_h5ad_sc.sbatch`         | H5AD SC processing                | 32 CPU, 256G |
| `process_h5ad_sm.sbatch`         | H5AD SM processing                | 64 CPU, 512G |
| `s3_sync_sample.sbatch`          | Sync SC output to S3              | 1 CPU, 4G    |
| `s3_sync_sm.sbatch`              | Sync SM output to S3 (parallel)   | 8 CPU, 4G    |
| `sync_mapping.sbatch`            | Sync all mapping.json files       | 1 CPU, 1G    |

### Python Scripts (called by sbatch jobs)

Located in `scripts/` (parent directory):

| Script                       | Purpose                                         |
| ---------------------------- | ----------------------------------------------- |
| `combine_slices_v3.py`       | Combines multi-slice MERSCOPE samples           |
| `process_spatial_data.py`    | Converts H5AD/Xenium/MERSCOPE to chunked binary |
| `process_single_molecule.py` | Converts transcripts to per-gene binary files   |
| `map_my_cell.py`             | Cell type annotation via Allen Brain MapMyCells |
| `update_mapping_prefix.py`   | Re-prefixes S3 URLs in mapping.json             |

---

## Processing Steps In Detail

### combine_slices_v3.py

Combines multiple MERSCOPE tissue slices into a single unified dataset. Also works for single slices (passes through).

**What it does:**

1. **Discovers sample directories** — BFS from each top-level child of the input directory to find `cell_metadata.csv` + `cell_by_gene.csv` (fuzzy matched)

2. **Standardizes cell IDs** — For each sample:

   - `cell_metadata.csv`: searches for ID column in order: `EntityID` → `id` → `cell_id`. If none found, falls back to the first column (renamed to `cell_id`)
   - `cell_by_gene.csv`: searches for `cell` column. If not found, uses first column
   - IDs are converted to strings for consistent matching
   - Cell IDs must match between metadata and expression files — rows are aligned by inner join on common IDs

3. **Separates spatial coordinates** — Shifts each slice so they don't overlap, with configurable padding (default: 1000 units) between slices. Arranges slices in a grid layout.

4. **Adds tracking columns:**

   - `_source_file` — absolute path to the original `cell_metadata.csv` for each cell
   - `_sample_id` — top-level directory name for each slice (used for SM linking)

5. **Generates artifact mask** — Computes per-cell total gene expression (row sums) and generates artifact mask CSVs at multiple percentile thresholds (5th, 10th, 15th, ..., 80th in 5% increments). Cells below the threshold are flagged as artifacts (low-expression border cells).

**Output (`combined_output/`):**

```
combined_output/
  cell_metadata.csv              # Combined metadata with shifted coordinates
  cell_by_gene.csv               # Combined expression matrix
  check_spatial.png              # Spatial layout sanity check
  check_expression.png           # Gene expression sanity check
  check_artifact_thresholds.png  # Comparison of artifact filtering at different percentiles
  artifact_mask_p5.csv           # Artifact mask at 5th percentile
  artifact_mask_p10.csv          # ... 10th
  artifact_mask_p15.csv          # ... 15th
  ...                            # (every 5th percentile)
  artifact_mask_p75.csv          # ... 75th
  artifact_mask_p80.csv          # ... 80th
```

Each `artifact_mask_p*.csv` has two columns: `cell_id` and `is_artifact` (boolean). You choose which percentile to use when running `process_spatial_data.py` with the `--mask` flag.

**Mask-only mode:** If you've already combined the data but want to regenerate the artifact masks:

```bash
sbatch combine_slices.sbatch --mask-only /path/to/combined_output /path/to/combined_output
```

---

### map_my_cell.py

Annotates cells with cell type labels using the Allen Brain Cell Atlas taxonomy via the [MapMyCells algorithm](https://github.com/AllenInstitute/cell_type_mapper/blob/main/examples/mapping_to_subset_of_abc_atlas_data.ipynb).

**What it does:**

1. **Loads input** — reads either:

   - `combined_output/` directory (CSV mode): loads `cell_metadata.csv` + `cell_by_gene.csv`
   - `.h5ad` file (H5AD mode): loads expression from `adata.X` (or `adata.obsm['X_raw']` if X appears normalized)

2. **Translates gene names** — converts gene symbols to Ensembl IDs using the reference `gene.csv` mapping. Genes without an Ensembl match are dropped.

3. **Builds query H5AD** — creates a temporary `query.h5ad` file in the format required by the MapMyCells mapper

4. **Runs hierarchical mapping** — assigns each cell to a taxonomy node using the reference precomputed statistics and marker genes. Uses parallel processing (default: 50 workers, 1 thread each).

5. **Generates output** — joins mapping results back to the original metadata and produces sanity check plots

**Output (`mmc_output/`):**

```
mmc_output/
  mapping_output.csv             # Cell type annotations (one row per cell)
  query.h5ad                     # Intermediate query file (can be deleted)
  check_celltype_*.png           # Spatial plots of top 3 most common cell types
```

The `mapping_output.csv` contains columns like `class_name`, `subclass_name`, `supertype_name`, `cluster_label`, plus bootstrapping probability scores for each level.

**Reference directory structure:**

The reference data must be at `/bil/data/meyes/mapmycells-reference` (or set via `--reference_dir` or `MERFISHEYES_REFERENCE_DIR` env var):

```
mapmycells-reference/
  mouse/
    precomputed_stats_ABC_revision_230821.h5
    mouse_markers_230821.json
    gene.csv
  human/
    precomputed_stats.siletti.training.h5
    query_markers.n10.20240221800.json
    gene.csv
```

Download from:

- Mouse: https://knowledge.brain-map.org/data/LVDBJAW34Y7YOLTLWKGM/summary
- Human: https://knowledge.brain-map.org/data/Y4E2MJPILJNA6BMIP5W/summary

**Species:** Default is `mouse`. Pass `--species human` for human data. The species determines which reference files and taxonomy configuration to use.

---

### process_spatial_data.py

Converts spatial transcriptomics data into the chunked binary format used by the MERFISH Eyes web viewer. This is the final processing step before S3 upload.

**What it does:**

1. **Auto-detects input format** — H5AD, Xenium folder, or MERSCOPE folder (can be overridden with `--format`)
2. **Rounds spatial coordinates** to 2 decimal places (raw microns, no normalization)
3. **Processes expression matrix** — creates chunked sparse binary files for on-demand gene loading
4. **Processes observation columns** — creates compressed JSON files for each metadata/cluster column
5. **Generates color palettes** — for categorical columns

**Can be run without map_my_cell:**

`process_spatial_data.py` does NOT require MapMyCells output. The `--mmc-csv` flag is optional. Without it, the viewer will only show the metadata columns that are already in the input data (e.g., columns from `cell_metadata.csv` or H5AD `obs`).

```bash
# Without MapMyCells (no cell type annotations)
python process_spatial_data.py /path/to/combined_output /path/to/meyes_output --chunk-size 1 --workers 16

# With MapMyCells (adds cell type columns to viewer)
python process_spatial_data.py /path/to/combined_output /path/to/meyes_output --chunk-size 1 --workers 16 \
  --mmc-csv /path/to/mmc_output/mapping_output.csv

# With artifact mask (filters low-expression border cells)
python process_spatial_data.py /path/to/combined_output /path/to/meyes_output --chunk-size 1 --workers 16 \
  --mmc-csv /path/to/mmc_output/mapping_output.csv \
  --mask /path/to/combined_output/artifact_mask_p65.csv
```

**Optional flags:**

| Flag           | Description                                                             | Default                                      |
| -------------- | ----------------------------------------------------------------------- | -------------------------------------------- |
| `--mmc-csv`    | MapMyCells output CSV — adds cell type annotation columns to the viewer | None (skip)                                  |
| `--mask`       | Artifact mask CSV — filters cells where `is_artifact=true`              | None (keep all cells)                        |
| `--mask-col`   | Column name in mask CSV to filter on                                    | `is_artifact`                                |
| `--mask-keep`  | Value to keep (case-insensitive)                                        | `false` (keep cells where is_artifact=false) |
| `--chunk-size` | Genes per expression chunk                                              | Auto-determined                              |
| `--workers`    | Parallel workers for chunk writing                                      | 1                                            |

**Mask behavior:** The `--mask` flag expects a CSV with a boolean column (default: `is_artifact`). Cells where the column value equals `--mask-keep` (default: `false`, meaning NOT an artifact) are kept. All other cells are removed. Choose a percentile mask from the artifact masks generated by `combine_slices` (e.g., `artifact_mask_p65.csv` filters at the 65th percentile).

**MapMyCells behavior:** The `--mmc-csv` flag adds all non-empty columns from the MapMyCells output as observation columns in the viewer. The `cell_id` column is used for alignment and then dropped. If `--mask` is also used, the mask is applied to the MMC data as well.

**Output (`meyes_output/`):**

```
meyes_output/
  manifest.json          # Dataset metadata (cell count, gene count, dimensions, etc.)
  coords/
    spatial.bin.gz       # Spatial coordinates (binary, compressed)
  expr/
    index.json           # Gene → chunk mapping
    chunk_00000.bin.gz   # Expression data chunks (sparse format)
    chunk_00001.bin.gz
    ...
  obs/
    metadata.json        # Column types (categorical/numerical)
    class_name.json.gz   # Cluster values (from MMC or metadata)
    subclass_name.json.gz
    ...
  palettes/
    class_name.json      # Color palettes for categorical columns
    ...
```
