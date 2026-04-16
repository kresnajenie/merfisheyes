# BIL SLURM Pipeline Documentation

This directory contains SLURM batch scripts and pipeline launchers for processing spatial transcriptomics data on the BIL HPC cluster. The processed output is used by the [MERFISH Eyes](https://merfisheyes.com) web viewer.

## Table of Contents

- [Quick Start](#quick-start)
- [Pipelines](#pipelines)
  - [Pipeline A: MERSCOPE Multi-Slice](#pipeline-a-merscope-multi-slice)
  - [Pipeline B: Single Molecule](#pipeline-b-single-molecule)
  - [Pipeline C: H5AD (Single Cell + Single Molecule)](#pipeline-c-h5ad-single-cell--single-molecule)
- [Input Formats](#input-formats)
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

# For MERSCOPE data (requires combined_output/ to exist):
./launch_pipeline.sh samples.csv human --no-sync

# For single molecule data:
./launch_sm_pipeline.sh samples.csv --no-sync

# For H5AD data (single cell + single molecule):
./launch_h5ad_pipeline.sh h5ad-samples.csv --no-sync

# Monitor jobs:
squeue -u $USER
```

All launch scripts accept `--no-sync` to skip S3 upload (recommended for testing).

---

## Pipelines

### Pipeline A: MERSCOPE Multi-Slice

**Purpose:** Processes combined MERSCOPE multi-slice data into chunked binary format for the web viewer.

**Prerequisite:** The raw multi-slice data must already be combined into `combined_output/`. If not, run `launch_combine_mmc.sh` first.

**Input:** A directory containing `combined_output/` with:
- `cell_metadata.csv`
- `cell_by_gene.csv`

**Command:**
```bash
./launch_pipeline.sh samples.csv [species] [--no-sync]
```

**Arguments:**
- `samples.csv` — CSV file with `sample_name,input_path` columns (default: `samples.csv` in this directory)
- `species` — `human` (default) or `mouse`
- `--no-sync` — skip S3 upload after processing

**Job Chain:**
```
Step 1: map_my_cell     → Cell type annotation (64 CPUs, 512G, ~2h)
Step 2: process_spatial  → Chunked binary output (32 CPUs, 256G, ~1-4h)
Step 3: s3_sync          → Upload to S3 (1 CPU, 4G, varies)
```

Each step depends on the previous. Steps 1-2 run for each sample; step 3 is optional (`--no-sync`).

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

**Input:** A parent directory containing sample sub-folders, each with a `detected_transcripts.csv` (auto-detected via fuzzy matching).

**Command:**
```bash
./launch_sm_pipeline.sh samples.csv [--no-sync]
```

**Job Chain:**
```
Step 1: process_single_molecule → Per-gene binary files (64 CPUs, 512G, ~1-8h)
Step 2: copy mapping.json       → Links SM data to SC viewer (1 CPU, 1G, <1min)
Step 3: s3_sync_sm              → Upload to S3 with 8 parallel uploads (8 CPUs, 4G, varies)
```

**Output:** `{output_base}/sm_output/` containing:
```
sm_output/
  mapping.json              # Maps _sample_id → S3 URLs
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
./launch_h5ad_pipeline.sh h5ad-samples.csv [--no-sync]
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
- `sample_name` is used for the output directory name under `/bil/data/meyes/`
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

After processing, data is organized under `/bil/data/meyes/{sample_name}/`:

```
/bil/data/meyes/{sample_name}/
  combined_output/           # From combine_slices (MERSCOPE pipeline only)
    cell_metadata.csv
    cell_by_gene.csv
    check_spatial.png        # Sanity check plots
    check_expression.png
    artifact_mask_p25.csv    # Artifact filter mask
  mmc_output/                # From map_my_cell
    mapping_output.csv       # Cell type annotations
  meyes_output/              # From process_spatial (→ synced to S3)
    manifest.json
    coords/
    expr/
    obs/
    palettes/
    mapping.json             # Copied from sm_output (links to SM viewer)
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

By default, all launch scripts include S3 sync as the final step. To skip sync (for testing or when S3 is not configured), pass `--no-sync`:

```bash
./launch_pipeline.sh samples.csv human --no-sync
./launch_sm_pipeline.sh samples.csv --no-sync
./launch_h5ad_pipeline.sh h5ad-samples.csv --no-sync
```

### Configuring S3

The sync scripts use these defaults:
- **Bucket:** `merfisheyes-bil`
- **Prefix:** `bil-psc-data2/{sample_name}/`
- **Region:** `us-west-2`

To change these, edit the following files:
- `s3_sync_sample.sbatch` — lines with `SOURCE` and `DESTINATION`
- `s3_sync_sm.sbatch` — lines with `SOURCE` and `DESTINATION`
- `launch_h5ad_pipeline.sh` — `S3_BUCKET` and `S3_PREFIX_BASE` variables
- `launch_sm_pipeline.sh` — `S3_HTTPS_BASE` variable

AWS credentials must be configured. Run `aws configure` or set environment variables:
```bash
export AWS_ACCESS_KEY_ID=...
export AWS_SECRET_ACCESS_KEY=...
export AWS_DEFAULT_REGION=us-west-2
```

---

## Verifying Output

### Single Cell Output

After `process_spatial` completes, verify:
```bash
ls /bil/data/meyes/{sample_name}/meyes_output/
# Expected: manifest.json  coords/  expr/  obs/  palettes/

# Check manifest
cat /bil/data/meyes/{sample_name}/meyes_output/manifest.json | python -m json.tool | head -20
# Should show: normalized: false, statistics with total_cells and total_genes
```

### Single Molecule Output

After `process_single_molecule` completes, verify:
```bash
ls /bil/data/meyes/{sample_name}/sm_output/
# Expected: mapping.json  {sample_id_1}/  {sample_id_2}/ ...

# Check mapping
cat /bil/data/meyes/{sample_name}/sm_output/mapping.json | python -m json.tool
# Should show: linkColumn and links with sample_id → S3 URL mappings

# Check a sample directory
ls /bil/data/meyes/{sample_name}/sm_output/{sample_id}/
# Expected: manifest.json.gz  genes/

ls /bil/data/meyes/{sample_name}/sm_output/{sample_id}/genes/ | head -5
# Expected: GENE1.bin.gz  GENE2.bin.gz ...
```

### MapMyCells Output

```bash
ls /bil/data/meyes/{sample_name}/mmc_output/
# Expected: mapping_output.csv

head -2 /bil/data/meyes/{sample_name}/mmc_output/mapping_output.csv
# Should show cell type annotation columns
```

---

## Resource Requirements

| Job | CPUs | Memory | Typical Duration | Max Time |
|-----|------|--------|-----------------|----------|
| combine_slices | 2 | 64G | 30min–2h | 2 days |
| map_my_cell | 64 | 512G | 1–4h | 2 days |
| process_spatial | 32 | 256G | 1–4h | 2 days |
| process_single_molecule | 64 | 512G | 1–8h | 2 days |
| s3_sync_sample | 1 | 4G | 10min–2h | 2 days |
| s3_sync_sm | 8 | 4G | 10min–4h | 2 days |

**Cluster partitions:**
- `compute` — 7 nodes, 80 CPUs each, ~2.8 TB RAM per node, 2-day max walltime
- `applications` — 1 node, 80 CPUs, ~2.9 TB RAM, 8-hour max walltime

All jobs use the `compute` partition by default.

---

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

| Script | Purpose | Input |
|--------|---------|-------|
| `launch_pipeline.sh` | MERSCOPE SC pipeline (map_my_cell → process_spatial → sync) | samples.csv |
| `launch_sm_pipeline.sh` | Single molecule pipeline (process → copy mapping → sync) | samples.csv |
| `launch_h5ad_pipeline.sh` | H5AD SC + SM pipeline (parallel processing → sync) | h5ad-samples.csv |
| `launch_combine_mmc.sh` | Combine multi-slice + artifact mask + map_my_cell | samples.csv |
| `launch_process_sync.sh` | Process + sync only (for already-combined data) | samples.csv |
| `launch_sync.sh` | Re-sync to S3 only | CSV with sample names |

### SBATCH Jobs (called by launch scripts)

| Script | Purpose | Resources |
|--------|---------|-----------|
| `combine_slices.sbatch` | Combine multi-slice MERSCOPE data | 2 CPU, 64G |
| `map_my_cell.sbatch` | Cell type annotation | 64 CPU, 512G |
| `process_spatial.sbatch` | SC chunked binary output | 32 CPU, 256G |
| `process_single_molecule.sbatch` | SM per-gene binary output | 64 CPU, 512G |
| `process_h5ad_sc.sbatch` | H5AD SC processing | 32 CPU, 256G |
| `process_h5ad_sm.sbatch` | H5AD SM processing | 64 CPU, 512G |
| `s3_sync_sample.sbatch` | Sync SC output to S3 | 1 CPU, 4G |
| `s3_sync_sm.sbatch` | Sync SM output to S3 (parallel) | 8 CPU, 4G |
| `sync_mapping.sbatch` | Sync all mapping.json files | 1 CPU, 1G |

### Python Scripts (called by sbatch jobs)

Located in `scripts/` (parent directory):

| Script | Purpose |
|--------|---------|
| `combine_slices_v3.py` | Combines multi-slice MERSCOPE samples |
| `process_spatial_data.py` | Converts H5AD/Xenium/MERSCOPE to chunked binary |
| `process_single_molecule.py` | Converts transcripts to per-gene binary files |
| `map_my_cell.py` | Cell type annotation via Allen Brain MapMyCells |
| `update_mapping_prefix.py` | Re-prefixes S3 URLs in mapping.json |
