# BIL HPC Environment Setup

## Python Environment

```bash
# Load modules
module load gcc/11.2.0
module load python/3.10.2    # or whatever version is available

# Create venv
python3 -m venv ~/merfisheyes_env

# Activate
source ~/merfisheyes_env/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install dependencies
pip install numpy pandas scipy matplotlib

# h5py needs prebuilt binary (no HDF5 module on BIL)
pip install h5py --only-binary=:all:

# anndata (for H5AD support, pulls in h5py)
pip install anndata

# Verify
python -c "import numpy, pandas, scipy, matplotlib, anndata; print('All imports OK')"
```

## Cluster Resources

**compute partition** (default):
- 7 nodes, 80 CPUs each, ~2.8 TB RAM per node
- Max walltime: 2 days
- Jobs can share nodes

**applications partition**:
- 1 node, 80 CPUs, ~2.9 TB RAM
- Max walltime: 8 hours
- Exclusive node access

No caps on max jobs or submissions.

## Pipeline sbatch settings

### Single Cell Pipeline

| Job | CPUs | Memory | Time |
|-----|------|--------|------|
| combine_slices | 2 | 64G | 2 days |
| process_spatial | 32 | 256G | 2 days |
| s3_sync_sample | 1 | 4G | 2 days |

### Single Molecule Pipeline

| Job | CPUs | Memory | Time |
|-----|------|--------|------|
| process_single_molecule | 64 | 512G | 2 days |
| s3_sync_sm | 1 | 4G | 2 days |

## Log Directory

Logs go to `/bil/users/ijenie/meyes_process_logs/`. Create it if needed:

```bash
mkdir -p /bil/users/ijenie/meyes_process_logs
```

## Python Scripts

**Active scripts** (used by the pipeline):

| Script | Purpose |
|--------|---------|
| `combine_slices_v3.py` | Combines multi-slice MERSCOPE samples into one dataset |
| `process_spatial_data.py` | Converts H5AD/Xenium/MERSCOPE to chunked binary for viewer |
| `process_single_molecule.py` | Converts detected_transcripts.csv to per-gene binary files |
| `map_my_cell.py` | Runs MapMyCells cell type annotation |
| `update_mapping_prefix.py` | Re-prefixes S3 URLs in an existing mapping.json |

**Legacy scripts** (superseded, not used):

| Script | Superseded by |
|--------|---------------|
| `combine_slices.py` | `combine_slices_v3.py` |
| `combine_slices_v2.py` | `combine_slices_v3.py` |
| `process_h5ad.py` | `process_spatial_data.py` |

## Running the Pipelines

### Single Cell Pipeline

Chains: `combine_slices` → `process_spatial` → `s3_sync_sample`

```bash
cd ~/merfisheyes/scripts/bil-scripts

# Edit samples.csv with your samples (comma-separated: sample_name,input_path)
# Then launch
./launch_pipeline.sh

# Monitor
squeue -u $USER
```

### Single Molecule Pipeline

Chains: `process_single_molecule` → `s3_sync_sm` + copy mapping.json to meyes_output

```bash
cd ~/merfisheyes/scripts/bil-scripts

# Uses the same samples.csv
./launch_sm_pipeline.sh

# Monitor
squeue -u $USER
```

The SM pipeline runs independently from the single cell pipeline — both can be launched at the same time.

**Standalone sbatch usage:**
```bash
# Process a single sample's detected_transcripts.csv files
sbatch process_single_molecule.sbatch /bil/data/input/ /bil/data/output/ https://bucket.s3.../prefix

# Sync SM output to S3
sbatch s3_sync_sm.sbatch ace-dip-use

# Re-prefix an existing mapping.json (no sbatch needed, runs instantly)
python ~/merfisheyes/scripts/update_mapping_prefix.py /bil/data/meyes/sample/sm_output/mapping.json \
  --new-prefix https://new-bucket.s3.us-west-2.amazonaws.com/new_path
```

### Parallelism (Single Molecule)

`process_single_molecule.sbatch` supports two levels of parallelism:

- `--sample-workers N` — process N samples concurrently (default: 4)
- `--workers N` — write gene files in parallel within each sample (default: 8)

Default: 4 × 8 = 32 cores. Adjust via 4th and 5th sbatch args:
```bash
# 8 samples at a time, 8 write workers each = 64 cores
sbatch process_single_molecule.sbatch INPUT OUTPUT S3_PREFIX 8 8
```

## Data Layout

```
/bil/data/meyes/{sample_name}/
  combined_output/          # from combine_slices
    cell_metadata.csv
    cell_by_gene.csv
    check_spatial.png
    check_expression.png
  meyes_output/             # from process_spatial (synced to S3)
    manifest.json
    coords/
    expr/
    obs/
    palettes/
    mapping.json            # copied from sm_output (links _sample_id → SM S3 URLs)
  sm_output/                # from process_single_molecule (synced to S3)
    mapping.json            # maps _sample_id → S3 URLs for right-click SM viewer
    {sample_id_1}/
      manifest.json.gz
      genes/
        GENE1.bin.gz
        GENE2.bin.gz
        ...
    {sample_id_2}/
      ...
```

## S3 Destinations

- Single cell: `s3://merfisheyes-bil/bil-psc-data2/{sample_name}/`
  - Syncs with `--size-only --exclude "*.csv" --exclude "mmc_output/*"`
- Single molecule: `s3://merfisheyes-bil/bil-psc-data2/{sample_name}/sm_output/`
  - Syncs with `--size-only`
