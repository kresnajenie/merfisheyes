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

| Job | CPUs | Memory | Time |
|-----|------|--------|------|
| combine_slices | 2 | 64G | 2 days |
| process_spatial | 32 | 256G | 2 days |
| s3_sync | 1 | 4G | 4 hours |

## Log Directory

Logs go to `/bil/users/ijenie/meyes_process_logs/`. Create it if needed:

```bash
mkdir -p /bil/users/ijenie/meyes_process_logs
```

## Running the Pipeline

```bash
# From the repo
cd ~/merfisheyes/scripts/bil-scripts

# Edit samples.tsv with your samples (tab-separated: input_path  sample_name)
# Then launch
./launch_pipeline.sh

# Monitor
squeue -u $USER
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
```

## S3 Destination

`s3://merfisheyes-bil/bil-psc-data2/{sample_name}/`

Syncs with `--size-only --exclude "*.csv"` (skips large raw CSVs).
