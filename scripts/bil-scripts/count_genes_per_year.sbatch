#!/bin/bash
#SBATCH --job-name=count_genes_per_year
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/count_genes_per_year_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
# Usage: sbatch count_genes_per_year.sbatch [extra args passed to the python script]
#   e.g. sbatch count_genes_per_year.sbatch \
#            --out-dir /bil/data/meyes/_genes_inventory \
#            --walk-timeout 120

# ── Logging ──────────────────────────────────────────────────
echo "========================================"
echo "  count_genes_per_year — SLURM Job $SLURM_JOB_ID"
echo "========================================"
echo "Node:    $(hostname)"
echo "CPUs:    $SLURM_CPUS_PER_TASK"
echo "Memory:  $SLURM_MEM_PER_NODE MB"
echo "Args:    $@"
echo "Started: $(date)"
echo "========================================"
echo ""

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

echo "Python:  $(python --version)"
echo ""

# ── Run ──────────────────────────────────────────────────────
# Script reads --datasets-csv relative to CWD; cd to repo root so the
# default `scripts/bil-scripts/all_datasets.csv` resolves correctly.
cd ~/merfisheyes

START_TIME=$SECONDS

python scripts/bil-scripts/count_genes_per_year.py "$@"
EXIT_CODE=$?

ELAPSED=$(( SECONDS - START_TIME ))
MINUTES=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "========================================"
echo "  Job Complete"
echo "========================================"
echo "Exit code: $EXIT_CODE"
echo "Runtime:   ${MINUTES}m ${SECS}s"
echo "Finished:  $(date)"
echo "========================================"

exit $EXIT_CODE
