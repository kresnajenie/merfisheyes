#!/bin/bash
#SBATCH --job-name=collect_stats
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/collect_stats_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
# Usage: sbatch collect_stats.sbatch [args passed to collect_dataset_stats.py]
#   e.g. sbatch collect_stats.sbatch --datasets ace-irk-sag --species human \
#            --out /bil/data/meyes/_dataset_stats/dataset_stats.csv

# ── Logging ──────────────────────────────────────────────────
echo "========================================"
echo "  collect_dataset_stats — SLURM Job $SLURM_JOB_ID"
echo "========================================"
echo "Node:    $(hostname)"
echo "Memory:  $SLURM_MEM_PER_NODE MB"
echo "Args:    $@"
echo "Started: $(date)"
echo "========================================"
echo ""

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes

START_TIME=$SECONDS

python scripts/bil-scripts/analysis-scripts/collect_dataset_stats.py "$@"
EXIT_CODE=$?

ELAPSED=$(( SECONDS - START_TIME ))
echo ""
echo "========================================"
echo "  Job Complete"
echo "Exit code: $EXIT_CODE"
echo "Runtime:   $(( ELAPSED / 60 ))m $(( ELAPSED % 60 ))s"
echo "Finished:  $(date)"
echo "========================================"

exit $EXIT_CODE
