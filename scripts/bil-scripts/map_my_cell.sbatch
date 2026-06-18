#!/bin/bash
#SBATCH --job-name=map_my_cell
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/map_my_cell_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"    # e.g. /bil/data/meyes/ace-dip-use/combined_output
OUTPUT_BASE="$2"   # e.g. /bil/data/meyes/ace-dip-use/mmc_output
SPECIES="${3:-mouse}"  # "mouse" (default) or "human"
REFERENCE_DIR="/bil/data/meyes/mapmycells-reference"

# ── Prevent thread oversubscription ──────────────────────────
# Each worker process should use 1 thread to avoid contention
# (50 processes × 4 MKL threads = 200 threads fighting for 64 CPUs)
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

# ── Logging ──────────────────────────────────────────────────
echo "========================================"
echo "  MapMyCells — SLURM Job $SLURM_JOB_ID"
echo "========================================"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Memory:     $SLURM_MEM_PER_NODE MB"
echo "TMPDIR:     ${TMPDIR:-/tmp}"
echo "Input:      $INPUT_PATH"
echo "Output:     $OUTPUT_BASE"
echo "Species:    $SPECIES"
echo "Reference:  $REFERENCE_DIR"
echo "Started:    $(date)"
echo "========================================"
echo ""

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

echo "Python:     $(python --version)"
echo "NumPy:      $(python -c 'import numpy; print(numpy.__version__)')"
echo ""

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

START_TIME=$SECONDS

python map_my_cell.py \
    "$INPUT_PATH" \
    "$OUTPUT_BASE" \
    --species "$SPECIES" \
    --n_processors 50 \
    --reference_dir "$REFERENCE_DIR"

EXIT_CODE=$?
ELAPSED=$(( SECONDS - START_TIME ))
MINUTES=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "========================================"
echo "  Job Complete"
echo "========================================"
echo "Exit code:  $EXIT_CODE"
echo "Runtime:    ${MINUTES}m ${SECS}s"
echo "Finished:   $(date)"
echo "========================================"
