#!/bin/bash
#SBATCH --job-name=process_h5ad_sc
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/process_h5ad_sc_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch process_h5ad_sc.sbatch H5AD_PATH OUTPUT_PATH MMC_CSV"
    echo ""
    echo "Runs process_spatial_data.py on an h5ad file with MapMyCells annotations."
    echo ""
    echo "Arguments:"
    echo "  H5AD_PATH    Path to cell_by_gene.h5ad"
    echo "  OUTPUT_PATH  Output directory for chunked data (meyes_output)"
    echo "  MMC_CSV      MapMyCells output CSV (mapping_output.csv)"
    echo ""
    echo "Examples:"
    echo "  sbatch process_h5ad_sc.sbatch /bil/data/meyes/ace-low-cop/cell_by_gene.h5ad /bil/data/meyes/ace-low-cop/meyes_output /bil/data/meyes/ace-low-cop/mmc_output/mapping_output.csv"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
H5AD_PATH="$1"
OUTPUT_PATH="$2"
MMC_CSV="$3"

echo "========================================"
echo "  Process H5AD Single Cell — SLURM Job $SLURM_JOB_ID"
echo "========================================"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Memory:     $SLURM_MEM_PER_NODE MB"
echo "H5AD:       $H5AD_PATH"
echo "Output:     $OUTPUT_PATH"
echo "MMC CSV:    $MMC_CSV"
echo "Started:    $(date)"
echo "========================================"
echo ""

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

echo "Python:     $(python --version)"
echo ""

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

START_TIME=$SECONDS

python process_spatial_data.py \
    "$H5AD_PATH" \
    "$OUTPUT_PATH" \
    --chunk-size 1 \
    --workers 16 \
    --mmc-csv "$MMC_CSV"

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
