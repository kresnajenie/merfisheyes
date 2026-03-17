#!/bin/bash
#SBATCH --job-name=process_spatial
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/process_spatial_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"    # e.g. /bil/data/meyes/ace-dud-vex/combined_output
OUTPUT_PATH="$2"   # e.g. /bil/data/meyes/ace-dud-vex/meyes_output

echo "=== process_spatial ==="
echo "Input:  $INPUT_PATH"
echo "Output: $OUTPUT_PATH"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes/.venv/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts
python process_spatial_data.py \
    "$INPUT_PATH" \
    "$OUTPUT_PATH" \
    --chunk-size 1 \
    --workers 16

echo "Job finished at $(date)"
