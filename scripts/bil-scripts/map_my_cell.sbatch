#!/bin/bash
#SBATCH --job-name=map_my_cell
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/map_my_cell_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"    # e.g. /bil/data/meyes/ace-dip-use/combined_output
OUTPUT_BASE="$2"   # e.g. /bil/data/meyes/ace-dip-use/mmc_output

echo "=== map my cells ==="
echo "Input:  $INPUT_PATH"
echo "Output: $OUTPUT_BASE"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts
python map_my_cell.py \
    "$INPUT_PATH" \
    "$OUTPUT_BASE"

echo "Job finished at $(date)"