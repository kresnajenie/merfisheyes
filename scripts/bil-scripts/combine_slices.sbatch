#!/bin/bash
#SBATCH --job-name=combine_slices
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/combine_slices_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
# Normal mode:    sbatch combine_slices.sbatch <input_path> <output_base>
# Mask-only mode: sbatch combine_slices.sbatch --mask-only <combined_output_dir> [extra_percentiles]
#   extra_percentiles: optional comma-separated list (e.g. "80") added to the
#   standard 5-75 sweep, for percentiles outside the default range.

MODE="$1"

if [ "$MODE" == "--mask-only" ]; then
    COMBINED_DIR="$2"
    EXTRA_PERCENTILES="${3:-}"
    echo "=== combine_slices (mask-only mode) ==="
    echo "Combined dir: $COMBINED_DIR"
    [ -n "$EXTRA_PERCENTILES" ] && echo "Extra percentiles: $EXTRA_PERCENTILES"
else
    INPUT_PATH="$1"
    OUTPUT_BASE="$2"
    echo "=== combine_slices ==="
    echo "Input:  $INPUT_PATH"
    echo "Output: $OUTPUT_BASE"
fi

echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

if [ "$MODE" == "--mask-only" ]; then
    CMD=(python combine_slices_v3.py --mask-only "$COMBINED_DIR" "$COMBINED_DIR")
    [ -n "$EXTRA_PERCENTILES" ] && CMD+=(--percentiles "$EXTRA_PERCENTILES")
    "${CMD[@]}"
else
    python combine_slices_v3.py \
        "$INPUT_PATH" \
        "$OUTPUT_BASE"
fi

echo "Job finished at $(date)"
