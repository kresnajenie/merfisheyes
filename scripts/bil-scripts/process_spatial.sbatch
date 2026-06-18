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

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch process_spatial.sbatch INPUT_PATH OUTPUT_PATH [MMC_CSV] [MASK_CSV]"
    echo ""
    echo "Arguments:"
    echo "  INPUT_PATH   Path to input data (H5AD file or Xenium/MERSCOPE folder)"
    echo "  OUTPUT_PATH  Output directory for chunked data"
    echo "  MMC_CSV      (optional) MapMyCells output CSV (e.g. mapping_output.csv)"
    echo "  MASK_CSV     (optional) Artifact mask CSV with boolean column"
    echo ""
    echo "Examples:"
    echo "  sbatch process_spatial.sbatch /bil/data/meyes/ace-dip-use/combined_output /bil/data/meyes/ace-dip-use/meyes_output"
    echo "  sbatch process_spatial.sbatch /input /output /path/to/mapping_output.csv"
    echo "  sbatch process_spatial.sbatch /input /output /path/to/mapping_output.csv /path/to/mask.csv"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"    # e.g. /bil/data/meyes/ace-dud-vex/combined_output
OUTPUT_PATH="$2"   # e.g. /bil/data/meyes/ace-dud-vex/meyes_output
MMC_CSV="$3"       # (optional) e.g. /bil/data/meyes/ace-dud-vex/mmc_output/mapping_output.csv
MASK_CSV="$4"      # (optional) e.g. /bil/data/meyes/ace-dud-vex/artifact_mask.csv

echo "=== process_spatial ==="
echo "Input:  $INPUT_PATH"
echo "Output: $OUTPUT_PATH"
echo "MMC:    ${MMC_CSV:-<none>}"
echo "Mask:   ${MASK_CSV:-<none>}"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Build command ───────────────────────────────────────────
cd ~/merfisheyes/scripts
CMD=(python process_spatial_data.py "$INPUT_PATH" "$OUTPUT_PATH" --chunk-size 1 --workers 16)

if [ -n "$MMC_CSV" ]; then
    CMD+=(--mmc-csv "$MMC_CSV")
fi

if [ -n "$MASK_CSV" ]; then
    CMD+=(--mask "$MASK_CSV")
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

echo "Job finished at $(date)"
