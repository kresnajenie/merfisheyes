#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_combine_mmc.sh
#
# Runs combine_slices -> artifact mask -> map_my_cell.
# Use launch_process_sync.sh afterwards for process_spatial + s3_sync.
#
# Usage:
#   ./launch_combine_mmc.sh ace-dip-use /bil/data/18/aa/.../input   # single sample
#   ./launch_combine_mmc.sh samples.csv                              # from file (sample_name,input_path)
#   ./launch_combine_mmc.sh                                          # uses samples.csv in same dir
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"

# Parse arguments
if [ $# -eq 0 ]; then
    # No args — read from default samples.csv
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
    SINGLE_MODE=false
elif [ $# -eq 1 ] && [ -f "$1" ]; then
    # One arg that's a file — read from it
    SAMPLE_FILE="$1"
    SINGLE_MODE=false
elif [ $# -eq 2 ]; then
    # Two args — single sample mode (sample_name, input_path)
    SAMPLE_FILE=$(mktemp)
    echo "$1,$2" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
    SINGLE_MODE=true
else
    echo "Usage:"
    echo "  $0 <sample_name> <input_path>   # single sample"
    echo "  $0 <samples.csv>                # from CSV (sample_name,input_path)"
    echo "  $0                              # uses samples.csv in same dir"
    exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  Combine + Mask + MapMyCells"
echo "============================================"
echo "Source: $SAMPLE_FILE"
echo ""

count=0

while IFS=',' read -r sample_name input_path; do
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "$input_path" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    mmc_output="${output_base}/mmc_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:    ${input_path}"
    echo "  Combined: ${combined_output}"
    echo "  MMC:      ${mmc_output}"

    # Step 1: combine_slices
    combine_job=$(sbatch --parsable \
        --job-name="combine_${sample_name}" \
        "${SCRIPT_DIR}/combine_slices.sbatch" \
        "$input_path" \
        "$output_base")
    echo "  [1/3] combine_slices  -> Job ${combine_job}"

    # Step 2: artifact mask (after combine)
    filter_job=$(sbatch --parsable \
        --dependency=afterok:${combine_job} \
        --job-name="filter_${sample_name}" \
        "${SCRIPT_DIR}/combine_slices.sbatch" \
        --mask-only \
        "$combined_output")
    echo "  [2/3] artifact_mask   -> Job ${filter_job} (after ${combine_job})"

    # Step 3: map_my_cell (after mask)
    mmc_job=$(sbatch --parsable \
        --dependency=afterok:${filter_job} \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$combined_output" \
        "$mmc_output")
    echo "  [3/3] map_my_cell     -> Job ${mmc_job} (after ${filter_job})"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
