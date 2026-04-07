#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_process_sync.sh
#
# Runs process_spatial + s3_sync only.
# Assumes combine_slices, map_my_cell, and mask are already done.
#
# Usage:
#   ./launch_process_sync.sh ace-dip-use          # single sample
#   ./launch_process_sync.sh samples.csv          # from file
#   ./launch_process_sync.sh                      # uses samples.csv in same dir
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"

# If argument is a file, read samples from it; otherwise treat as a single sample name
if [ $# -eq 0 ]; then
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
elif [ -f "$1" ]; then
    SAMPLE_FILE="$1"
else
    SAMPLE_FILE=$(mktemp)
    echo "$1" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  Process Spatial + S3 Sync"
echo "============================================"
echo "Source: $SAMPLE_FILE"
echo ""

count=0

while IFS=',' read -r sample_name _rest; do
    sample_name="$(echo "$sample_name" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    mmc_output="${output_base}/mmc_output"
    meyes_output="${output_base}/meyes_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Combined: ${combined_output}"
    echo "  MMC:      ${mmc_output}/mapping_output.csv"
    echo "  Mask:     ${combined_output}/artifact_mask_p25.csv"
    echo "  Output:   ${meyes_output}"

    # Step 1: process_spatial (with MMC + mask)
    process_job=$(sbatch --parsable \
        --job-name="process_${sample_name}" \
        "${SCRIPT_DIR}/process_spatial.sbatch" \
        "$combined_output" \
        "$meyes_output" \
        "${mmc_output}/mapping_output.csv" \
        "${combined_output}/artifact_mask_p25.csv")
    echo "  [1/2] process_spatial -> Job ${process_job}"

    # Step 2: s3 sync (waits for process_spatial)
    sync_job=$(sbatch --parsable \
        --dependency=afterok:${process_job} \
        --job-name="sync_${sample_name}" \
        "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
        "$sample_name")
    echo "  [2/2] s3_sync         -> Job ${sync_job} (after ${process_job})"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
