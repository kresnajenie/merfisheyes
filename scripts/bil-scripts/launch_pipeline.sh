#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_pipeline.sh
#
# Reads samples.csv (one sample_name per line) and for each sample submits:
#   1. process_spatial  (with MMC + mask already available)
#   2. s3_sync_sample   (after process_spatial finishes)
#
# Assumes combine_slices, map_my_cell, and mask are already done.
#
# Usage:
#   ./launch_pipeline.sh                  # uses samples.csv in same dir
#   ./launch_pipeline.sh my_samples.csv   # custom sample list
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.csv}"
MEYES_BASE="/bil/data/meyes"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  MERFISH Eyes Pipeline Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "Output base: $MEYES_BASE"
echo ""

count=0

while IFS= read -r sample_name; do
    # Trim whitespace
    sample_name="$(echo "$sample_name" | xargs)"

    # Skip comments and empty lines
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
