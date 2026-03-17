#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_pipeline.sh
#
# Reads samples.tsv and for each sample submits:
#   1. combine_slices  (no dependency)
#   2. process_spatial  (after combine_slices finishes)
#   3. s3_sync_sample   (after process_spatial finishes)
#
# Usage:
#   ./launch_pipeline.sh                  # uses samples.tsv in same dir
#   ./launch_pipeline.sh my_samples.tsv   # custom sample list
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.tsv}"
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

while IFS=$'\t' read -r input_path sample_name; do
    # Skip comments and empty lines
    [[ "$input_path" =~ ^#.*$ ]] && continue
    [[ -z "$input_path" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    meyes_output="${output_base}/meyes_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:  ${input_path}"
    echo "  Output: ${output_base}"

    # Step 1: combine_slices
    combine_job=$(sbatch --parsable \
        --job-name="combine_${sample_name}" \
        "${SCRIPT_DIR}/combine_slices.sbatch" \
        "$input_path" \
        "$output_base")
    echo "  [1/3] combine_slices  -> Job ${combine_job}"

    # Step 2: process_spatial (waits for combine_slices)
    process_job=$(sbatch --parsable \
        --dependency=afterok:${combine_job} \
        --job-name="process_${sample_name}" \
        "${SCRIPT_DIR}/process_spatial.sbatch" \
        "$combined_output" \
        "$meyes_output")
    echo "  [2/3] process_spatial -> Job ${process_job} (after ${combine_job})"

    # Step 3: s3 sync (waits for process_spatial)
    sync_job=$(sbatch --parsable \
        --dependency=afterok:${process_job} \
        --job-name="sync_${sample_name}" \
        "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
        "$sample_name")
    echo "  [3/3] s3_sync         -> Job ${sync_job} (after ${process_job})"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
