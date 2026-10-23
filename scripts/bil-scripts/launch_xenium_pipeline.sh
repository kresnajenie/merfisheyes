#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_xenium_pipeline.sh
#
# Reads xenium-samples.csv (sample_name,input_dir) and for each sample submits:
#   1. process_spatial_data    (xenium folder → meyes_output)
#   2. s3 sync                 (meyes_output, after step 1)
#
# Expected input directory structure (auto-detected by process_spatial_data.py):
#   {input_dir}/cells.csv.gz  (or cells.csv)
#   {input_dir}/cell_feature_matrix/
#
# Usage:
#   ./launch_xenium_pipeline.sh                          # uses xenium-samples.csv in same dir
#   ./launch_xenium_pipeline.sh my_samples.csv           # custom sample list
#   ./launch_xenium_pipeline.sh my_samples.csv --sync    # enable S3 sync
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/xenium-samples.csv}"
SYNC=false
for arg in "$@"; do
    if [ "$arg" = "--sync" ]; then
        SYNC=true
    fi
done
MEYES_BASE="/bil/data/meyes"
S3_BUCKET="merfisheyes-bil"
S3_PREFIX_BASE="bil-psc-data2"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  Xenium Pipeline Launcher"
echo "============================================"
echo "Sample file:  $SAMPLE_FILE"
echo "S3 sync:      $(if $SYNC; then echo 'enabled (--sync)'; else echo 'DISABLED (pass --sync to enable)'; fi)"
echo "Output base:  $MEYES_BASE"
echo ""

count=0
errors=0

while IFS=',' read -r sample_name input_dir; do
    # Trim leading/trailing whitespace
    sample_name="$(echo "$sample_name" | xargs)"
    input_dir="$(echo "$input_dir" | xargs)"

    # Skip comments, empty lines, and header
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    # Validate input directory
    if [ ! -d "$input_dir" ]; then
        echo "ERROR: Input directory not found: $input_dir"
        errors=$((errors + 1))
        continue
    fi

    count=$((count + 1))

    # Output paths
    output_base="${MEYES_BASE}/${sample_name}"
    meyes_output="${output_base}/meyes_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input dir:  ${input_dir}"
    echo "  Output:     ${meyes_output}"

    # Step 1: process_spatial (no MMC, no mask for Xenium)
    process_job=$(sbatch --parsable \
        --job-name="process_${sample_name}" \
        "${SCRIPT_DIR}/process_spatial.sbatch" \
        "$input_dir" \
        "$meyes_output")
    echo "  [1/2] process_spatial -> Job ${process_job}"

    # Step 2: s3 sync (waits for process_spatial)
    if $SYNC; then
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${process_job} \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [2/2] s3_sync         -> Job ${sync_job} (after ${process_job})"
    else
        echo "  [2/2] s3_sync         -> SKIPPED (pass --sync to enable)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
if [ "$errors" -gt 0 ]; then
    echo "  Skipped ${errors} sample(s) due to missing directories"
fi
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
