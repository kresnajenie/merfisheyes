#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sm_pipeline.sh
#
# Reads samples.csv (sample_name,input_path) and for each sample submits:
#   1. process_single_molecule  (no dependency)
#   2. s3_sync_sm               (after process_single_molecule finishes)
#
# Usage:
#   ./launch_sm_pipeline.sh                          # uses samples.csv in same dir
#   ./launch_sm_pipeline.sh my_samples.csv           # custom sample list
#   ./launch_sm_pipeline.sh my_samples.csv --no-sync # skip S3 sync
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.csv}"
NO_SYNC=false
for arg in "$@"; do
    if [ "$arg" = "--no-sync" ]; then
        NO_SYNC=true
    fi
done
MEYES_BASE="/bil/data/meyes"
S3_HTTPS_BASE="https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  Single Molecule Pipeline Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "S3 sync:     $(if $NO_SYNC; then echo 'DISABLED (--no-sync)'; else echo 'enabled'; fi)"
echo "Output base: $MEYES_BASE"
echo "S3 base:     $S3_HTTPS_BASE"
echo ""

count=0

while IFS=',' read -r sample_name input_path; do
    # Trim leading/trailing whitespace
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "$input_path" | xargs)"

    # Skip comments, empty lines, and header
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    sm_output="${output_base}/sm_output"
    meyes_output="${output_base}/meyes_output"
    s3_prefix="${S3_HTTPS_BASE}/${sample_name}/sm_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:      ${input_path}"
    echo "  SM output:  ${sm_output}"
    echo "  S3 prefix:  ${s3_prefix}"

    # Step 1: process_single_molecule
    process_job=$(sbatch --parsable \
        --job-name="sm_${sample_name}" \
        "${SCRIPT_DIR}/process_single_molecule.sbatch" \
        "$input_path" \
        "$sm_output" \
        "$s3_prefix")
    echo "  [1/3] process_sm -> Job ${process_job}"

    # Step 2: copy mapping.json to meyes_output (after process completes)
    copy_job=$(sbatch --parsable \
        --dependency=afterok:${process_job} \
        --job-name="cpmap_${sample_name}" \
        --wrap="mkdir -p '${meyes_output}' && cp '${sm_output}/mapping.json' '${meyes_output}/mapping.json' && echo 'Copied mapping.json to ${meyes_output}'" \
        --output="/bil/users/ijenie/meyes_process_logs/cpmap_${sample_name}_%j.log" \
        --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 --partition=compute)
    echo "  [2/3] copy mapping.json -> Job ${copy_job} (after ${process_job})"

    # Step 3: s3 sync (waits for copy)
    if ! $NO_SYNC; then
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${copy_job} \
            --job-name="sync_sm_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sm.sbatch" \
            "$sample_name")
        echo "  [3/3] s3_sync_sm -> Job ${sync_job} (after ${copy_job})"
    else
        echo "  [3/3] s3_sync_sm -> SKIPPED (--no-sync)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
