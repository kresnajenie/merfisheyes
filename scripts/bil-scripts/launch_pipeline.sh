#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_pipeline.sh
#
# Reads samples.csv (sample_name,input_path) and for each sample submits:
#   1. map_my_cell      (no dependency — expects combined_output already exists)
#   2. process_spatial   (after map_my_cell finishes)
#   3. s3_sync_sample    (after process_spatial finishes, only with --sync)
#
# Usage:
#   ./launch_pipeline.sh                              # uses samples.csv, default species=mouse
#   ./launch_pipeline.sh my_samples.csv               # custom sample list
#   ./launch_pipeline.sh my_samples.csv human          # custom sample list + species
#   ./launch_pipeline.sh my_samples.csv mouse --sync   # enable S3 sync
#
# MapMyCells Reference Data:
#   The map_my_cell step requires reference taxonomy data.
#   Currently stored at: /bil/data/meyes/mapmycells-reference
#   To set up, download the Allen Brain Cell Atlas taxonomy files:
#     - Mouse: https://knowledge.brain-map.org/data/LVDBJAW34Y7YOLTLWKGM/summary
#     - Human: https://knowledge.brain-map.org/data/Y4E2MJPILJNA6BMIP5W/summary
#   Place the downloaded files in your reference directory and update the
#   --reference_dir path in map_my_cell.sbatch if using a different location.
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.csv}"
SPECIES="${2:-mouse}"
SYNC=false
for arg in "$@"; do
    if [ "$arg" = "--sync" ]; then
        SYNC=true
    fi
done
MEYES_BASE="/bil/data/meyes"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  MERFISH Eyes Pipeline Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "Species:     $SPECIES"
echo "S3 sync:     $(if $SYNC; then echo 'enabled (--sync)'; else echo 'DISABLED (pass --sync to enable)'; fi)"
echo "Output base: $MEYES_BASE"
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
    combined_output="${output_base}/combined_output"
    mmc_output="${output_base}/mmc_output"
    mmc_csv="${mmc_output}/mapping_output.csv"
    meyes_output="${output_base}/meyes_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:  ${input_path}"
    echo "  Output: ${output_base}"

    # # Step 1: combine_slices (skipped — combined_output already exists)
    # combine_job=$(sbatch --parsable \
    #     --job-name="combine_${sample_name}" \
    #     "${SCRIPT_DIR}/combine_slices.sbatch" \
    #     "$input_path" \
    #     "$output_base")
    # echo "  [1/4] combine_slices  -> Job ${combine_job}"

    # Step 1: map_my_cell (no dependency — combined_output already exists)
    mmc_job=$(sbatch --parsable \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$combined_output" \
        "$mmc_output" \
        "$SPECIES")
    echo "  [1/3] map_my_cell     -> Job ${mmc_job}"

    # Step 2: process_spatial (waits for map_my_cell)
    process_job=$(sbatch --parsable \
        --dependency=afterok:${mmc_job} \
        --job-name="process_${sample_name}" \
        "${SCRIPT_DIR}/process_spatial.sbatch" \
        "$combined_output" \
        "$meyes_output" \
        "$mmc_csv")
    echo "  [2/3] process_spatial -> Job ${process_job} (after ${mmc_job})"

    # Step 3: s3 sync (waits for process_spatial)
    if $SYNC; then
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${process_job} \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [3/3] s3_sync         -> Job ${sync_job} (after ${process_job})"
    else
        echo "  [3/3] s3_sync         -> SKIPPED (pass --sync to enable)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
