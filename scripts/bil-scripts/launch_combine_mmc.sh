#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_combine_mmc.sh
#
# Runs combine_slices (which also writes artifact masks) -> map_my_cell
# (both hierarchical AND correlation methods) -> collect_dataset_stats.
# Use launch_process_sync.sh afterwards for process_spatial + s3_sync.
#
# Per sample: mmc_output/ holds the hierarchical mapping, mmc_output_corr/
# holds the correlation mapping. Both depend only on combine and run in parallel.
#
# Note: combine_slices_v3.py generates the artifact masks as part of its normal
# run, so no separate mask job is submitted here. To regenerate masks for an
# already-combined dataset, run combine_slices.sbatch with --mask-only manually.
#
# Stats for every sample accumulate in one CSV (see STATS_CSV below), one row
# per dataset; re-running a sample replaces its row.
#
# Usage:
#   ./launch_combine_mmc.sh ace-dip-use /bil/data/18/aa/.../input [species]  # single sample
#   ./launch_combine_mmc.sh samples.csv                                       # from CSV (sample_name,input_path,species)
#   ./launch_combine_mmc.sh                                                   # uses samples.csv in same dir
#
# species column is optional and defaults to "mouse" when absent.
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"
STATS_CSV="${MEYES_BASE}/_dataset_stats/dataset_stats.csv"

# Parse arguments
if [ $# -eq 0 ]; then
    # No args — read from default samples.csv
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
    SINGLE_MODE=false
elif [ $# -eq 1 ] && [ -f "$1" ]; then
    # One arg that's a file — read from it
    SAMPLE_FILE="$1"
    SINGLE_MODE=false
elif [ $# -eq 2 ] || [ $# -eq 3 ]; then
    # Single sample mode (sample_name, input_path, [species])
    SAMPLE_FILE=$(mktemp)
    echo "$1,$2,${3:-mouse}" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
    SINGLE_MODE=true
else
    echo "Usage:"
    echo "  $0 <sample_name> <input_path> [species]   # single sample"
    echo "  $0 <samples.csv>                          # from CSV (sample_name,input_path,species)"
    echo "  $0                                        # uses samples.csv in same dir"
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

while IFS=',' read -r sample_name input_path species; do
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "$input_path" | xargs)"
    species="$(echo "${species:-}" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue
    # Default species when the column is absent or empty
    species="${species:-mouse}"

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    mmc_output="${output_base}/mmc_output"            # hierarchical method
    mmc_corr_output="${output_base}/mmc_output_corr"  # correlation method

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:    ${input_path}"
    echo "  Combined: ${combined_output}"
    echo "  MMC hier: ${mmc_output}"
    echo "  MMC corr: ${mmc_corr_output}"
    echo "  Species:  ${species}"

    # Step 1: combine_slices (also writes artifact masks)
    combine_job=$(sbatch --parsable \
        --job-name="combine_${sample_name}" \
        "${SCRIPT_DIR}/combine_slices.sbatch" \
        "$input_path" \
        "$output_base")
    echo "  [1/4] combine_slices  -> Job ${combine_job}"

    # Step 2: map_my_cell — hierarchical (after combine)
    mmc_job=$(sbatch --parsable \
        --dependency=afterok:${combine_job} \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$combined_output" \
        "$mmc_output" \
        "$species" \
        "hierarchical")
    echo "  [2/4] mmc hierarchical -> Job ${mmc_job} (after ${combine_job})"

    # Step 3: map_my_cell — correlation (after combine, parallel to hierarchical)
    mmc_corr_job=$(sbatch --parsable \
        --dependency=afterok:${combine_job} \
        --job-name="mmccorr_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$combined_output" \
        "$mmc_corr_output" \
        "$species" \
        "correlation")
    echo "  [3/4] mmc correlation  -> Job ${mmc_corr_job} (after ${combine_job})"

    # Step 4: collect dataset stats (after hierarchical MMC, which it reads)
    stats_job=$(sbatch --parsable \
        --dependency=afterok:${mmc_job} \
        --job-name="stats_${sample_name}" \
        "${SCRIPT_DIR}/collect_stats.sbatch" \
        --datasets "$sample_name" \
        --species "$species" \
        --out "$STATS_CSV")
    echo "  [4/4] collect_stats    -> Job ${stats_job} (after ${mmc_job})"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
