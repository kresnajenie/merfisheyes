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
# Single-sample h5ad inputs (input is an .h5ad or a dir with cell_by_gene.h5ad)
# SKIP combine and feed the h5ad straight to MMC — they need no slice merging.
#
# Note: combine_slices_v3.py generates the artifact masks as part of its normal
# run, so no separate mask job is submitted here. To regenerate masks for an
# already-combined dataset, run combine_slices.sbatch with --mask-only manually.
#
# Stats accumulate in two CSVs under ${MEYES_BASE}/_dataset_stats/:
#   dataset_stats.csv  one row per dataset (counts at every taxonomy level,
#                      mean bootstrapping + mean correlation)
#   group_stats.csv    one row per (dataset, method, level, group) with
#                      per-group cell count + mean bootstrapping/correlation
# Re-running a sample replaces its rows in both.
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
GROUP_CSV="${MEYES_BASE}/_dataset_stats/group_stats.csv"

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

    # Detect single-sample h5ad input — these skip combine and feed MMC directly
    # (input is either an .h5ad file or a dir containing cell_by_gene.h5ad).
    h5ad_input=""
    if [[ "$input_path" == *.h5ad ]]; then
        h5ad_input="$input_path"
    elif [ -f "${input_path}/cell_by_gene.h5ad" ]; then
        h5ad_input="${input_path}/cell_by_gene.h5ad"
    fi

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:    ${input_path}"
    if [ -n "$h5ad_input" ]; then
        echo "  Type:     h5ad (combine skipped) → ${h5ad_input}"
    else
        echo "  Type:     csv → Combined: ${combined_output}"
    fi
    echo "  MMC hier: ${mmc_output}"
    echo "  MMC corr: ${mmc_corr_output}"
    echo "  Species:  ${species}"

    # Step 1: combine_slices (also writes artifact masks) — skipped for h5ad
    mmc_dep=()
    if [ -n "$h5ad_input" ]; then
        mmc_input="$h5ad_input"
        echo "  [1/4] combine_slices  -> SKIPPED (single-sample h5ad)"
    else
        combine_job=$(sbatch --parsable \
            --job-name="combine_${sample_name}" \
            "${SCRIPT_DIR}/combine_slices.sbatch" \
            "$input_path" \
            "$output_base")
        mmc_input="$combined_output"
        mmc_dep=(--dependency=afterok:${combine_job})
        echo "  [1/4] combine_slices  -> Job ${combine_job}"
    fi

    # Step 2: map_my_cell — hierarchical (after combine if any)
    mmc_job=$(sbatch --parsable \
        ${mmc_dep[@]+"${mmc_dep[@]}"} \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$mmc_input" \
        "$mmc_output" \
        "$species" \
        "hierarchical")
    echo "  [2/4] mmc hierarchical -> Job ${mmc_job}"

    # Step 3: map_my_cell — correlation (parallel to hierarchical)
    mmc_corr_job=$(sbatch --parsable \
        ${mmc_dep[@]+"${mmc_dep[@]}"} \
        --job-name="mmccorr_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$mmc_input" \
        "$mmc_corr_output" \
        "$species" \
        "correlation")
    echo "  [3/4] mmc correlation  -> Job ${mmc_corr_job}"

    # Step 4: collect dataset stats — needs hierarchical to succeed and
    # correlation to finish (afterany, so a flaky correlation run doesn't block
    # stats). For h5ad, pass the file so combine-derived metrics come from it.
    stats_h5ad=()
    [ -n "$h5ad_input" ] && stats_h5ad=(--h5ad "$h5ad_input")
    stats_job=$(sbatch --parsable \
        --dependency=afterok:${mmc_job},afterany:${mmc_corr_job} \
        --job-name="stats_${sample_name}" \
        "${SCRIPT_DIR}/collect_stats.sbatch" \
        --datasets "$sample_name" \
        --species "$species" \
        --out "$STATS_CSV" \
        --group-out "$GROUP_CSV" \
        ${stats_h5ad[@]+"${stats_h5ad[@]}"})
    echo "  [4/4] collect_stats    -> Job ${stats_job} (after hier+corr)"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
