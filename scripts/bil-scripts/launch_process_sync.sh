#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_process_sync.sh
#
# Runs process_spatial + s3_sync only. Assumes map_my_cell is already done
# (mmc_output/mapping_output.csv exists). Auto-detects input type per
# sample, same as launch_mmc.sh:
#   - ${MEYES_BASE}/${sample_name}/combined_output   (CSV/MERSCOPE)
#       -> process_spatial.sbatch (combined_output + artifact mask)
#   - ${input_path}/cell_by_gene.h5ad                 (BIL h5ad)
#       -> process_h5ad_sm.sbatch + process_h5ad_sc.sbatch + copy mapping.json
#          (h5ad has no combine/mask step -- CSV/MERSCOPE-only)
# Both then run s3_sync_sample.sbatch.
#
# Usage:
#   ./launch_process_sync.sh ace-dip-use          # single sample (combined_output)
#   ./launch_process_sync.sh samples.csv          # from file (sample_name,input_path[,species[,percentile]])
#   ./launch_process_sync.sh                      # uses samples.csv in same dir
#
# percentile column is optional and defaults to 25 when absent (artifact_mask_p<N>.csv).
# Only applies to combined_output (CSV/MERSCOPE) samples; ignored for h5ad (no mask).
# The mask CSV for that percentile must already exist -- generate it first with
# combine_slices.sbatch --mask-only <combined_output_dir> <percentile> if it's
# outside the standard 5-75 sweep (e.g. 80).
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"
S3_BUCKET="merfisheyes-bil"
S3_PREFIX_BASE="bil-psc-data2"
S3_HTTPS_BASE="https://${S3_BUCKET}.s3.us-west-2.amazonaws.com/${S3_PREFIX_BASE}"

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
errors=0

while IFS=',' read -r sample_name input_path _species percentile; do
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "${input_path:-}" | xargs)"
    percentile="$(echo "${percentile:-}" | xargs)"
    percentile="${percentile:-25}"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    mmc_output="${output_base}/mmc_output"
    meyes_output="${output_base}/meyes_output"
    mask_csv="${combined_output}/artifact_mask_p${percentile}.csv"

    echo "── Sample ${count}: ${sample_name} ──"

    if [ -d "$combined_output" ]; then
        # ── CSV/MERSCOPE: process_spatial with mask ──────────────
        echo "  Kind:     combined_output (CSV/MERSCOPE)"
        echo "  MMC:      ${mmc_output}/mapping_output.csv"
        echo "  Mask:     ${mask_csv} (p${percentile})"
        echo "  Output:   ${meyes_output}"

        if [ ! -f "$mask_csv" ]; then
            echo "  ERROR: mask not found: $mask_csv"
            echo "         Generate it first: sbatch combine_slices.sbatch --mask-only '${combined_output}' ${percentile}"
            errors=$((errors + 1)); echo ""; continue
        fi

        process_job=$(sbatch --parsable \
            --job-name="process_${sample_name}" \
            "${SCRIPT_DIR}/process_spatial.sbatch" \
            "$combined_output" \
            "$meyes_output" \
            "${mmc_output}/mapping_output.csv" \
            "$mask_csv")
        echo "  [1/2] process_spatial -> Job ${process_job}"

        sync_job=$(sbatch --parsable \
            --dependency=afterok:${process_job} \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [2/2] s3_sync         -> Job ${sync_job} (after ${process_job})"

    elif [ -n "$input_path" ] && [ -f "${input_path}/cell_by_gene.h5ad" ]; then
        # ── BIL h5ad: no combine/mask, process sc + sm then sync ──
        h5ad_path="${input_path}/cell_by_gene.h5ad"
        csv_path="${input_path}/segmented_spot_table.csv"
        sm_output="${output_base}/sm_output"
        sm_s3_url="${S3_HTTPS_BASE}/${sample_name}/sm_output"
        mmc_csv="${mmc_output}/mapping_output.csv"

        echo "  Kind:     h5ad (no mask)"
        echo "  H5AD:     ${h5ad_path}"
        echo "  MMC:      ${mmc_csv}"
        echo "  Output:   ${meyes_output}"

        if [ ! -f "$csv_path" ]; then
            echo "  ERROR: segmented_spot_table.csv not found: $csv_path"
            errors=$((errors + 1)); echo ""; continue
        fi
        if [ ! -f "$mmc_csv" ]; then
            echo "  ERROR: mapping_output.csv not found (run launch_mmc.sh first): $mmc_csv"
            errors=$((errors + 1)); echo ""; continue
        fi

        sm_job=$(sbatch --parsable --job-name="sm_${sample_name}" \
            "${SCRIPT_DIR}/process_h5ad_sm.sbatch" "$csv_path" "$sm_output" "$sm_s3_url")
        echo "  [1/4] process_sm       -> Job ${sm_job}"

        sc_job=$(sbatch --parsable --job-name="sc_${sample_name}" \
            "${SCRIPT_DIR}/process_h5ad_sc.sbatch" "$h5ad_path" "$meyes_output" "$mmc_csv")
        echo "  [2/4] process_spatial  -> Job ${sc_job}"

        copy_job=$(sbatch --parsable \
            --dependency=afterok:${sc_job}:${sm_job} \
            --job-name="cpmap_${sample_name}" \
            --wrap="cp '${sm_output}/mapping.json' '${meyes_output}/mapping.json'" \
            --output="/bil/users/ijenie/meyes_process_logs/cpmap_h5ad_${sample_name}_%j.log" \
            --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 --partition=compute)
        echo "  [3/4] copy mapping.json -> Job ${copy_job}"

        sync_job=$(sbatch --parsable \
            --dependency=afterok:${copy_job} \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [4/4] s3_sync          -> Job ${sync_job} (after ${copy_job})"

    else
        echo "  ⚠️  SKIP: no combined_output dir and no cell_by_gene.h5ad at input_path"
        errors=$((errors + 1))
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
if [ "$errors" -gt 0 ]; then
    echo "  Skipped/errored ${errors} sample(s)"
fi
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
