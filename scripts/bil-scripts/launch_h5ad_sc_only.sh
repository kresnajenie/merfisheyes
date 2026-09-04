#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_h5ad_sc_only.sh
#
# Single-cell-only H5AD processing, for datasets that already have SM
# output from a prior run (e.g. ace-owl: SM was processed first, and
# cell_by_gene.h5ad was only confirmed present afterward). Does NOT
# re-run process_single_molecule -- it reuses the existing sm_output.
#
# Also skips MapMyCells entirely: MMC only has mouse/human reference
# taxonomies (see scripts/map_my_cell.py TAXONOMY_CONFIG), so it can't
# run for other species (e.g. marmoset) regardless. The dataset still
# processes fine without MMC cluster columns.
#
# For each sample: process_spatial_data (h5ad → meyes_output) → copy
# the existing sm_output/mapping.json into meyes_output → s3 sync
# (optional, whole sample dir; --size-only makes the already-synced
# sm_output half a no-op).
#
# Reads samples.csv (sample_name,input_dir) -- same format as
# xenium-samples.csv. Fuzzy-matches *cell_by_gene.h5ad within input_dir.
#
# Requires: ${MEYES_BASE}/{sample_name}/sm_output/mapping.json must
# already exist (from a prior SM-only run). Samples missing it are
# skipped with an error.
#
# Usage:
#   ./launch_h5ad_sc_only.sh xenium-samples.csv
#   ./launch_h5ad_sc_only.sh xenium-samples.csv --sync
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.csv}"
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
echo "  H5AD Single-Cell-Only Pipeline Launcher"
echo "============================================"
echo "Sample file:  $SAMPLE_FILE"
echo "S3 sync:      $(if $SYNC; then echo 'enabled (--sync)'; else echo 'DISABLED (pass --sync to enable)'; fi)"
echo "Output base:  $MEYES_BASE"
echo "MapMyCells:   SKIPPED (not supported for this pipeline's species)"
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

    # Find h5ad file: exact "cell_by_gene.h5ad", or fuzzy match any file ending
    # in "cell_by_gene.h5ad" (e.g. Xenium exports prefixed with the run ID, like
    # "output-XETG00220__.../..._cell_by_gene.h5ad")
    if [ -f "${input_dir}/cell_by_gene.h5ad" ]; then
        h5ad_path="${input_dir}/cell_by_gene.h5ad"
    else
        h5ad_path=$(find "$input_dir" -maxdepth 1 -type f -iname "*cell_by_gene.h5ad" | head -1)
    fi
    if [ -z "$h5ad_path" ] || [ ! -f "$h5ad_path" ]; then
        echo "ERROR: No *cell_by_gene.h5ad file found in: $input_dir"
        errors=$((errors + 1))
        continue
    fi

    output_base="${MEYES_BASE}/${sample_name}"
    meyes_output="${output_base}/meyes_output"
    sm_mapping="${output_base}/sm_output/mapping.json"

    if [ ! -f "$sm_mapping" ]; then
        echo "ERROR: No existing sm_output/mapping.json for ${sample_name} (expected ${sm_mapping})"
        errors=$((errors + 1))
        continue
    fi

    count=$((count + 1))

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  H5AD:       ${h5ad_path}"
    echo "  SC output:  ${meyes_output}"
    echo "  SM mapping: ${sm_mapping} (reused, not regenerated)"

    # Step 1: process_spatial_data (no MMC)
    sc_job=$(sbatch --parsable \
        --job-name="sc_${sample_name}" \
        "${SCRIPT_DIR}/process_h5ad_sc.sbatch" \
        "$h5ad_path" \
        "$meyes_output")
    echo "  [1/3] process_spatial   -> Job ${sc_job}"

    # Step 2: copy existing sm_output/mapping.json into the new meyes_output
    copy_job=$(sbatch --parsable \
        --dependency=afterok:${sc_job} \
        --job-name="cpmap_${sample_name}" \
        --wrap="cp '${sm_mapping}' '${meyes_output}/mapping.json' && echo 'Copied mapping.json from ${sm_mapping} to ${meyes_output}'" \
        --output="/bil/users/ijenie/meyes_process_logs/cpmap_h5ad_sc_${sample_name}_%j.log" \
        --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 --partition=compute)
    echo "  [2/3] copy mapping.json -> Job ${copy_job} (after ${sc_job})"

    # Step 3: s3 sync whole sample dir (size-only, so the already-synced
    # sm_output half is a no-op; only meyes_output actually uploads)
    if $SYNC; then
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${copy_job} \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [3/3] s3_sync           -> Job ${sync_job} (after ${copy_job})"
    else
        echo "  [3/3] s3_sync           -> SKIPPED (pass --sync to enable)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
if [ "$errors" -gt 0 ]; then
    echo "  Skipped ${errors} sample(s) due to missing files"
fi
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
