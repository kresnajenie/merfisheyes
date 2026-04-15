#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_h5ad_pipeline.sh
#
# Reads samples.csv (sample_name,input_dir) and for each sample submits:
#   1. map_my_cell            (h5ad → mmc_output/mapping_output.csv)
#   2. process_single_molecule (csv → sm_output, runs in parallel with step 1)
#   3. process_spatial_data    (h5ad + mmc csv → meyes_output, after step 1)
#   4. copy mapping.json       (sm_output → meyes_output, after steps 2+3)
#   5. s3 sync                 (both meyes_output + sm_output, after step 4)
#
# Expected input directory structure:
#   {input_dir}/cell_by_gene.h5ad
#   {input_dir}/segmented_spot_table.csv
#
# Usage:
#   ./launch_h5ad_pipeline.sh                  # uses samples.csv in same dir
#   ./launch_h5ad_pipeline.sh my_samples.csv   # custom sample list
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="${1:-${SCRIPT_DIR}/samples.csv}"
MEYES_BASE="/bil/data/meyes"
REFERENCE_DIR="/bil/data/meyes/mapmycells-reference"
S3_BUCKET="merfisheyes-bil"
S3_PREFIX_BASE="bil-psc-data2"
S3_HTTPS_BASE="https://${S3_BUCKET}.s3.us-west-2.amazonaws.com/${S3_PREFIX_BASE}"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  H5AD Pipeline Launcher"
echo "============================================"
echo "Sample file:  $SAMPLE_FILE"
echo "Output base:  $MEYES_BASE"
echo "Reference:    $REFERENCE_DIR"
echo "S3 base:      $S3_HTTPS_BASE"
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

    # Validate input files
    h5ad_path="${input_dir}/cell_by_gene.h5ad"
    csv_path="${input_dir}/segmented_spot_table.csv"

    if [ ! -f "$h5ad_path" ]; then
        echo "ERROR: cell_by_gene.h5ad not found: $h5ad_path"
        errors=$((errors + 1))
        continue
    fi
    if [ ! -f "$csv_path" ]; then
        echo "ERROR: segmented_spot_table.csv not found: $csv_path"
        errors=$((errors + 1))
        continue
    fi

    count=$((count + 1))

    # Output paths
    output_base="${MEYES_BASE}/${sample_name}"
    mmc_output="${output_base}/mmc_output"
    meyes_output="${output_base}/meyes_output"
    sm_output="${output_base}/sm_output"
    sm_s3_url="${S3_HTTPS_BASE}/${sample_name}/sm_output"

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input dir:    ${input_dir}"
    echo "  H5AD:         ${h5ad_path}"
    echo "  CSV:          ${csv_path}"
    echo "  MMC output:   ${mmc_output}"
    echo "  SC output:    ${meyes_output}"
    echo "  SM output:    ${sm_output}"
    echo "  SM S3 URL:    ${sm_s3_url}"

    # Step 1: map_my_cell (h5ad → mmc_output)
    mmc_job=$(sbatch --parsable \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$h5ad_path" \
        "$mmc_output" \
        "$REFERENCE_DIR")
    echo "  [1/5] map_my_cell      -> Job ${mmc_job}"

    # Step 2: process_single_molecule (runs in parallel with step 1)
    sm_job=$(sbatch --parsable \
        --job-name="sm_${sample_name}" \
        "${SCRIPT_DIR}/process_h5ad_sm.sbatch" \
        "$csv_path" \
        "$sm_output" \
        "$sm_s3_url")
    echo "  [2/5] process_sm       -> Job ${sm_job} (parallel)"

    # Step 3: process_spatial_data (after map_my_cell finishes)
    sc_job=$(sbatch --parsable \
        --dependency=afterok:${mmc_job} \
        --job-name="sc_${sample_name}" \
        "${SCRIPT_DIR}/process_h5ad_sc.sbatch" \
        "$h5ad_path" \
        "$meyes_output" \
        "${mmc_output}/mapping_output.csv")
    echo "  [3/5] process_spatial  -> Job ${sc_job} (after ${mmc_job})"

    # Step 4: copy mapping.json from sm_output to meyes_output (after both SC and SM)
    copy_job=$(sbatch --parsable \
        --dependency=afterok:${sc_job}:${sm_job} \
        --job-name="cpmap_${sample_name}" \
        --wrap="cp '${sm_output}/mapping.json' '${meyes_output}/mapping.json' && echo 'Copied mapping.json from ${sm_output} to ${meyes_output}'" \
        --output="/bil/users/ijenie/meyes_process_logs/cpmap_h5ad_${sample_name}_%j.log" \
        --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 --partition=compute)
    echo "  [4/5] copy mapping.json -> Job ${copy_job} (after ${sc_job} + ${sm_job})"

    # Step 5: s3 sync both meyes_output and sm_output (after copy)
    sync_job=$(sbatch --parsable \
        --dependency=afterok:${copy_job} \
        --job-name="sync_${sample_name}" \
        --output="/bil/users/ijenie/meyes_process_logs/s3_sync_h5ad_${sample_name}_%j.log" \
        --ntasks=1 --cpus-per-task=8 --mem=4G --time=2-00:00:00 --partition=compute \
        --wrap="
module load aws-cli
export AWS_MAX_CONCURRENT_REQUESTS=50
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.multipart_chunksize 16MB

echo '=== S3 sync for ${sample_name} ==='
echo 'Started at \$(date)'

echo ''
echo '--- Syncing meyes_output ---'
aws s3 sync '${meyes_output}/' 's3://${S3_BUCKET}/${S3_PREFIX_BASE}/${sample_name}/meyes_output/' --size-only
echo 'meyes_output sync done at \$(date)'

echo ''
echo '--- Syncing sm_output ---'
aws s3 sync '${sm_output}/' 's3://${S3_BUCKET}/${S3_PREFIX_BASE}/${sample_name}/sm_output/' --size-only --exclude '*.csv'
echo 'sm_output sync done at \$(date)'

echo ''
echo 'All syncs complete at \$(date)'
")
    echo "  [5/5] s3_sync          -> Job ${sync_job} (after ${copy_job})"

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
if [ "$errors" -gt 0 ]; then
    echo "  Skipped ${errors} sample(s) due to missing files"
fi
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
