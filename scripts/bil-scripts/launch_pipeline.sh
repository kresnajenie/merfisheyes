#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_pipeline.sh
#
# End-to-end MERSCOPE multi-slice pipeline.
# For each sample: combine → map_my_cell → process_spatial → s3_sync
#
# Usage:
#   CSV mode:
#     ./launch_pipeline.sh                                              # uses samples.csv
#     ./launch_pipeline.sh my_samples.csv                               # custom CSV
#     ./launch_pipeline.sh my_samples.csv /output/base                  # custom output base
#     ./launch_pipeline.sh my_samples.csv /output/base human            # + species
#
#   Single sample mode (input_path, output_path):
#     ./launch_pipeline.sh /path/to/input /path/to/output               # output is absolute
#     ./launch_pipeline.sh /path/to/input /path/to/output human         # + species
#
#   S3 sync (disabled by default):
#     ./launch_pipeline.sh samples.csv --sync s3://bucket/prefix
#     ./launch_pipeline.sh samples.csv --sync https://bucket.s3.us-west-2.amazonaws.com/prefix
#     ./launch_pipeline.sh samples.csv --sync s3://bucket/prefix --region us-east-1
#
# Output:
#   CSV mode:    {output_base}/{sample_name}/  (default output_base: /bil/data/meyes)
#   Single mode: exactly the output_path you provide
#
# MapMyCells Reference Data:
#   Requires reference taxonomy data at /bil/data/meyes/mapmycells-reference
#   See README.md for download instructions.
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"
SYNC=false
SYNC_PREFIX=""
REGION="us-west-2"
SPECIES="mouse"

# Parse --sync and --region flags
ARGS=("$@")
i=0
while [ $i -lt ${#ARGS[@]} ]; do
    case "${ARGS[$i]}" in
        --sync)
            SYNC=true
            i=$((i + 1))
            if [ $i -lt ${#ARGS[@]} ] && [[ "${ARGS[$i]}" != "--"* ]]; then
                SYNC_PREFIX="${ARGS[$i]}"
            else
                echo "ERROR: --sync requires an S3 prefix argument"
                echo "  Example: --sync s3://bucket/prefix"
                echo "  Example: --sync https://bucket.s3.us-west-2.amazonaws.com/prefix"
                exit 1
            fi
            ;;
        --region)
            i=$((i + 1))
            if [ $i -lt ${#ARGS[@]} ]; then
                REGION="${ARGS[$i]}"
            fi
            ;;
    esac
    i=$((i + 1))
done

# Convert S3 prefix between s3:// and https:// formats
S3_URL=""
HTTPS_URL=""
if $SYNC; then
    SYNC_PREFIX="${SYNC_PREFIX%/}"  # strip trailing slash
    if [[ "$SYNC_PREFIX" == s3://* ]]; then
        # s3://bucket/prefix → https://bucket.s3.region.amazonaws.com/prefix
        S3_URL="$SYNC_PREFIX"
        BUCKET="$(echo "$SYNC_PREFIX" | sed 's|s3://||' | cut -d/ -f1)"
        PREFIX="$(echo "$SYNC_PREFIX" | sed 's|s3://[^/]*/||')"
        HTTPS_URL="https://${BUCKET}.s3.${REGION}.amazonaws.com/${PREFIX}"
    elif [[ "$SYNC_PREFIX" == https://* ]]; then
        # https://bucket.s3.region.amazonaws.com/prefix → s3://bucket/prefix
        HTTPS_URL="$SYNC_PREFIX"
        BUCKET="$(echo "$SYNC_PREFIX" | sed 's|https://||' | cut -d. -f1)"
        PREFIX="$(echo "$SYNC_PREFIX" | sed "s|https://[^/]*/||")"
        REGION="$(echo "$SYNC_PREFIX" | sed 's|https://[^.]*\.s3\.||' | cut -d. -f1)"
        S3_URL="s3://${BUCKET}/${PREFIX}"
    else
        echo "ERROR: --sync prefix must start with s3:// or https://"
        exit 1
    fi
fi

# Collect positional args (skip flags and their values)
POSITIONAL=()
skip_next=false
for arg in "$@"; do
    if $skip_next; then
        skip_next=false
        continue
    fi
    if [ "$arg" = "--sync" ] || [ "$arg" = "--region" ]; then
        skip_next=true
        continue
    fi
    [[ "$arg" == "--"* ]] && continue
    POSITIONAL+=("$arg")
done

# Determine mode based on positional args
if [ ${#POSITIONAL[@]} -eq 0 ]; then
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
elif [ ${#POSITIONAL[@]} -eq 1 ] && [ -f "${POSITIONAL[0]}" ]; then
    SAMPLE_FILE="${POSITIONAL[0]}"
elif [ ${#POSITIONAL[@]} -ge 2 ] && [ -f "${POSITIONAL[0]}" ]; then
    SAMPLE_FILE="${POSITIONAL[0]}"
    if [[ "${POSITIONAL[1]}" == /* ]]; then
        MEYES_BASE="${POSITIONAL[1]}"
        [ ${#POSITIONAL[@]} -ge 3 ] && SPECIES="${POSITIONAL[2]}"
    else
        SPECIES="${POSITIONAL[1]}"
    fi
elif [ ${#POSITIONAL[@]} -ge 2 ] && [ ! -f "${POSITIONAL[0]}" ]; then
    # Single sample mode: input_path output_path [species]
    INPUT_PATH="${POSITIONAL[0]}"
    OUTPUT_PATH="${POSITIONAL[1]}"
    SAMPLE_NAME="$(basename "$OUTPUT_PATH")"
    MEYES_BASE="$(dirname "$OUTPUT_PATH")"
    SAMPLE_FILE=$(mktemp)
    echo "${SAMPLE_NAME},${INPUT_PATH}" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
    [ ${#POSITIONAL[@]} -ge 3 ] && SPECIES="${POSITIONAL[2]}"
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    echo ""
    echo "Usage:"
    echo "  $0                                                      # uses samples.csv"
    echo "  $0 my_samples.csv                                       # custom CSV"
    echo "  $0 my_samples.csv /output/base                          # + output base"
    echo "  $0 /path/to/input /path/to/output                       # single sample"
    echo "  $0 /path/to/input /path/to/output human                 # + species"
    echo "  $0 samples.csv --sync s3://bucket/prefix                # + S3 sync"
    echo "  $0 samples.csv --sync s3://bucket/prefix --region us-east-1"
    exit 1
fi

echo "============================================"
echo "  MERFISH Eyes Pipeline Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "Species:     $SPECIES"
echo "Output base: $MEYES_BASE"
if $SYNC; then
    echo "S3 sync:     ENABLED"
    echo "  S3 URL:    $S3_URL"
    echo "  HTTPS URL: $HTTPS_URL"
    echo "  Region:    $REGION"
else
    echo "S3 sync:     DISABLED (pass --sync <prefix> to enable)"
fi
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

    # Step 1: combine_slices (combines multi-slice MERSCOPE data)
    combine_job=$(sbatch --parsable \
        --job-name="combine_${sample_name}" \
        "${SCRIPT_DIR}/combine_slices.sbatch" \
        "$input_path" \
        "$output_base")
    echo "  [1/4] combine_slices  -> Job ${combine_job}"

    # Step 2: map_my_cell (waits for combine_slices)
    mmc_job=$(sbatch --parsable \
        --dependency=afterok:${combine_job} \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$combined_output" \
        "$mmc_output" \
        "$SPECIES")
    echo "  [2/4] map_my_cell     -> Job ${mmc_job} (after ${combine_job})"

    # Step 3: process_spatial (waits for map_my_cell)
    process_job=$(sbatch --parsable \
        --dependency=afterok:${mmc_job} \
        --job-name="process_${sample_name}" \
        "${SCRIPT_DIR}/process_spatial.sbatch" \
        "$combined_output" \
        "$meyes_output" \
        "$mmc_csv")
    echo "  [3/4] process_spatial -> Job ${process_job} (after ${mmc_job})"

    # Step 4: s3 sync (waits for process_spatial)
    if $SYNC; then
        s3_dest="${S3_URL}/${sample_name}/"
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${process_job} \
            --job-name="sync_${sample_name}" \
            --output="/bil/users/ijenie/meyes_process_logs/s3_sync_${sample_name}_%j.log" \
            --ntasks=1 --cpus-per-task=1 --mem=4G --time=2-00:00:00 --partition=compute \
            --mail-type=BEGIN,END,FAIL --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu \
            --wrap="
module load aws-cli
aws configure set default.s3.max_concurrent_requests 50
echo '=== S3 sync: ${sample_name} ==='
echo 'Source: ${output_base}/'
echo 'Dest:   ${s3_dest}'
echo 'Started at \$(date)'
aws s3 sync '${output_base}/' '${s3_dest}' --size-only --exclude '*.csv' --exclude 'mmc_output/*'
echo 'Sync complete at \$(date)'
")
        echo "  [4/4] s3_sync         -> Job ${sync_job} (after ${process_job})"
        echo "        Dest: ${s3_dest}"
    else
        echo "  [4/4] s3_sync         -> SKIPPED (pass --sync <prefix> to enable)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
