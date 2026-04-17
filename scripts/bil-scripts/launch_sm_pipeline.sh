#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sm_pipeline.sh
#
# Single molecule pipeline.
# For each sample: process_single_molecule → copy mapping.json → s3_sync
#
# Usage:
#   CSV mode:
#     ./launch_sm_pipeline.sh                                           # uses samples.csv
#     ./launch_sm_pipeline.sh my_samples.csv                            # custom CSV
#     ./launch_sm_pipeline.sh my_samples.csv /output/base               # custom output base
#
#   Single sample mode (input_path, output_path):
#     ./launch_sm_pipeline.sh /path/to/input /path/to/output            # output is absolute
#
#   S3 sync (disabled by default):
#     ./launch_sm_pipeline.sh samples.csv --sync s3://bucket/prefix
#     ./launch_sm_pipeline.sh samples.csv --sync https://bucket.s3.us-west-2.amazonaws.com/prefix
#     ./launch_sm_pipeline.sh samples.csv --sync s3://bucket/prefix --region us-east-1
#
# Output:
#   CSV mode:    {output_base}/{sample_name}/sm_output/
#   Single mode: {output_path}/sm_output/
#
# When --sync is provided:
#   - mapping.json contains full S3 URLs for right-click SM viewer linking
#   - mapping.json is copied to meyes_output/ for SC viewer integration
#   - sm_output/ is uploaded to S3
# When --sync is NOT provided:
#   - mapping.json contains relative paths (can be updated later with update_mapping_prefix.py)
#   - No S3 upload, no copy to meyes_output
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"
SYNC=false
SYNC_PREFIX=""
REGION="us-west-2"

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
        S3_URL="$SYNC_PREFIX"
        BUCKET="$(echo "$SYNC_PREFIX" | sed 's|s3://||' | cut -d/ -f1)"
        PREFIX="$(echo "$SYNC_PREFIX" | sed 's|s3://[^/]*/||')"
        HTTPS_URL="https://${BUCKET}.s3.${REGION}.amazonaws.com/${PREFIX}"
    elif [[ "$SYNC_PREFIX" == https://* ]]; then
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
    fi
elif [ ${#POSITIONAL[@]} -ge 2 ] && [ ! -f "${POSITIONAL[0]}" ]; then
    # Single sample mode: input_path output_path
    INPUT_PATH="${POSITIONAL[0]}"
    OUTPUT_PATH="${POSITIONAL[1]}"
    SAMPLE_NAME="$(basename "$OUTPUT_PATH")"
    MEYES_BASE="$(dirname "$OUTPUT_PATH")"
    SAMPLE_FILE=$(mktemp)
    echo "${SAMPLE_NAME},${INPUT_PATH}" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    echo ""
    echo "Usage:"
    echo "  $0                                                      # uses samples.csv"
    echo "  $0 my_samples.csv                                       # custom CSV"
    echo "  $0 my_samples.csv /output/base                          # + output base"
    echo "  $0 /path/to/input /path/to/output                       # single sample"
    echo "  $0 samples.csv --sync s3://bucket/prefix                # + S3 sync"
    echo "  $0 samples.csv --sync s3://bucket/prefix --region us-east-1"
    exit 1
fi

echo "============================================"
echo "  Single Molecule Pipeline Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "Output base: $MEYES_BASE"
if $SYNC; then
    echo "S3 sync:     ENABLED"
    echo "  S3 URL:    $S3_URL"
    echo "  HTTPS URL: $HTTPS_URL"
    echo "  Region:    $REGION"
else
    echo "S3 sync:     DISABLED (pass --sync <prefix> to enable)"
    echo "  mapping.json will use relative paths"
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
    sm_output="${output_base}/sm_output"
    meyes_output="${output_base}/meyes_output"

    # Determine S3 prefix for mapping.json
    if $SYNC; then
        s3_prefix="${HTTPS_URL}/${sample_name}/sm_output"
    else
        # Relative paths — can be updated later with update_mapping_prefix.py
        s3_prefix=""
    fi

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:      ${input_path}"
    echo "  SM output:  ${sm_output}"
    if $SYNC; then
        echo "  S3 prefix:  ${s3_prefix}"
    fi

    # Step 1: process_single_molecule
    SM_CMD_ARGS=("$input_path" "$sm_output")
    if [ -n "$s3_prefix" ]; then
        SM_CMD_ARGS+=(--s3-prefix "$s3_prefix")
    fi
    SM_CMD_ARGS+=(--sample-workers 4 --workers 8)

    process_job=$(sbatch --parsable \
        --job-name="sm_${sample_name}" \
        "${SCRIPT_DIR}/process_single_molecule.sbatch" \
        "${SM_CMD_ARGS[@]}")
    echo "  [1/3] process_sm      -> Job ${process_job}"

    if $SYNC; then
        # Step 2: copy mapping.json to meyes_output (after process completes)
        copy_job=$(sbatch --parsable \
            --dependency=afterok:${process_job} \
            --job-name="cpmap_${sample_name}" \
            --wrap="mkdir -p '${meyes_output}' && cp '${sm_output}/mapping.json' '${meyes_output}/mapping.json' && echo 'Copied mapping.json to ${meyes_output}'" \
            --output="/bil/users/ijenie/meyes_process_logs/cpmap_${sample_name}_%j.log" \
            --ntasks=1 --cpus-per-task=1 --mem=1G --time=00:05:00 --partition=compute)
        echo "  [2/3] copy mapping    -> Job ${copy_job} (after ${process_job})"

        # Step 3: s3 sync (waits for copy)
        s3_dest="${S3_URL}/${sample_name}/sm_output/"
        sync_job=$(sbatch --parsable \
            --dependency=afterok:${copy_job} \
            --job-name="sync_sm_${sample_name}" \
            --output="/bil/users/ijenie/meyes_process_logs/s3_sync_sm_${sample_name}_%j.log" \
            --ntasks=1 --cpus-per-task=8 --mem=4G --time=2-00:00:00 --partition=compute \
            --mail-type=BEGIN,END,FAIL --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu \
            --wrap="
module load aws-cli
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.multipart_chunksize 16MB
echo '=== SM S3 sync: ${sample_name} ==='
echo 'Source: ${sm_output}/'
echo 'Dest:   ${s3_dest}'
echo 'Started at \$(date)'
aws s3 sync '${sm_output}/' '${s3_dest}' --size-only
echo 'Sync complete at \$(date)'
")
        echo "  [3/3] s3_sync_sm      -> Job ${sync_job} (after ${copy_job})"
        echo "        Dest: ${s3_dest}"
    else
        echo "  [2/3] copy mapping    -> SKIPPED (no --sync)"
        echo "  [3/3] s3_sync_sm      -> SKIPPED (pass --sync <prefix> to enable)"
    fi

    echo ""

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
