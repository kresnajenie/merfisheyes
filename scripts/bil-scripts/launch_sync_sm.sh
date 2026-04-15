#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sync_sm.sh
#
# Reads a samples CSV and submits s3_sync_sm jobs ONLY for samples
# whose sm_output/ prefix on S3 is missing or empty. Samples that
# already have sm_output synced are skipped.
#
# Usage:
#   ./launch_sync_sm.sh h5ad-samples.csv
#   ./launch_sync_sm.sh h5ad-samples.csv --dry-run   # list only, no submit
#   ./launch_sync_sm.sh h5ad-samples.csv --force     # submit for every sample
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ $# -lt 1 ]; then
    echo "Usage: ./launch_sync_sm.sh <samples.csv> [--dry-run|--force]"
    exit 1
fi

SAMPLE_FILE="$1"
MODE="${2:-}"

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

S3_BASE="s3://merfisheyes-bil/bil-psc-data2"

# Ensure aws is available on the login node
if ! command -v aws >/dev/null 2>&1; then
    module load aws-cli 2>/dev/null || {
        echo "ERROR: aws CLI not available and 'module load aws-cli' failed"
        exit 1
    }
fi

echo "============================================"
echo "  S3 sm_output Sync Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo "Mode:        ${MODE:-normal}"
echo ""
echo "Scanning S3 for missing sm_output/ prefixes..."
echo ""

checked=0
missing=0
skipped=0
submitted=0
MISSING_SAMPLES=()

while IFS=',' read -r sample_name _rest; do
    sample_name="$(echo "${sample_name:-}" | xargs)"

    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    checked=$((checked + 1))

    if [ "$MODE" = "--force" ]; then
        MISSING_SAMPLES+=("$sample_name")
        missing=$((missing + 1))
        printf "  [check] %-20s FORCE\n" "$sample_name"
        continue
    fi

    s3_prefix="${S3_BASE}/${sample_name}/sm_output/"

    if aws s3 ls "$s3_prefix" 2>/dev/null | grep -q .; then
        skipped=$((skipped + 1))
        printf "  [check] %-20s already synced, skipping\n" "$sample_name"
    else
        missing=$((missing + 1))
        MISSING_SAMPLES+=("$sample_name")
        printf "  [check] %-20s MISSING sm_output\n" "$sample_name"
    fi

done < "$SAMPLE_FILE"

echo ""
echo "Scan summary: $checked checked, $missing missing, $skipped already synced"
echo ""

if [ "$missing" -eq 0 ]; then
    echo "Nothing to sync. Exiting."
    exit 0
fi

if [ "$MODE" = "--dry-run" ]; then
    echo "Dry run — not submitting. Samples that would be synced:"
    for s in "${MISSING_SAMPLES[@]}"; do
        echo "  - $s"
    done
    exit 0
fi

echo "Submitting sm_output sync jobs..."
echo ""

for sample_name in "${MISSING_SAMPLES[@]}"; do
    submitted=$((submitted + 1))
    sync_job=$(sbatch --parsable \
        --job-name="sync_sm_${sample_name}" \
        "${SCRIPT_DIR}/s3_sync_sm.sbatch" \
        "$sample_name")
    echo "  [${submitted}/${missing}] ${sample_name} -> Job ${sync_job}"
done

echo ""
echo "============================================"
echo "  Submitted ${submitted} sm_output sync job(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
