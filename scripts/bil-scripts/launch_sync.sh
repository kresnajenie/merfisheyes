#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sync.sh
#
# Reads a samples CSV and submits s3_sync_sample jobs only.
# For re-syncing samples that already have processed output on BIL.
#
# Usage:
#   ./launch_sync.sh sync-remaining.csv
#   ./launch_sync.sh                      # uses sync-remaining.csv in same dir
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMPLE_FILE="$1"

if [ -z "$SAMPLE_FILE" ]; then
    echo "Usage: ./launch_sync.sh <samples.csv>"
    exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  S3 Sync Launcher"
echo "============================================"
echo "Sample file: $SAMPLE_FILE"
echo ""

count=0

while IFS=',' read -r sample_name _rest; do
    sample_name="$(echo "$sample_name" | xargs)"

    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))

    sync_job=$(sbatch --parsable \
        --job-name="sync_${sample_name}" \
        "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
        "$sample_name")
    echo "  [${count}] ${sample_name} -> Job ${sync_job}"

done < "$SAMPLE_FILE"

echo ""
echo "============================================"
echo "  Submitted ${count} sync job(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
