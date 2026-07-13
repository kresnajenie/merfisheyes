#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sync.sh
#
# Runs s3 sync ONLY (no process_spatial, no MMC) — for re-syncing datasets
# whose upload didn't complete or was interrupted. Submits
# s3_sync_sample.sbatch per sample with no dependency, so it's safe to run
# even while other jobs for the same sample are unrelated/finished.
#
# s3_sync_sample.sbatch compares size + mtime (no --size-only, removed
# because reprocessed files can land on the same byte size as a stale S3
# copy with different content), so re-running only re-uploads files that
# are actually missing or newer.
#
# Usage:
#   ./launch_sync.sh ace-dip-use     # single sample
#   ./launch_sync.sh samples.csv     # from file (only sample_name column used)
#   ./launch_sync.sh                 # uses samples.csv in same dir
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
echo "  S3 Sync (sync only)"
echo "============================================"
echo "Source: $SAMPLE_FILE"
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
    echo "── Sample ${count}: ${sample_name} -> Job ${sync_job}"

done < "$SAMPLE_FILE"

echo ""
echo "============================================"
echo "  Submitted ${count} sync job(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
