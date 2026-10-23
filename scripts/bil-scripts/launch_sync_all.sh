#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_sync_all.sh
#
# Submits s3_sync_sample.sbatch for every sample in one or more
# CSV files. Use this to (re-)sync all processed outputs to S3
# without re-running the processing pipeline.
#
# Usage:
#   ./launch_sync_all.sh                         # syncs csv-samples.csv + h5ad-samples.csv
#   ./launch_sync_all.sh my.csv [other.csv ...]  # syncs named CSVs
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CSV_FILES=()
for arg in "$@"; do
    if [ -f "$arg" ]; then
        CSV_FILES+=("$arg")
    else
        echo "WARNING: '$arg' is not a file, skipping"
    fi
done

if [ ${#CSV_FILES[@]} -eq 0 ]; then
    CSV_FILES=("${SCRIPT_DIR}/csv-samples.csv" "${SCRIPT_DIR}/h5ad-samples.csv")
fi

echo "============================================"
echo "  S3 Sync All Samples"
echo "============================================"
echo "CSV files:"
for f in "${CSV_FILES[@]}"; do
    echo "  $f"
done
echo ""

count=0
declare -A seen

for csv in "${CSV_FILES[@]}"; do
    if [ ! -f "$csv" ]; then
        echo "ERROR: $csv not found, skipping"
        continue
    fi

    while IFS=',' read -r sample_name _rest; do
        sample_name="$(echo "$sample_name" | xargs)"
        [[ "$sample_name" =~ ^#.*$ ]] && continue
        [[ -z "$sample_name" ]] && continue

        if [ "${seen[$sample_name]+_}" ]; then
            continue
        fi
        seen[$sample_name]=1

        source_dir="/bil/data/meyes/${sample_name}"
        if [ ! -d "$source_dir" ]; then
            echo "  [SKIP] $sample_name — $source_dir not found on cluster"
            continue
        fi

        count=$((count + 1))
        job_id=$(sbatch --parsable \
            --job-name="sync_${sample_name}" \
            "${SCRIPT_DIR}/s3_sync_sample.sbatch" \
            "$sample_name")
        echo "  [${count}] $sample_name -> Job ${job_id}"

    done < "$csv"
done

echo ""
echo "============================================"
echo "  Submitted ${count} sync job(s)"
echo "  Monitor with: squeue -u \$USER"
echo "  Logs: /bil/users/ijenie/meyes_process_logs/s3_sync_*.log"
echo "============================================"
