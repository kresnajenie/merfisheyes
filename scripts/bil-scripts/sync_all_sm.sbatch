#!/bin/bash
#SBATCH --job-name=sync_all_sm
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/sync_all_sm_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ═══════════════════════════════════════════════════════════════
# sync_all_sm.sbatch
#
# One-shot sync of every sample's sm_output/ directory under
# /bil/data/meyes/ to s3://merfisheyes-bil/bil-psc-data2/.
#
# `aws s3 sync` is idempotent — only missing/changed files are
# uploaded, so samples that are already fully synced are skipped
# automatically. No per-sample bookkeeping needed.
#
# Usage:
#   sbatch sync_all_sm.sbatch
# ═══════════════════════════════════════════════════════════════

SOURCE="/bil/data/meyes/"
DESTINATION="s3://merfisheyes-bil/bil-psc-data2/"

echo "=== sync_all_sm ==="
echo "Source: $SOURCE"
echo "Dest:   $DESTINATION"
echo "Filter: */sm_output/*"
echo "Started at $(date)"

module load aws-cli

export AWS_MAX_CONCURRENT_REQUESTS=50
export AWS_S3_MAX_CONCURRENT_REQUESTS=50
export AWS_S3_MULTIPART_CHUNKSIZE=16MB

aws s3 sync "$SOURCE" "$DESTINATION" \
    --exclude "*" \
    --include "*/sm_output/*" \
    --size-only

echo ""
echo "Sync complete at $(date)"
