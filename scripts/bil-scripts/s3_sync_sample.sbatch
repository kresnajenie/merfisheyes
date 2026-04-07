#!/bin/bash
#SBATCH --job-name=s3_sync_sample
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/s3_sync_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Arguments ────────────────────────────────────────────────
SAMPLE_NAME="$1"   # e.g. ace-dud-vex

SOURCE="/bil/data/meyes/${SAMPLE_NAME}/"
DESTINATION="s3://merfisheyes-bil/bil-psc-data2/${SAMPLE_NAME}/"

echo "=== s3_sync ==="
echo "Sample: $SAMPLE_NAME"
echo "Source: $SOURCE"
echo "Dest:   $DESTINATION"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load aws-cli

# Increase AWS S3 transfer concurrency (default is 10)
export AWS_MAX_CONCURRENT_REQUESTS=50
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.multipart_chunksize 16MB

# ── Preflight check ─────────────────────────────────────────
if [ ! -d "$SOURCE" ]; then
    echo "ERROR: Source directory $SOURCE does not exist"
    exit 1
fi

FILE_COUNT=$(find "$SOURCE" -type f | head -100 | wc -l)
if [ "$FILE_COUNT" -eq 0 ]; then
    echo "ERROR: Source directory appears empty"
    exit 1
fi

# ── Sync ─────────────────────────────────────────────────────
aws s3 sync "$SOURCE" "$DESTINATION" \
    --size-only \
    --exclude "*.csv" \
    --exclude "mmc_output/*"

echo "Sync complete at $(date)"
