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

# Increase AWS S3 transfer concurrency (default is 10). Use a private,
# per-job config file instead of `aws configure set` -- that mutates the
# shared ~/.aws/config, and since many samples' sync jobs run concurrently
# under this account with no dependency between them, concurrent
# read-modify-write on that one file was tearing it (intermittent
# "Unable to parse config file" + a sync that silently no-ops).
export AWS_MAX_CONCURRENT_REQUESTS=50
AWS_CONFIG_FILE="$(mktemp)"
export AWS_CONFIG_FILE
trap 'rm -f "$AWS_CONFIG_FILE"' EXIT
cat > "$AWS_CONFIG_FILE" <<'EOC'
[default]
s3 =
    max_concurrent_requests = 50
    multipart_chunksize = 16MB
EOC

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
# No --size-only: reprocessed files (e.g. manifest.json, de/*.bin.gz) can
# land on the same byte size as the stale S3 copy with different content.
# Default sync also checks mtime, so a freshly regenerated file still
# re-uploads even when its size happens to match.
aws s3 sync "$SOURCE" "$DESTINATION" \
    --exclude "*.csv" \
    --exclude "mmc_output/*"
SYNC_EXIT=$?

if [ "$SYNC_EXIT" -ne 0 ]; then
    echo "ERROR: aws s3 sync failed with exit code $SYNC_EXIT"
    exit "$SYNC_EXIT"
fi

echo "Sync complete at $(date)"
