#!/bin/bash
#SBATCH --job-name=s3_sync_sm_output
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/s3_sync_sm_output_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch s3_sync_sm_output.sbatch SOURCE_DIR S3_DEST"
    echo ""
    echo "Flat recursive 'aws s3 sync' of one sample's sm_output/ directory."
    echo "Run as a real sbatch script (not --wrap) so 'module load aws-cli' works —"
    echo "sbatch --wrap runs under /bin/sh, where the module function isn't defined."
    echo ""
    echo "Examples:"
    echo "  sbatch s3_sync_sm_output.sbatch /bil/data/meyes/ace-pie-bit/sm_output/ s3://merfisheyes-bil/bil-psc-data2/ace-pie-bit/sm_output/"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
SOURCE="$1"
DESTINATION="$2"

echo "=== s3_sync_sm_output ==="
echo "Source: $SOURCE"
echo "Dest:   $DESTINATION"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load aws-cli
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.multipart_chunksize 16MB

# ── Run ──────────────────────────────────────────────────────
aws s3 sync "$SOURCE" "$DESTINATION" --size-only

echo "Sync complete at $(date)"
