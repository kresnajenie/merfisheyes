#!/bin/bash
#SBATCH --job-name=s3_sync_sm
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/s3_sync_sm_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch s3_sync_sm.sbatch SAMPLE_NAME [--overwrite]"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME   Sample name (e.g. ace-dip-use)"
    echo "  --overwrite   (optional) Re-upload everything instead of syncing"
    echo ""
    echo "Paths:"
    echo "  Source: /bil/data/meyes/{SAMPLE_NAME}/sm_output/"
    echo "  Dest:   s3://merfisheyes-bil/bil-psc-data2/{SAMPLE_NAME}/sm_output/"
    echo ""
    echo "Examples:"
    echo "  sbatch s3_sync_sm.sbatch ace-dip-use"
    echo "  sbatch s3_sync_sm.sbatch ace-dip-use --overwrite"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
SAMPLE_NAME="$1"
SOURCE="/bil/data/meyes/${SAMPLE_NAME}/sm_output/"
DESTINATION="s3://merfisheyes-bil/bil-psc-data2/${SAMPLE_NAME}/sm_output/"

OVERWRITE=false
if [ "${2:-}" = "--overwrite" ]; then
    OVERWRITE=true
fi

echo "=== s3_sync_sm ==="
echo "Sample: $SAMPLE_NAME"
echo "Source: $SOURCE"
echo "Dest:   $DESTINATION"
echo "Mode:   $(if $OVERWRITE; then echo 'overwrite (cp --recursive)'; else echo 'sync (new/changed only)'; fi)"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load aws-cli

# ── Preflight check ─────────────────────────────────────────
if [ ! -d "$SOURCE" ]; then
    echo "ERROR: Source directory $SOURCE does not exist"
    exit 1
fi

# ── Sync ────────────────────────────────────────────────────
if $OVERWRITE; then
    aws s3 cp "$SOURCE" "$DESTINATION" --recursive
else
    aws s3 sync "$SOURCE" "$DESTINATION" --size-only
fi

echo ""
echo "Sync complete at $(date)"
