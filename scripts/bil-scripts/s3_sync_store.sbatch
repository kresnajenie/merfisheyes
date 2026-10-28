#!/bin/bash
#SBATCH --job-name=s3_sync_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/s3_sync_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The finished store -> the bucket the viewer reads by URL.
#
#   sbatch s3_sync_store.sbatch <sample_name> <store_dir>
#
# <store_dir> is the folder holding <name>.zarr, mapping.json and
# stats.json: all three go up, because the viewer needs the store and the
# catalog needs the sidecars.
#
# A store is tens of thousands of small files, so concurrency matters far
# more than chunk size. Set through the environment, never with
# `aws configure set` — that writes ~/.aws/config and two jobs doing it at
# once corrupt each other's (fixed three times already).

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

SAMPLE_NAME="$1"
STORE_DIR="$2"
PREFIX="${MEYES_S3_PREFIX:+${MEYES_S3_PREFIX%/}/}"
DESTINATION="s3://${MEYES_S3_BUCKET}/${PREFIX}${SAMPLE_NAME}/"

meyes_banner "s3 sync"
echo "Sample: $SAMPLE_NAME"
echo "Source: ${STORE_DIR}/"
echo "Dest:   $DESTINATION"

if [ ! -d "$STORE_DIR" ]; then
    echo "ERROR: $STORE_DIR does not exist"
    meyes_done 3
fi
if [ -z "$(find "$STORE_DIR" -type f -print -quit)" ]; then
    echo "ERROR: $STORE_DIR is empty"
    meyes_done 3
fi

module load aws-cli

export AWS_MAX_CONCURRENT_REQUESTS=50
export AWS_S3_MAX_CONCURRENT_REQUESTS=50
export AWS_S3_MULTIPART_CHUNKSIZE=16MB
export AWS_DEFAULT_REGION="$MEYES_S3_REGION"

echo "Files to consider: $(find "$STORE_DIR" -type f | wc -l)"
aws s3 sync "${STORE_DIR}/" "$DESTINATION" --size-only --no-progress
code=$?

if [ $code -eq 0 ]; then
    echo
    echo "The store is at:"
    echo "  https://${MEYES_S3_BUCKET}.s3.${MEYES_S3_REGION}.amazonaws.com/${PREFIX}${SAMPLE_NAME}/"
    echo "The bucket needs public read and CORS for the viewer to open that URL"
    echo "(docs/specs/store-layout.md): a public URL cannot be listed, so the"
    echo "store's consolidated metadata is the only way the viewer learns what"
    echo "is in it."
fi
meyes_done $code
