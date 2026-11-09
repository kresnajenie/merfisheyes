#!/bin/bash
#SBATCH --job-name=register_lamin
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/register_lamin_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The synced store -> a record in the Lamin catalog, at the S3 path it
# already has (nothing is copied), then checked: the record's path, a fresh
# S3 listing against the record's file count and bytes, the public URL with
# CORS, and a read back through Lamin. Non-zero when any check fails.
#
#   sbatch register_lamin.sbatch <sample_name> [register_lamin.py options]
#
# Runs after s3_sync_store.sbatch (launch_spatialdata_pipeline.sh --sync --lamin).
# Needs $MEYES_LAMIN_PYTHON with lamindb, `lamin login` done once on the
# cluster, and AWS credentials that can list the prefix.

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

SAMPLE_NAME="$1"; shift

meyes_banner "lamin register"
echo "Sample:   $SAMPLE_NAME"
echo "Where:    s3://${MEYES_S3_BUCKET}/${MEYES_S3_PREFIX}/${SAMPLE_NAME}/"
echo "Catalog:  $MEYES_LAMIN_INSTANCE"

module load aws-cli 2>/dev/null || true
export AWS_DEFAULT_REGION="$MEYES_S3_REGION"

"$MEYES_LAMIN_PYTHON" -I "${MEYES_BIL_SCRIPTS:-$MEYES_REPO/scripts/bil-scripts}/register_lamin.py" \
    "$SAMPLE_NAME" \
    --format spatialdata \
    --bucket "$MEYES_S3_BUCKET" \
    --prefix "$MEYES_S3_PREFIX" \
    --region "$MEYES_S3_REGION" \
    --instance "$MEYES_LAMIN_INSTANCE" \
    "$@"
meyes_done $?
