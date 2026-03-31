#!/usr/bin/env bash
# Sync processed single molecule datasets to S3
#
# Usage:
#   ./scripts/sync_to_s3.sh                  # sync only new/changed files
#   ./scripts/sync_to_s3.sh --overwrite      # re-upload everything

set -euo pipefail

LOCAL_DIR="/data/merfisheyes-test-data/test-data-sizes-run/single-molecule"
S3_BUCKET="s3://merfisheyes-bil/test-data-sizes-run"

OVERWRITE=false
if [[ "${1:-}" == "--overwrite" ]]; then
  OVERWRITE=true
fi

if [[ ! -d "$LOCAL_DIR" ]]; then
  echo "Error: Local directory not found: $LOCAL_DIR"
  exit 1
fi

# List all dataset folders
datasets=("$LOCAL_DIR"/*)

echo "Syncing ${#datasets[@]} datasets to $S3_BUCKET/"
echo "Mode: $(if $OVERWRITE; then echo 'overwrite (cp --recursive)'; else echo 'sync (new/changed only)'; fi)"
echo ""

for dataset_dir in "${datasets[@]}"; do
  if [[ ! -d "$dataset_dir" ]]; then
    continue
  fi

  folder_name=$(basename "$dataset_dir")
  s3_path="$S3_BUCKET/$folder_name/"

  echo "--- $folder_name ---"

  if $OVERWRITE; then
    aws s3 cp "$dataset_dir/" "$s3_path" --recursive
  else
    aws s3 sync "$dataset_dir/" "$s3_path"
  fi

  echo ""
done

echo "Done."
