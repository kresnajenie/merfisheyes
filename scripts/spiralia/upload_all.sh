#!/usr/bin/env bash
# Publish every built labelled-molecule dataset to S3 and print the viewer links.
#
#   ./upload_all.sh [S3_PREFIX] [APP_ORIGIN]
#
# Defaults target the spiralia prefix and the dev site. Data files keep normal
# caching; the two index files that name which data files exist are marked
# no-cache, so replacing a dataset in place can't leave a client asking for a
# file that no longer exists.
set -euo pipefail

SRC=/home/data/yiqun-spiralia/Sep2026/merfisheyes_export
BUCKET=merfisheyes-bil
PREFIX=${1:-yiqun-spiralia}
APP=${2:-https://dev.merfisheyes.com}
REGION=us-west-2
BASE="https://${BUCKET}.s3.${REGION}.amazonaws.com/${PREFIX}"

for dir in "$SRC"/*_lm; do
  [ -d "$dir" ] || continue
  name=$(basename "$dir")
  aws s3 sync "$dir/" "s3://${BUCKET}/${PREFIX}/${name}/" --delete --only-show-errors
  for key in manifest.json obs/metadata.json; do
    aws s3 cp "s3://${BUCKET}/${PREFIX}/${name}/${key}" \
              "s3://${BUCKET}/${PREFIX}/${name}/${key}" \
      --metadata-directive REPLACE --cache-control no-cache \
      --content-type application/json --only-show-errors
  done
  echo "uploaded ${name}"
done

echo
echo "links:"
for dir in "$SRC"/*_lm; do
  [ -d "$dir" ] || continue
  name=$(basename "$dir")
  echo "${APP}/lm-viewer/from-s3?url=${BASE}/${name}"
done
