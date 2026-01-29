#!/bin/bash
#
# Sync Dataset from Production to Staging
#
# This script copies a dataset from production S3 bucket to staging S3 bucket.
# Use this to test production datasets in staging environment.
#
# ⚠️  ONE-WAY ONLY: Production → Staging (never staging → production)
#
# Usage:
#   ./scripts/sync-prod-to-staging.sh <dataset-id>
#
# Example:
#   ./scripts/sync-prod-to-staging.sh cm2abc123xyz456
#

set -e  # Exit on error

# Configuration
PROD_BUCKET="merfisheyes-production"
STAGING_BUCKET="merfisheyes-staging"
AWS_REGION="us-east-1"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
info() {
  echo -e "${BLUE}ℹ${NC} $1"
}

success() {
  echo -e "${GREEN}✓${NC} $1"
}

warning() {
  echo -e "${YELLOW}⚠${NC} $1"
}

error() {
  echo -e "${RED}✗${NC} $1"
}

# Check arguments
if [ -z "$1" ]; then
  error "Missing dataset ID"
  echo ""
  echo "Usage: $0 <dataset-id>"
  echo ""
  echo "Example:"
  echo "  $0 cm2abc123xyz456"
  echo ""
  echo "This will copy the dataset from:"
  echo "  s3://$PROD_BUCKET/datasets/<dataset-id>"
  echo "To:"
  echo "  s3://$STAGING_BUCKET/datasets/<dataset-id>"
  echo ""
  exit 1
fi

DATASET_ID=$1

# Display sync plan
echo ""
info "Dataset Sync: Production → Staging"
echo ""
echo "  Dataset ID: $DATASET_ID"
echo "  From:       s3://$PROD_BUCKET/datasets/$DATASET_ID"
echo "  To:         s3://$STAGING_BUCKET/datasets/$DATASET_ID"
echo "  Region:     $AWS_REGION"
echo ""

# Check if AWS CLI is installed
if ! command -v aws &> /dev/null; then
  error "AWS CLI not found"
  echo ""
  echo "Install AWS CLI:"
  echo "  brew install awscli              # macOS"
  echo "  apt-get install awscli           # Ubuntu/Debian"
  echo "  pip install awscli               # pip"
  echo ""
  exit 1
fi

# Check if dataset exists in production
info "Checking if dataset exists in production..."

if ! aws s3 ls "s3://$PROD_BUCKET/datasets/$DATASET_ID/" --region $AWS_REGION > /dev/null 2>&1; then
  error "Dataset not found in production bucket"
  echo ""
  echo "Verify dataset ID and try again."
  echo ""
  echo "To list all production datasets:"
  echo "  aws s3 ls s3://$PROD_BUCKET/datasets/"
  echo ""
  exit 1
fi

success "Dataset found in production"

# Get dataset size
info "Calculating dataset size..."
DATASET_SIZE=$(aws s3 ls "s3://$PROD_BUCKET/datasets/$DATASET_ID/" \
  --recursive \
  --summarize \
  --region $AWS_REGION \
  2>/dev/null | grep "Total Size" | awk '{print $3}')

if [ -n "$DATASET_SIZE" ]; then
  SIZE_MB=$((DATASET_SIZE / 1024 / 1024))
  SIZE_GB=$((SIZE_MB / 1024))

  if [ $SIZE_GB -gt 0 ]; then
    info "Dataset size: ${SIZE_GB}GB"
  else
    info "Dataset size: ${SIZE_MB}MB"
  fi
fi

# Confirm sync
echo ""
warning "This will overwrite any existing staging dataset with the same ID"
read -p "Continue? (yes/no): " confirm

if [ "$confirm" != "yes" ]; then
  echo ""
  info "Sync cancelled"
  exit 0
fi

# Perform sync
echo ""
info "Syncing files..."

aws s3 sync \
  "s3://$PROD_BUCKET/datasets/$DATASET_ID" \
  "s3://$STAGING_BUCKET/datasets/$DATASET_ID" \
  --region $AWS_REGION \
  --no-progress

echo ""
success "S3 sync complete!"

# Reminder to update database
echo ""
warning "IMPORTANT: Also update staging database with dataset metadata"
echo ""
echo "Steps:"
echo ""
echo "  1. Export dataset record from production database:"
echo ""
echo "     psql \$PRODUCTION_DATABASE_URL -c \\"
echo "       \"SELECT * FROM datasets WHERE id = '$DATASET_ID';\""
echo ""
echo "  2. Copy the record to staging database:"
echo ""
echo "     # Adjust values as needed"
echo "     psql \$STAGING_DATABASE_URL -c \\"
echo "       \"INSERT INTO datasets (id, fingerprint, title, ...) \\"
echo "       \"VALUES ('$DATASET_ID', ..., ...);\""
echo ""
echo "  3. Update manifestUrl to point to staging bucket:"
echo ""
echo "     psql \$STAGING_DATABASE_URL -c \\"
echo "       \"UPDATE datasets SET manifest_url = \\"
echo "       \"REPLACE(manifest_url, '$PROD_BUCKET', '$STAGING_BUCKET') \\"
echo "       \"WHERE id = '$DATASET_ID';\""
echo ""
echo "See docs-private/S3.md for detailed instructions."
echo ""

exit 0
