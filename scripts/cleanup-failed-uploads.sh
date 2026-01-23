#!/bin/bash
#
# Cleanup Failed/Old Uploads from S3
#
# This script removes datasets marked as FAILED in the database from S3.
# Use this to free up storage and reduce costs.
#
# ⚠️  CAUTION: This permanently deletes files from S3
#
# Usage:
#   ./scripts/cleanup-failed-uploads.sh [staging|production] <dataset-id>
#
# Example:
#   ./scripts/cleanup-failed-uploads.sh staging cm2abc123xyz456
#   ./scripts/cleanup-failed-uploads.sh production cm2xyz789abc012
#

set -e

# Configuration
STAGING_BUCKET="merfisheyes-staging"
PROD_BUCKET="merfisheyes-production"
AWS_REGION="us-east-1"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

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
  error "Missing environment argument"
  echo ""
  echo "Usage: $0 [staging|production] <dataset-id>"
  echo ""
  echo "Examples:"
  echo "  $0 staging cm2abc123xyz456"
  echo "  $0 production cm2xyz789abc012"
  echo ""
  exit 1
fi

if [ -z "$2" ]; then
  error "Missing dataset ID"
  echo ""
  echo "Usage: $0 [staging|production] <dataset-id>"
  echo ""
  echo "To find failed datasets, query the database:"
  echo "  SELECT id, title, status, created_at"
  echo "  FROM datasets"
  echo "  WHERE status = 'FAILED'"
  echo "  ORDER BY created_at DESC;"
  echo ""
  exit 1
fi

ENVIRONMENT=$1
DATASET_ID=$2

# Select bucket
if [ "$ENVIRONMENT" = "production" ] || [ "$ENVIRONMENT" = "prod" ]; then
  BUCKET=$PROD_BUCKET
  ENV_NAME="Production"
elif [ "$ENVIRONMENT" = "staging" ] || [ "$ENVIRONMENT" = "stage" ] || [ "$ENVIRONMENT" = "dev" ]; then
  BUCKET=$STAGING_BUCKET
  ENV_NAME="Staging"
else
  error "Invalid environment: $ENVIRONMENT"
  echo ""
  echo "Valid environments: staging, production"
  echo ""
  exit 1
fi

# Display cleanup plan
echo ""
warning "Dataset Cleanup: $ENV_NAME Environment"
echo ""
echo "  Dataset ID: $DATASET_ID"
echo "  Bucket:     s3://$BUCKET/datasets/$DATASET_ID"
echo "  Action:     DELETE (permanent)"
echo ""

# Extra confirmation for production
if [ "$ENVIRONMENT" = "production" ] || [ "$ENVIRONMENT" = "prod" ]; then
  error "WARNING: You are about to delete from PRODUCTION bucket!"
  echo ""
  read -p "Type 'DELETE PRODUCTION' to confirm: " confirm

  if [ "$confirm" != "DELETE PRODUCTION" ]; then
    echo ""
    info "Cleanup cancelled"
    exit 0
  fi
else
  read -p "Type 'DELETE' to confirm: " confirm

  if [ "$confirm" != "DELETE" ]; then
    echo ""
    info "Cleanup cancelled"
    exit 0
  fi
fi

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

# Check if dataset exists
info "Checking if dataset exists in S3..."

if ! aws s3 ls "s3://$BUCKET/datasets/$DATASET_ID/" --region $AWS_REGION > /dev/null 2>&1; then
  warning "Dataset not found in S3 bucket"
  echo ""
  echo "Already deleted or never uploaded."
  exit 0
fi

success "Dataset found in S3"

# Get dataset size
info "Calculating dataset size..."
DATASET_SIZE=$(aws s3 ls "s3://$BUCKET/datasets/$DATASET_ID/" \
  --recursive \
  --summarize \
  --region $AWS_REGION \
  2>/dev/null | grep "Total Size" | awk '{print $3}')

if [ -n "$DATASET_SIZE" ]; then
  SIZE_MB=$((DATASET_SIZE / 1024 / 1024))
  SIZE_GB=$((SIZE_MB / 1024))

  if [ $SIZE_GB -gt 0 ]; then
    info "Dataset size: ${SIZE_GB}GB (will be freed)"
  else
    info "Dataset size: ${SIZE_MB}MB (will be freed)"
  fi
fi

# Delete from S3
echo ""
info "Deleting files from S3..."

aws s3 rm "s3://$BUCKET/datasets/$DATASET_ID" \
  --recursive \
  --region $AWS_REGION

echo ""
success "Dataset deleted from S3"

# Reminder to update database
echo ""
warning "IMPORTANT: Also delete dataset record from database"
echo ""
echo "Run this SQL command:"
echo ""
echo "  DELETE FROM datasets WHERE id = '$DATASET_ID';"
echo ""
echo "Or via psql:"
echo ""
if [ "$ENVIRONMENT" = "production" ] || [ "$ENVIRONMENT" = "prod" ]; then
  echo "  psql \$PRODUCTION_DATABASE_URL -c \\"
else
  echo "  psql \$STAGING_DATABASE_URL -c \\"
fi
echo "    \"DELETE FROM datasets WHERE id = '$DATASET_ID';\""
echo ""
echo "Verify deletion:"
echo "  SELECT * FROM datasets WHERE id = '$DATASET_ID';"
echo "  -- Should return 0 rows"
echo ""

exit 0
