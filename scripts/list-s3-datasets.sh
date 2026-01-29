#!/bin/bash
#
# List All Datasets in S3 Bucket
#
# This script lists all datasets in staging or production S3 bucket
# with their sizes and file counts.
#
# Usage:
#   ./scripts/list-s3-datasets.sh [staging|production]
#
# Example:
#   ./scripts/list-s3-datasets.sh staging
#   ./scripts/list-s3-datasets.sh production
#

# Configuration
STAGING_BUCKET="merfisheyes-staging"
PROD_BUCKET="merfisheyes-production"
AWS_REGION="us-east-1"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
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

# Get environment argument
ENVIRONMENT=${1:-staging}

# Select bucket based on environment
if [ "$ENVIRONMENT" = "production" ] || [ "$ENVIRONMENT" = "prod" ]; then
  BUCKET=$PROD_BUCKET
  ENV_NAME="Production"
elif [ "$ENVIRONMENT" = "staging" ] || [ "$ENVIRONMENT" = "stage" ] || [ "$ENVIRONMENT" = "dev" ]; then
  BUCKET=$STAGING_BUCKET
  ENV_NAME="Staging"
else
  error "Invalid environment: $ENVIRONMENT"
  echo ""
  echo "Usage: $0 [staging|production]"
  echo ""
  echo "Examples:"
  echo "  $0 staging"
  echo "  $0 production"
  echo ""
  exit 1
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

# Display header
echo ""
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${CYAN}  Datasets in $ENV_NAME${NC}"
echo -e "${CYAN}  Bucket: $BUCKET${NC}"
echo -e "${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check if bucket exists and is accessible
if ! aws s3 ls "s3://$BUCKET/datasets/" --region $AWS_REGION > /dev/null 2>&1; then
  error "Cannot access bucket: $BUCKET"
  echo ""
  echo "Possible issues:"
  echo "  - AWS credentials not configured (run: aws configure)"
  echo "  - Insufficient permissions"
  echo "  - Bucket does not exist"
  echo ""
  exit 1
fi

# List datasets
info "Fetching dataset list..."

# Get list of dataset directories
DATASETS=$(aws s3 ls "s3://$BUCKET/datasets/" --region $AWS_REGION | awk '{print $2}' | sed 's/\///')

if [ -z "$DATASETS" ]; then
  warning "No datasets found in $BUCKET"
  echo ""
  exit 0
fi

# Count datasets
DATASET_COUNT=$(echo "$DATASETS" | wc -l | xargs)
info "Found $DATASET_COUNT dataset(s)"
echo ""

# Display header for table
printf "%-30s %-15s %-10s %-30s\n" "Dataset ID" "Size" "Files" "Last Modified"
echo "────────────────────────────────────────────────────────────────────────────────────"

# Iterate through datasets
TOTAL_SIZE=0

while IFS= read -r dataset_id; do
  # Skip empty lines
  [ -z "$dataset_id" ] && continue

  # Get dataset details
  DETAILS=$(aws s3 ls "s3://$BUCKET/datasets/$dataset_id/" \
    --recursive \
    --summarize \
    --region $AWS_REGION \
    2>/dev/null)

  # Extract size and file count
  SIZE=$(echo "$DETAILS" | grep "Total Size" | awk '{print $3}')
  FILE_COUNT=$(echo "$DETAILS" | grep "Total Objects" | awk '{print $3}')

  # Convert size to human-readable format
  if [ -n "$SIZE" ]; then
    SIZE_MB=$((SIZE / 1024 / 1024))
    SIZE_GB=$((SIZE / 1024 / 1024 / 1024))

    if [ $SIZE_GB -gt 0 ]; then
      SIZE_HUMAN="${SIZE_GB}GB"
    elif [ $SIZE_MB -gt 0 ]; then
      SIZE_HUMAN="${SIZE_MB}MB"
    else
      SIZE_KB=$((SIZE / 1024))
      SIZE_HUMAN="${SIZE_KB}KB"
    fi

    TOTAL_SIZE=$((TOTAL_SIZE + SIZE))
  else
    SIZE_HUMAN="0KB"
  fi

  # Get last modified date of manifest.json.gz
  LAST_MODIFIED=$(aws s3 ls "s3://$BUCKET/datasets/$dataset_id/manifest.json.gz" \
    --region $AWS_REGION \
    2>/dev/null | awk '{print $1, $2}')

  if [ -z "$LAST_MODIFIED" ]; then
    LAST_MODIFIED="N/A"
  fi

  # Display row
  printf "%-30s %-15s %-10s %-30s\n" \
    "$dataset_id" \
    "$SIZE_HUMAN" \
    "${FILE_COUNT:-0}" \
    "$LAST_MODIFIED"

done <<< "$DATASETS"

# Display total
echo "────────────────────────────────────────────────────────────────────────────────────"

TOTAL_SIZE_MB=$((TOTAL_SIZE / 1024 / 1024))
TOTAL_SIZE_GB=$((TOTAL_SIZE_MB / 1024))

if [ $TOTAL_SIZE_GB -gt 0 ]; then
  TOTAL_SIZE_HUMAN="${TOTAL_SIZE_GB}GB"
else
  TOTAL_SIZE_HUMAN="${TOTAL_SIZE_MB}MB"
fi

printf "%-30s %-15s\n" "Total ($DATASET_COUNT datasets)" "$TOTAL_SIZE_HUMAN"

echo ""

# Display cost estimate (approximate)
if [ $TOTAL_SIZE_GB -gt 0 ]; then
  MONTHLY_COST=$(echo "scale=2; $TOTAL_SIZE_GB * 0.023" | bc)
  info "Estimated monthly storage cost: \$${MONTHLY_COST} (at \$0.023/GB)"
fi

echo ""

exit 0
