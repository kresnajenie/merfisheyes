#!/bin/bash

# Configuration - Edit these values
REMOTE_USER="your_username"
REMOTE_HOST="your_server.com"
REMOTE_PATH="/var/www/app"
PEM_FILE="$HOME/.ssh/your-key.pem"  # Path to your .pem file
BUILD_NAME="build-$(date +%Y%m%d-%H%M%S).tar.gz"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}Starting deployment process...${NC}"

# Build the project
echo -e "${GREEN}Building project...${NC}"
npm run build

if [ $? -ne 0 ]; then
    echo -e "${RED}Build failed!${NC}"
    exit 1
fi

# Create tarball of necessary files for deployment
echo -e "${GREEN}Creating deployment archive...${NC}"
echo -e "${YELLOW}Note: .env.local will NOT be transferred. Production server should have its own .env.local${NC}"

tar -czf "$BUILD_NAME" \
    .next \
    public \
    package.json \
    package-lock.json \
    next.config.js \
    prisma

if [ $? -ne 0 ]; then
    echo -e "${RED}Failed to create archive!${NC}"
    exit 1
fi

echo -e "${GREEN}Archive created: $BUILD_NAME${NC}"

# Transfer to remote server
echo -e "${GREEN}Transferring to remote server...${NC}"
scp -i "$PEM_FILE" "$BUILD_NAME" "${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PATH}/"

if [ $? -ne 0 ]; then
    echo -e "${RED}Transfer failed!${NC}"
    exit 1
fi

# SSH into server and extract
echo -e "${GREEN}Extracting and setting up on remote server...${NC}"
ssh -i "$PEM_FILE" "${REMOTE_USER}@${REMOTE_HOST}" << EOF
    cd ${REMOTE_PATH}

    echo "Extracting build files..."
    tar -xzf ${BUILD_NAME}

    echo "Installing production dependencies..."
    npm ci --production

    echo "Regenerating Prisma Client for production..."
    npx prisma generate

    echo "Restarting application..."
    pm2 restart all || npm start

    echo "Cleaning up build archive..."
    rm ${BUILD_NAME}
EOF

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Deployment successful!${NC}"
    echo -e "${YELLOW}Remember to ensure .env.local on production has:${NC}"
    echo -e "${YELLOW}  - DATABASE_URL${NC}"
    echo -e "${YELLOW}  - AWS_S3_BUCKET (or AWS_S3_BUCKET)${NC}"
    echo -e "${YELLOW}  - AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_REGION${NC}"
    echo -e "${YELLOW}  - NEXT_PUBLIC_BASE_URL=https://demo.merfisheyes.com${NC}"
    echo -e "${YELLOW}  - NEXT_PUBLIC_APP_URL=https://demo.merfisheyes.com${NC}"
    echo -e "${YELLOW}  - SENDGRID_API_KEY, SENDGRID_FROM_EMAIL${NC}"
    # Clean up local archive
    rm "$BUILD_NAME"
else
    echo -e "${RED}Deployment failed!${NC}"
    exit 1
fi
