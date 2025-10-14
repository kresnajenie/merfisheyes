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
tar -czf "$BUILD_NAME" \
    .next \
    public \
    package.json \
    package-lock.json \
    next.config.js \
    .env.local

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
echo -e "${GREEN}Extracting on remote server...${NC}"
ssh -i "$PEM_FILE" "${REMOTE_USER}@${REMOTE_HOST}" << EOF
    cd ${REMOTE_PATH}
    tar -xzf ${BUILD_NAME}
    npm ci --production
    pm2 restart all || npm start
    rm ${BUILD_NAME}
EOF

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Deployment successful!${NC}"
    # Clean up local archive
    rm "$BUILD_NAME"
else
    echo -e "${RED}Deployment failed!${NC}"
    exit 1
fi
