#!/bin/bash
#SBATCH --job-name=sync_mapping
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/sync_mapping_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --partition=compute

echo "=== Sync mapping.json files ==="
echo "Source: /bil/data/meyes/"
echo "Dest:   s3://merfisheyes-bil/bil-psc-data2/"
echo "Started at $(date)"

module load aws-cli

aws s3 sync /bil/data/meyes/ s3://merfisheyes-bil/bil-psc-data2/ \
    --exclude "*" \
    --include "*/mapping.json" \
    --include "mapping.json" \
    --size-only

echo "Done at $(date)"
