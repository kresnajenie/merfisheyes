#!/bin/bash
#SBATCH --job-name=collect_stats
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/collect_stats_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=compute

module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

python3 ~/merfisheyes/scripts/bil-scripts/analysis-scripts/collect_dataset_stats.py "$@"
