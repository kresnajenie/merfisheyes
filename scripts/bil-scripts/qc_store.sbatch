#!/bin/bash
#SBATCH --job-name=qc_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/qc_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# This run's QC numbers -> qc.json beside the store, and the dataset's row
# in the fleet accumulator (dataset_stats.csv). Sections, counts per cell,
# genes detected, how many cells the artifact mask flags, MapMyCells types
# per level. Synced with the store; register_lamin.py reads it.
#
#   sbatch qc_store.sbatch <sample_name> <store_dir> <species> <mask_percentile>
#
# Runs after prepare (which empties <store_dir>), beside points/validate.

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

SAMPLE_NAME="$1"; STORE_DIR="$2"; SPECIES="$3"; PCT="$4"

meyes_banner "qc stats"
echo "Sample:  $SAMPLE_NAME"
echo "Writes:  ${STORE_DIR}/qc.json"
echo "Fleet:   $MEYES_STATS_CSV"

[ -d "$STORE_DIR" ] || { echo "ERROR: $STORE_DIR does not exist — did prepare run?"; meyes_done 3; }

meyes_run bil-scripts/analysis-scripts/collect_dataset_stats.py \
    --datasets "$SAMPLE_NAME" \
    --meyes-base "$MEYES_BASE" \
    --species "$SPECIES" \
    --mask-percentile "$PCT" \
    --json "${STORE_DIR}/qc.json" \
    --out "$MEYES_STATS_CSV"
meyes_done $?
