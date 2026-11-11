#!/bin/bash
#SBATCH --job-name=catalog_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/catalog_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The dataset's catalog entry -> catalog.json beside the store: title,
# description, species, tissue, platform, BIL link, investigator and the rest,
# from the BIL API (falling back on $MEYES_BIL_CATALOG, then the scraped CSV).
# Synced with the store; register_lamin.py and build_bil_catalog.py read it.
#
#   sbatch catalog_store.sbatch <sample_name> <store_dir>
#
# Runs after prepare (which empties <store_dir>).

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

SAMPLE_NAME="$1"; STORE_DIR="$2"

meyes_banner "catalog entry"
echo "Sample:  $SAMPLE_NAME"
echo "Writes:  ${STORE_DIR}/catalog.json"

ARGS=("$SAMPLE_NAME" "$STORE_DIR")
[ -f "$MEYES_BIL_CATALOG" ] && ARGS+=(--catalog "$MEYES_BIL_CATALOG")
meyes_run bil-scripts/catalog_entry.py "${ARGS[@]}"
meyes_done $?
