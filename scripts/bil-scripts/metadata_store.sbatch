#!/bin/bash
#SBATCH --job-name=metadata_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/metadata_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The dataset's metadata into the store's root zarr.json (docs/specs/metadata.md §4):
# SpatialData-DB's metadata_general and per-technology block, meyes_metadata (ids, counts,
# sections, the count-matrix hash) and meyes_provenance — from catalog.json, qc.json,
# stats.json and the store itself.
#
#   sbatch metadata_store.sbatch <store_dir>
#
# Runs after points, qc and catalog (it reads what they wrote, and the hash must see the
# final table), before validate and sync. Not validated here: the bionty check is on the
# Lamin side (scripts/metadata/validate_metadata.py, before register_lamin.py).

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

STORE_DIR="$1"

meyes_banner "metadata"
echo "Folder:  $STORE_DIR"

if [ ! -d "$STORE_DIR" ]; then
    echo "ERROR: $STORE_DIR does not exist"
    meyes_done 3
fi

meyes_run metadata/build_metadata.py "$STORE_DIR"
meyes_done $?
