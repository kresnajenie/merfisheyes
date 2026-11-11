#!/bin/bash
#SBATCH --job-name=validate_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/validate_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The gate: nothing syncs unless this exits 0. It applies the same rules
# the viewer applies when it opens a store.
#
#   sbatch validate_store.sbatch <store.zarr>

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

STORE="$1"
MAPPING="$(dirname "$STORE")/mapping.json"

meyes_banner "validate"
echo "Store:   $STORE"

if [ ! -d "$STORE" ]; then
    echo "ERROR: $STORE does not exist"
    meyes_done 3
fi

# mapping.json names what is drawn beside the cells; prepare.py writes it
# next to the store, so check it too when it is there
ARGS=("$STORE")
if [ -f "$MAPPING" ]; then
    echo "Mapping: $MAPPING"
    ARGS+=(--mapping "$MAPPING")
fi

meyes_run validate_store.py "${ARGS[@]}"
meyes_done $?
