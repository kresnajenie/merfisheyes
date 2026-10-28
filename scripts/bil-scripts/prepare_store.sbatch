#!/bin/bash
#SBATCH --job-name=prepare_store
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/prepare_store_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# The combined cells -> one SpatialData store. The molecules are added
# afterwards by slices_to_points.sbatch, into this same store.
#
#   sbatch prepare_store.sbatch <combined_output> <store_dir> <name> [labels.csv ...]
#
# Every label CSV is one row per cell keyed by cell_id, merged into the
# table's obs in the order given — MapMyCells' mapping_output.csv first,
# then the combine's artifact_mask_p25.csv.
#
# Measured: 407,859 cells x 300 genes in 17 s, ~3 GB. 128 G is for the
# largest BIL items, where cell_by_gene.csv is read whole.

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

COMBINED="$1"; STORE_DIR="$2"; NAME="$3"; shift 3

meyes_banner "prepare (cells)"
echo "Combined: $COMBINED"
echo "Store:    ${STORE_DIR}/${NAME}.zarr"

if [ ! -f "${COMBINED}/cell_metadata.csv" ]; then
    echo "ERROR: no cell_metadata.csv in ${COMBINED} — did combine_slices run?"
    meyes_done 3
fi

ARGS=("$COMBINED" "$STORE_DIR" --name "$NAME")
for csv in "$@"; do
    if [ -f "$csv" ]; then
        echo "Labels:   $csv"
        ARGS+=(--labels "$csv")
    else
        echo "WARNING:  $csv not found, skipped"
    fi
done

# prepare.py refuses to write into a folder that is not empty
rm -rf "$STORE_DIR"
meyes_run prepare.py "${ARGS[@]}"
meyes_done $?
