#!/bin/bash
#SBATCH --job-name=slices_to_points
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/slices_to_points_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu
#
# Every slice's molecules -> points/molecules of the store the cells are
# in, each shifted onto the combined frame by the (dx, dy) the combine
# recorded in sample_offsets.json.
#
#   sbatch slices_to_points.sbatch <combined_output> <store.zarr>
#
# Measured: 433,023,021 molecules from 3 slices (40.6 GB of CSV) in
# 8 m 28 s, peak RSS 4.9 GB. Memory does not grow with the dataset —
# nothing is concatenated, each slice streams in batches — so 128 G is
# headroom, not a requirement.
#
# DISK: the writer spills to <store_dir>/.molecules.tmp before building
# the levels, roughly 16 bytes per molecule (7 GB for 433 M, ~180 GB at
# the 11.24 B scale). That is on the store's filesystem, not TMPDIR.

set -uo pipefail
source "${MEYES_BIL_SCRIPTS:-$HOME/merfisheyes/scripts/bil-scripts}/env.sh"

COMBINED="$1"; STORE="$2"

meyes_banner "points (molecules)"
echo "Combined: $COMBINED"
echo "Store:    $STORE"
echo "Free on store filesystem: $(df -h "$(dirname "$STORE")" | awk 'NR==2{print $4}')"

if [ ! -f "${COMBINED}/sample_offsets.json" ]; then
    echo "ERROR: no sample_offsets.json in ${COMBINED}."
    echo "       It is written by combine_slices_v3.py and is the only record of"
    echo "       where each slice's molecules are and how far that slice moved."
    echo "       An older combined_output predates it: re-run the combine."
    meyes_done 3
fi

meyes_run slices_to_points.py "$COMBINED" "$STORE"
meyes_done $?
