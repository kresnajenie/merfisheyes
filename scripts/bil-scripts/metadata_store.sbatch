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
# final table), before validate and sync. Then the terms are checked with bionty in the
# Lamin environment ($MEYES_LAMIN_PYTHON, M3: bionty is the Lamin side's, not the
# pipeline image's) and the verdict written into meyes_provenance.validation — the same
# two steps as metadata_backfill.sh, so a store from this pipeline needs no backfill. A
# failed or impossible check is reported and does not stop the store: it syncs and
# registers with validation fail / not_validated, and stays out of the index's store rows
# until its terms are fixed.

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
rc=$?
[ $rc -eq 0 ] || meyes_done $rc

STORE="$(ls -d "${STORE_DIR%/}"/*.zarr 2>/dev/null | head -1)"
VALIDATE="${MEYES_BIL_SCRIPTS:-$MEYES_REPO/scripts/bil-scripts}/../metadata/validate_metadata.py"
if [ -z "$STORE" ]; then
    echo "WARNING: no .zarr in $STORE_DIR — terms not checked"
elif [ ! -x "$MEYES_LAMIN_PYTHON" ]; then
    echo "WARNING: no Lamin python at $MEYES_LAMIN_PYTHON — terms not checked (validation not_validated)"
else
    "$MEYES_LAMIN_PYTHON" -I "$VALIDATE" "$STORE" --write
    v=$?
    # 0 pass · 3 incomplete (organism or tissue unknown) · 4 invalid — the verdict is in the store
    [ $v -eq 0 ] || echo "WARNING: bionty check exit $v — the store syncs as validation fail/not_validated; fix its catalog.json or metadata-overrides.json and rerun this stage"
fi
meyes_done 0
