#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════
# run_dry_run.sh — the BIL pipeline on one dataset, on this machine
#
# The same stages the cluster runs, in the same order, without SLURM:
# scripts/bil-scripts/launch_combine_mmc.sh chains them with
# `sbatch --dependency=afterok`; here they are chained by `set -e`.
#
#   combine   many slices  -> one cell table + sample_offsets.json
#   mask      the artifact masks (combine --mask-only)
#   annotate  MapMyCells, when a reference is given (else --labels is reused)
#   prepare   the combined cells -> one SpatialData store
#   points    every slice's molecules -> points/molecules of that same store
#   validate  the store, by the same gate the viewer applies
#
# Usage:
#   tests/bil/run_dry_run.sh INPUT OUT NAME [--labels mapping_output.csv]
#                                           [--reference DIR --species mouse]
#                                           [--method correlation|hierarchical]
#                                           [--mask-percentile N]   (0 = none)
#                                           [--skip-combine]
#
#   INPUT   the sample folder holding the slices (read only)
#   OUT     where everything is written (never inside INPUT)
#   NAME    the store's name, OUT/store/NAME.zarr
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

if [ $# -lt 3 ]; then
    sed -n '2,25p' "$0"
    exit 1
fi

INPUT="$1"; OUT="$2"; NAME="$3"; shift 3
LABELS=""; REFERENCE=""; SPECIES="mouse"; SKIP_COMBINE=false; METHOD="correlation"; MASK_PCT=25
while [ $# -gt 0 ]; do
    case "$1" in
        --labels)    LABELS="$2"; shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --species)   SPECIES="$2"; shift 2 ;;
        --method)    METHOD="$2"; shift 2 ;;
        --mask-percentile) MASK_PCT="$2"; shift 2 ;;
        --skip-combine) SKIP_COMBINE=true; shift ;;
        *) echo "unknown option: $1"; exit 1 ;;
    esac
done

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY="${MEYES_PYTHON:-python3}"
COMBINED="${OUT}/combined_output"
STORE_DIR="${OUT}/store"
STORE="${STORE_DIR}/${NAME}.zarr"
MMC="${OUT}/mmc_output"

mkdir -p "$OUT"
echo "============================================"
echo "  BIL dry run: ${NAME}"
echo "  input:  ${INPUT}"
echo "  output: ${OUT}"
echo "  python: $($PY -V 2>&1)"
echo "============================================"

stage () { echo; echo "── $1 ──"; date '+%H:%M:%S'; }

if ! $SKIP_COMBINE; then
    stage "1/6 combine"
    $PY -I "${ROOT}/scripts/combine_slices_v3.py" "$INPUT" "$OUT"

    stage "2/6 artifact mask"
    $PY -I "${ROOT}/scripts/combine_slices_v3.py" --mask-only "$COMBINED" "$COMBINED"
else
    stage "1-2/6 combine + mask — skipped"
fi

if [ -n "$REFERENCE" ]; then
    stage "3/6 annotate (MapMyCells)"
    # correlation needs no cell_type_mapper; hierarchical does
    $PY -I "${ROOT}/scripts/map_my_cell.py" "$COMBINED" "$MMC" \
        --species "$SPECIES" --reference_dir "$REFERENCE" --method "$METHOD"
    LABELS="${MMC}/mapping_output.csv"
elif [ -n "$LABELS" ]; then
    stage "3/6 annotate — reusing ${LABELS}"
else
    stage "3/6 annotate — skipped (no --reference, no --labels)"
fi

stage "4/6 prepare (cells)"
rm -rf "$STORE_DIR"
PREPARE_ARGS=("$COMBINED" "$STORE_DIR" --name "$NAME")
[ -n "$LABELS" ] && PREPARE_ARGS+=(--labels "$LABELS")
# the artifact mask rides the same path, as an is_artifact column
MASK="${COMBINED}/artifact_mask_p${MASK_PCT}.csv"
[ "$MASK_PCT" != "0" ] && [ -f "$MASK" ] && PREPARE_ARGS+=(--labels "$MASK") && echo "  mask: p${MASK_PCT}"
$PY -I "${ROOT}/scripts/prepare.py" "${PREPARE_ARGS[@]}"

stage "5/6 points (molecules of every slice)"
$PY -I "${ROOT}/scripts/slices_to_points.py" "$COMBINED" "$STORE"

stage "6/6 validate"
$PY -I "${ROOT}/scripts/validate_store.py" "$STORE"

echo
echo "============================================"
echo "  done: ${STORE}"
du -sh "$STORE"
echo "============================================"
