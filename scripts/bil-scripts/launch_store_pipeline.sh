#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_store_pipeline.sh
#
# The whole pipeline, one SLURM chain per sample:
#
#   combine -> mask -> annotate (MapMyCells) -> prepare (cells)
#           -> points (molecules) -> validate -> sync
#
# Each stage is submitted with --dependency=afterok on the one before,
# so a failure stops that sample's chain and leaves the others alone.
#
# The output is ONE SpatialData store per dataset: the combined cells and
# every slice's molecules in one coordinate frame, which is what the new
# bucket holds and what the viewer opens by URL. It replaces the per-gene
# .bin.gz folders of the old format.
#
# Usage:
#   ./launch_store_pipeline.sh                      # samples.csv beside this script
#   ./launch_store_pipeline.sh my-samples.csv
#   ./launch_store_pipeline.sh ace-dip-tug /bil/data/18/aa/.../input
#
#   --sync              also upload to s3://$MEYES_S3_BUCKET (off by default)
#   --species human     MapMyCells taxonomy (default mouse)
#   --method correlation   MapMyCells method (default hierarchical)
#   --skip-combine      the combined_output is already there
#   --skip-mmc          reuse mmc_output/mapping_output.csv if present
#   --dry-run           print the sbatch commands, submit nothing
#
# The sample file is `sample_name,input_path` — the same shape
# launch_combine_mmc.sh takes, so samples.csv is reusable as it is.
#
# Settings (bucket, paths, interpreter) live in env.sh.
# ═══════════════════════════════════════════════════════════════
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MEYES_BIL_SCRIPTS="$SCRIPT_DIR"
source "${SCRIPT_DIR}/env.sh"

SYNC=false
SPECIES="mouse"
METHOD="hierarchical"
SKIP_COMBINE=false
SKIP_MMC=false
DRY_RUN=false
POSITIONAL=()

while [ $# -gt 0 ]; do
    case "$1" in
        --sync)         SYNC=true; shift ;;
        --species)      SPECIES="$2"; shift 2 ;;
        --method)       METHOD="$2"; shift 2 ;;
        --skip-combine) SKIP_COMBINE=true; shift ;;
        --skip-mmc)     SKIP_MMC=true; shift ;;
        --dry-run)      DRY_RUN=true; shift ;;
        -h|--help)      sed -n '2,40p' "$0"; exit 0 ;;
        --*)            echo "unknown option: $1"; exit 1 ;;
        *)              POSITIONAL+=("$1"); shift ;;
    esac
done

# ── what to run over ─────────────────────────────────────────
if [ ${#POSITIONAL[@]} -eq 0 ]; then
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
elif [ ${#POSITIONAL[@]} -eq 1 ] && [ -f "${POSITIONAL[0]}" ]; then
    SAMPLE_FILE="${POSITIONAL[0]}"
elif [ ${#POSITIONAL[@]} -eq 2 ]; then
    SAMPLE_FILE=$(mktemp)
    echo "${POSITIONAL[0]},${POSITIONAL[1]}" > "$SAMPLE_FILE"
    trap 'rm -f "$SAMPLE_FILE"' EXIT
else
    echo "Usage: $0 [samples.csv | <sample_name> <input_path>] [options]"
    exit 1
fi
[ -f "$SAMPLE_FILE" ] || { echo "ERROR: no sample file: $SAMPLE_FILE"; exit 1; }

submit () {  # submit <description> <dependency|""> <script> [args...]
    local desc="$1" dep="$2"; shift 2
    local args=(--parsable --job-name="${desc}_${sample_name}")
    [ -n "$dep" ] && args+=(--dependency=afterok:"$dep")
    if $DRY_RUN; then
        echo "       sbatch ${args[*]} ${SCRIPT_DIR}/$*" >&2
        echo "<${desc}>"
    else
        sbatch "${args[@]}" "${SCRIPT_DIR}/$@"
    fi
}

echo "============================================"
echo "  SpatialData pipeline"
echo "============================================"
echo "Samples:   $SAMPLE_FILE"
echo "Output:    ${MEYES_BASE}/<sample>"
echo "Bucket:    $(if $SYNC; then echo "s3://${MEYES_S3_BUCKET}/${MEYES_S3_PREFIX}"; else echo 'DISABLED (pass --sync)'; fi)"
echo "MapMyCells: ${METHOD}, ${SPECIES}, ${MEYES_REFERENCE}"
echo "Python:    ${MEYES_IMAGE:-$MEYES_PYTHON}"
$DRY_RUN && echo "MODE:      DRY RUN — nothing is submitted"
echo ""

count=0
while IFS=',' read -r sample_name input_path || [ -n "$sample_name" ]; do
    sample_name="$(echo "${sample_name:-}" | xargs)"
    input_path="$(echo "${input_path:-}" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue
    [[ "$sample_name" == "sample_name" ]] && continue      # a header row

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined="${output_base}/combined_output"
    mmc="${output_base}/mmc_output"
    store_dir="${output_base}/store"
    store="${store_dir}/${sample_name//-/_}.zarr"

    echo "── ${count}: ${sample_name} ──"
    echo "   input:  ${input_path:-<none>}"
    echo "   store:  ${store}"

    dep=""
    if ! $SKIP_COMBINE; then
        [ -n "$input_path" ] || { echo "   ERROR: no input_path for ${sample_name}"; continue; }
        dep=$(submit combine "" combine_slices.sbatch "$input_path" "$output_base")
        echo "   [1/7] combine         ${dep}"
        dep=$(submit mask "$dep" combine_slices.sbatch --mask-only "$combined")
        echo "   [2/7] mask            ${dep}"
    else
        echo "   [1-2/7] combine+mask  skipped"
    fi

    if $SKIP_MMC; then
        echo "   [3/7] annotate        skipped (reusing ${mmc}/mapping_output.csv)"
    else
        dep=$(submit mmc "$dep" map_my_cell.sbatch "$combined" "$mmc" "$SPECIES" "$METHOD")
        echo "   [3/7] annotate        ${dep}"
    fi

    # MapMyCells' labels and the combine's artifact mask, both one row per
    # cell keyed by cell_id, merged into the table's obs in this order
    dep=$(submit prepare "$dep" prepare_store.sbatch "$combined" "$store_dir" "${sample_name//-/_}" \
                 "${mmc}/mapping_output.csv" "${combined}/artifact_mask_p25.csv")
    echo "   [4/7] prepare (cells) ${dep}"

    dep=$(submit points "$dep" slices_to_points.sbatch "$combined" "$store")
    echo "   [5/7] points          ${dep}"

    dep=$(submit validate "$dep" validate_store.sbatch "$store")
    echo "   [6/7] validate        ${dep}"

    if $SYNC; then
        dep=$(submit sync "$dep" s3_sync_store.sbatch "$sample_name" "$store_dir")
        echo "   [7/7] sync            ${dep}"
    else
        echo "   [7/7] sync            skipped (pass --sync)"
    fi
    echo ""
done < "$SAMPLE_FILE"

echo "============================================"
echo "  ${count} sample(s)$(if $DRY_RUN; then echo ' — DRY RUN, nothing submitted'; else echo ' submitted'; fi)"
$DRY_RUN || echo "  squeue -u \$USER"
echo "============================================"
