#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_mmc.sh
#
# Runs map_my_cell on already-combined datasets under
#   ${MEYES_BASE}/${sample_name}/combined_output
#
# Usage:
#   ./launch_mmc.sh [--method hierarchical|correlation] [--species mouse|human] \
#                   <sample_name> <input_path>   # single sample
#   ./launch_mmc.sh [--method ...] [--species ...] <samples.csv>
#   ./launch_mmc.sh [--method ...] [--species ...]    # uses samples.csv in same dir
#
# Examples:
#   ./launch_mmc.sh --method correlation samples.csv
#   ./launch_mmc.sh --method correlation             # all rows in default samples.csv
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"

METHOD="hierarchical"
SPECIES="mouse"

# ── Parse optional flags ─────────────────────────────────────
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --method)
            METHOD="$2"
            shift 2
            ;;
        --method=*)
            METHOD="${1#*=}"
            shift
            ;;
        --species)
            SPECIES="$2"
            shift 2
            ;;
        --species=*)
            SPECIES="${1#*=}"
            shift
            ;;
        -h|--help)
            sed -n '2,16p' "$0"
            exit 0
            ;;
        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done
set -- "${POSITIONAL[@]}"

if [[ "$METHOD" != "hierarchical" && "$METHOD" != "correlation" ]]; then
    echo "ERROR: --method must be 'hierarchical' or 'correlation' (got: $METHOD)"
    exit 1
fi

# ── Resolve sample source ────────────────────────────────────
if [ $# -eq 0 ]; then
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
elif [ $# -eq 1 ] && [ -f "$1" ]; then
    SAMPLE_FILE="$1"
elif [ $# -eq 2 ]; then
    SAMPLE_FILE=$(mktemp)
    echo "$1,$2" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
else
    echo "Usage:"
    echo "  $0 [--method correlation|hierarchical] [--species mouse|human] <sample_name> <input_path>"
    echo "  $0 [--method ...] [--species ...] <samples.csv>"
    echo "  $0 [--method ...] [--species ...]"
    exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  MapMyCells"
echo "============================================"
echo "Source:  $SAMPLE_FILE"
echo "Method:  $METHOD"
echo "Species: $SPECIES"
echo ""

count=0

while IFS=',' read -r sample_name input_path; do
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "${input_path:-}" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"
    h5ad_input="${input_path}/cell_by_gene.h5ad"

    # Method-specific output dir so correlation never clobbers hierarchical
    if [ "$METHOD" = "correlation" ]; then
        mmc_output="${output_base}/mmc_output_corr"
    else
        mmc_output="${output_base}/mmc_output"
    fi

    # Pick MMC input: combined_output (Xenium/MERSCOPE) or h5ad (BIL h5ad pipeline)
    if [ -d "$combined_output" ]; then
        mmc_input="$combined_output"
        input_kind="combined_output"
    elif [ -n "$input_path" ] && [ -f "$h5ad_input" ]; then
        mmc_input="$h5ad_input"
        input_kind="h5ad"
    else
        mmc_input=""
        input_kind="none"
    fi

    echo "── Sample ${count}: ${sample_name} ──"
    echo "  Input:    ${mmc_input:-<none>} (${input_kind})"
    echo "  MMC:      ${mmc_output}"

    if [ -z "$mmc_input" ]; then
        echo "  ⚠️  SKIP: no combined_output dir and no cell_by_gene.h5ad at \$input_path"
        continue
    fi

    mmc_job=$(sbatch --parsable \
        --job-name="mmc_${sample_name}" \
        "${SCRIPT_DIR}/map_my_cell.sbatch" \
        "$mmc_input" \
        "$mmc_output" \
        "$SPECIES" \
        "$METHOD")
    echo "  map_my_cell  -> Job ${mmc_job}"

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
