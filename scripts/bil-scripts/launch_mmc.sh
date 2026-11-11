#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# launch_mmc.sh
#
# Runs map_my_cell ONLY (no combine, no downstream) on datasets whose MMC input
# is already in place — for re-running just the mapping step.
#
# Per sample it auto-detects the input:
#   - ${MEYES_BASE}/${sample_name}/combined_output   (CSV/MERSCOPE), else
#   - ${input_path}/cell_by_gene.h5ad or a *.h5ad     (BIL h5ad pipeline)
#
# By default runs BOTH methods: hierarchical -> mmc_output/, correlation ->
# mmc_output_corr/. Use --method to run just one.
#
# Species is read per-row from the sample file (sample_name,input_path,species);
# --species sets the default for rows that omit it.
#
# Usage:
#   ./launch_mmc.sh [--method both|hierarchical|correlation] [--species mouse|human] \
#                   <sample_name> [input_path] [species]      # single sample
#   ./launch_mmc.sh [--method ...] [--species ...] <samples.csv>
#   ./launch_mmc.sh [--method ...] [--species ...]            # default samples.csv
#
# Examples:
#   ./launch_mmc.sh samples.csv                       # both methods, per-row species
#   ./launch_mmc.sh --method correlation samples.csv  # correlation only
#   ./launch_mmc.sh --species human ace-irk-sag       # single sample, human
# ═══════════════════════════════════════════════════════════════

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MEYES_BASE="/bil/data/meyes"

METHOD="both"
SPECIES="mouse"   # default species for rows that don't specify one

# ── Parse optional flags ─────────────────────────────────────
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --method)     METHOD="$2"; shift 2 ;;
        --method=*)   METHOD="${1#*=}"; shift ;;
        --species)    SPECIES="$2"; shift 2 ;;
        --species=*)  SPECIES="${1#*=}"; shift ;;
        -h|--help)    awk 'NR>1 && /^#/{print} NR>1 && !/^#/{exit}' "$0"; exit 0 ;;
        *)            POSITIONAL+=("$1"); shift ;;
    esac
done
set -- "${POSITIONAL[@]}"

if [[ "$METHOD" != "hierarchical" && "$METHOD" != "correlation" && "$METHOD" != "both" ]]; then
    echo "ERROR: --method must be 'both', 'hierarchical', or 'correlation' (got: $METHOD)"
    exit 1
fi
if [ "$METHOD" = "both" ]; then
    METHODS=(hierarchical correlation)
else
    METHODS=("$METHOD")
fi

# ── Resolve sample source ────────────────────────────────────
if [ $# -eq 0 ]; then
    SAMPLE_FILE="${SCRIPT_DIR}/samples.csv"
elif [ $# -eq 1 ] && [ -f "$1" ]; then
    SAMPLE_FILE="$1"
elif [ $# -ge 1 ] && [ $# -le 3 ]; then
    SAMPLE_FILE=$(mktemp)
    echo "${1},${2:-},${3:-}" > "$SAMPLE_FILE"
    trap "rm -f $SAMPLE_FILE" EXIT
else
    echo "Usage:"
    echo "  $0 [--method both|hierarchical|correlation] [--species mouse|human] <sample_name> [input_path] [species]"
    echo "  $0 [--method ...] [--species ...] <samples.csv>"
    echo "  $0 [--method ...] [--species ...]"
    exit 1
fi

if [ ! -f "$SAMPLE_FILE" ]; then
    echo "ERROR: Sample file not found: $SAMPLE_FILE"
    exit 1
fi

echo "============================================"
echo "  MapMyCells (MMC only)"
echo "============================================"
echo "Source:   $SAMPLE_FILE"
echo "Methods:  ${METHODS[*]}"
echo "Species:  per-row (default: $SPECIES)"
echo ""

count=0

while IFS=',' read -r sample_name input_path row_species; do
    sample_name="$(echo "$sample_name" | xargs)"
    input_path="$(echo "${input_path:-}" | xargs)"
    row_species="$(echo "${row_species:-}" | xargs)"
    [[ "$sample_name" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue
    eff_species="${row_species:-$SPECIES}"

    count=$((count + 1))
    output_base="${MEYES_BASE}/${sample_name}"
    combined_output="${output_base}/combined_output"

    # Pick MMC input: combined_output (CSV/MERSCOPE) or an h5ad (BIL h5ad pipeline)
    if [ -d "$combined_output" ]; then
        mmc_input="$combined_output"; input_kind="combined_output"
    elif [ -n "$input_path" ] && [ -f "${input_path}/cell_by_gene.h5ad" ]; then
        mmc_input="${input_path}/cell_by_gene.h5ad"; input_kind="h5ad"
    elif [[ "$input_path" == *.h5ad ]] && [ -f "$input_path" ]; then
        mmc_input="$input_path"; input_kind="h5ad"
    else
        mmc_input=""; input_kind="none"
    fi

    echo "── Sample ${count}: ${sample_name} (${eff_species}) ──"
    echo "  Input:  ${mmc_input:-<none>} (${input_kind})"

    if [ -z "$mmc_input" ]; then
        echo "  ⚠️  SKIP: no combined_output dir and no cell_by_gene.h5ad at input_path"
        continue
    fi

    for method in "${METHODS[@]}"; do
        if [ "$method" = "correlation" ]; then
            mmc_out="${output_base}/mmc_output_corr"; jn="mmccorr_${sample_name}"
        else
            mmc_out="${output_base}/mmc_output"; jn="mmc_${sample_name}"
        fi
        job=$(sbatch --parsable \
            --job-name="$jn" \
            "${SCRIPT_DIR}/map_my_cell.sbatch" \
            "$mmc_input" \
            "$mmc_out" \
            "$eff_species" \
            "$method")
        echo "  ${method}  -> Job ${job} (${mmc_out##*/})"
    done

done < "$SAMPLE_FILE"

echo "============================================"
echo "  Submitted ${count} sample(s) × ${#METHODS[@]} method(s)"
echo "  Monitor with: squeue -u \$USER"
echo "============================================"
