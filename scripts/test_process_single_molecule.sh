#!/usr/bin/env bash
#
# Test script for process_single_molecule.py
# Runs processing on all datasets in the test data directory.
#
# Usage:
#   ./scripts/test_process_single_molecule.sh
#   ./scripts/test_process_single_molecule.sh --workers 3
#
# Options:
#   --workers N   Process N files in parallel (default: 1, sequential)
#   --overwrite   Re-process datasets even if output already exists

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROCESSOR="$SCRIPT_DIR/process_single_molecule.py"

INPUT_BASE="/data/merfisheyes-test-data/test-data-sizes/single-molecule/merscope"
OUTPUT_BASE="/data/merfisheyes-test-data/test-data-sizes-run/single-molecule"

WORKERS=1
OVERWRITE=0

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        --workers)
            WORKERS="$2"
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            echo "Usage: $0 [--workers N] [--overwrite]" >&2
            exit 1
            ;;
    esac
done

# Collect datasets, skip already-processed unless --overwrite
DATASETS=()
SKIPPED=()
for dir in "$INPUT_BASE"/*/; do
    dataset_name="$(basename "$dir")"
    csv_file="$dir/detected_transcripts.csv"
    if [[ ! -f "$csv_file" ]]; then
        echo "WARNING: No detected_transcripts.csv in $dir, skipping"
        continue
    fi
    if [[ "$OVERWRITE" -eq 0 ]] && [[ -f "$OUTPUT_BASE/$dataset_name/manifest.json.gz" ]]; then
        SKIPPED+=("$dataset_name")
    else
        DATASETS+=("$dataset_name")
    fi
done

echo "========================================"
echo "Single Molecule Processing Test"
echo "========================================"
echo "Input:    $INPUT_BASE"
echo "Output:   $OUTPUT_BASE"
echo "To run:   ${#DATASETS[@]}"
echo "Skipped:  ${#SKIPPED[@]} (already processed)"
echo "Workers:  $WORKERS"
if [[ "$OVERWRITE" -eq 1 ]]; then
    echo "Overwrite: yes"
fi
echo "========================================"
for name in "${SKIPPED[@]}"; do
    echo "[SKIP]  $name"
done
echo ""

if [[ ${#DATASETS[@]} -eq 0 ]]; then
    echo "Nothing to process. Use --overwrite to re-run all."
    exit 0
fi

process_dataset() {
    local dataset_name="$1"
    local input_file="$INPUT_BASE/$dataset_name/detected_transcripts.csv"
    local output_dir="$OUTPUT_BASE/$dataset_name"

    echo "[START] $dataset_name"
    mkdir -p "$output_dir"

    python3 "$PROCESSOR" \
        "$input_file" \
        "$output_dir" \
        --dataset-type merscope \
        2>&1 | sed "s/^/  [$dataset_name] /"

    echo "[DONE]  $dataset_name"
    echo ""
}

export -f process_dataset
export INPUT_BASE OUTPUT_BASE PROCESSOR OVERWRITE

if [[ "$WORKERS" -le 1 ]]; then
    for dataset_name in "${DATASETS[@]}"; do
        process_dataset "$dataset_name"
    done
else
    printf '%s\n' "${DATASETS[@]}" | xargs -P "$WORKERS" -I {} bash -c 'process_dataset "$@"' _ {}
fi

echo "========================================"
echo "All datasets processed!"
echo "Output at: $OUTPUT_BASE"
echo "========================================"
