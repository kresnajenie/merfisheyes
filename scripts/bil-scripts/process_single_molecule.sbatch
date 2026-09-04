#!/bin/bash
#SBATCH --job-name=process_sm
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/process_sm_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch process_single_molecule.sbatch INPUT_PATH OUTPUT_PATH [PYTHON_ARGS...]"
    echo ""
    echo "Arguments:"
    echo "  INPUT_PATH       Path to input (directory with sample sub-folders, single .csv/.parquet,"
    echo "                   or a raw Xenium run folder containing transcripts.parquet directly)"
    echo "  OUTPUT_PATH      Output directory for processed files"
    echo "  PYTHON_ARGS      Any additional flags forwarded as-is to process_single_molecule.py"
    echo "                   (e.g. --s3-prefix, --sample-workers, --workers). --s3-prefix is"
    echo "                   required by the Python script when INPUT_PATH is a directory."
    echo ""
    echo "  If INPUT_PATH/transcripts.parquet exists (raw 10x Xenium output), this script runs"
    echo "  single-file mode with --dataset-type xenium and translates any --s3-prefix in"
    echo "  PYTHON_ARGS into --mapping-url (single-file mode's equivalent flag)."
    echo ""
    echo "Examples:"
    echo "  # Directory mode: samples each with detected_transcripts.csv"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/dataset/ /bil/data/meyes/dataset/sm_output --s3-prefix https://bucket.s3.us-west-2.amazonaws.com/prefix"
    echo ""
    echo "  # Custom parallelism: 8 samples at a time, 4 write workers each"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/dataset/ /bil/data/meyes/dataset/sm_output --s3-prefix https://bucket.s3.us-west-2.amazonaws.com/prefix --sample-workers 8 --workers 4"
    echo ""
    echo "  # Single file mode (no --s3-prefix needed)"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/detected_transcripts.csv /bil/data/meyes/sm_output"
    echo ""
    echo "  # Raw Xenium run folder (auto-detected via transcripts.parquet)"
    echo "  sbatch process_single_molecule.sbatch /bil/data/xenium/Run1_Slide2_Mus_M1_F1 /bil/data/meyes/ace-pie-bit/sm_output --s3-prefix https://bucket.s3.us-west-2.amazonaws.com/ace-pie-bit/sm_output"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"
OUTPUT_PATH="$2"
shift 2

echo "=== process_single_molecule ==="
echo "Input:          $INPUT_PATH"
echo "Output:         $OUTPUT_PATH"
echo "Extra args:     $*"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

XENIUM_TRANSCRIPTS="${INPUT_PATH%/}/transcripts.parquet"

if [ -f "$XENIUM_TRANSCRIPTS" ]; then
    # Raw Xenium run folder: single-file mode. Directory-mode's --s3-prefix
    # isn't read in single-file mode, so translate it to --mapping-url.
    echo "Detected raw Xenium folder (transcripts.parquet found) -> single-file mode"
    ARGS=("$@")
    NEW_ARGS=()
    i=0
    while [ $i -lt ${#ARGS[@]} ]; do
        if [ "${ARGS[$i]}" = "--s3-prefix" ]; then
            i=$((i + 1))
            NEW_ARGS+=(--mapping-url "${ARGS[$i]}")
        else
            NEW_ARGS+=("${ARGS[$i]}")
        fi
        i=$((i + 1))
    done

    python process_single_molecule.py \
        "$XENIUM_TRANSCRIPTS" \
        "$OUTPUT_PATH" \
        --dataset-type xenium \
        "${NEW_ARGS[@]}"
else
    python process_single_molecule.py \
        "$INPUT_PATH" \
        "$OUTPUT_PATH" \
        "$@"
fi

echo "Job finished at $(date)"
