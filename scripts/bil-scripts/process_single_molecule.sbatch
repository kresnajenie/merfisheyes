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
    echo "  INPUT_PATH       Path to input (directory with sample sub-folders, or single .csv/.parquet)"
    echo "  OUTPUT_PATH      Output directory for processed files"
    echo "  PYTHON_ARGS      Any additional flags forwarded as-is to process_single_molecule.py"
    echo "                   (e.g. --s3-prefix, --sample-workers, --workers). --s3-prefix is"
    echo "                   required by the Python script when INPUT_PATH is a directory."
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

python process_single_molecule.py \
    "$INPUT_PATH" \
    "$OUTPUT_PATH" \
    "$@"

echo "Job finished at $(date)"
