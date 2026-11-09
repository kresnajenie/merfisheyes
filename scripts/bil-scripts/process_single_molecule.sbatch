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
    echo "Usage: sbatch process_single_molecule.sbatch INPUT_PATH OUTPUT_PATH S3_PREFIX [SAMPLE_WORKERS] [WORKERS]"
    echo ""
    echo "Arguments:"
    echo "  INPUT_PATH       Path to input (directory with sample sub-folders, or single .csv/.parquet)"
    echo "  OUTPUT_PATH      Output directory for processed files"
    echo "  S3_PREFIX        S3 base URL prefix (required for directory mode)"
    echo "  SAMPLE_WORKERS   (optional) Number of samples to process concurrently (default: 4)"
    echo "  WORKERS          (optional) Number of gene file write workers per sample (default: 8)"
    echo ""
    echo "Examples:"
    echo "  # Directory mode: 10 samples, each with detected_transcripts.csv"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/dataset/ /bil/data/meyes/dataset/sm_output https://bucket.s3.us-west-2.amazonaws.com/prefix"
    echo ""
    echo "  # Custom parallelism: 8 samples at a time, 4 write workers each"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/dataset/ /bil/data/meyes/dataset/sm_output https://bucket.s3.us-west-2.amazonaws.com/prefix 8 4"
    echo ""
    echo "  # Single file mode (no S3_PREFIX needed, pass '-' as placeholder)"
    echo "  sbatch process_single_molecule.sbatch /bil/data/meyes/detected_transcripts.csv /bil/data/meyes/sm_output -"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
INPUT_PATH="$1"
OUTPUT_PATH="$2"
S3_PREFIX="$3"
SAMPLE_WORKERS="${4:-4}"
WORKERS="${5:-8}"

echo "=== process_single_molecule ==="
echo "Input:          $INPUT_PATH"
echo "Output:         $OUTPUT_PATH"
echo "S3 prefix:      $S3_PREFIX"
echo "Sample workers: $SAMPLE_WORKERS"
echo "Write workers:  $WORKERS"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

if [ -d "$INPUT_PATH" ]; then
    # Directory mode
    python process_single_molecule.py \
        "$INPUT_PATH" \
        "$OUTPUT_PATH" \
        --s3-prefix "$S3_PREFIX" \
        --sample-workers "$SAMPLE_WORKERS" \
        --workers "$WORKERS"
else
    # Single file mode
    python process_single_molecule.py \
        "$INPUT_PATH" \
        "$OUTPUT_PATH" \
        --workers "$WORKERS"
fi

echo "Job finished at $(date)"
