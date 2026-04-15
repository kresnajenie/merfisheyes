#!/bin/bash
#SBATCH --job-name=process_h5ad_sm
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/process_h5ad_sm_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=512G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch process_h5ad_sm.sbatch CSV_PATH OUTPUT_PATH MAPPING_URL"
    echo ""
    echo "Runs process_single_molecule.py on a segmented_spot_table.csv file."
    echo ""
    echo "Arguments:"
    echo "  CSV_PATH     Path to segmented_spot_table.csv"
    echo "  OUTPUT_PATH  Output directory for processed SM files (sm_output)"
    echo "  MAPPING_URL  S3 URL for this SM dataset (used to generate mapping.json)"
    echo ""
    echo "Examples:"
    echo "  sbatch process_h5ad_sm.sbatch /bil/data/meyes/ace-low-cop/segmented_spot_table.csv /bil/data/meyes/ace-low-cop/sm_output https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2/ace-low-cop/sm_output"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
CSV_PATH="$1"
OUTPUT_PATH="$2"
MAPPING_URL="$3"

echo "========================================"
echo "  Process H5AD Single Molecule — SLURM Job $SLURM_JOB_ID"
echo "========================================"
echo "Node:       $(hostname)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Memory:     $SLURM_MEM_PER_NODE MB"
echo "CSV:        $CSV_PATH"
echo "Output:     $OUTPUT_PATH"
echo "Mapping URL: $MAPPING_URL"
echo "Started:    $(date)"
echo "========================================"
echo ""

# ── Environment setup ────────────────────────────────────────
module load gcc/11.2.0
source ~/merfisheyes_env/bin/activate

echo "Python:     $(python --version)"
echo ""

# ── Run ──────────────────────────────────────────────────────
cd ~/merfisheyes/scripts

START_TIME=$SECONDS

python process_single_molecule.py \
    "$CSV_PATH" \
    "$OUTPUT_PATH" \
    --workers 10 \
    --x-col x \
    --y-col y \
    --z-col z \
    --gene-col gene_names \
    --mapping-url "$MAPPING_URL"

EXIT_CODE=$?
ELAPSED=$(( SECONDS - START_TIME ))
MINUTES=$(( ELAPSED / 60 ))
SECS=$(( ELAPSED % 60 ))

echo ""
echo "========================================"
echo "  Job Complete"
echo "========================================"
echo "Exit code:  $EXIT_CODE"
echo "Runtime:    ${MINUTES}m ${SECS}s"
echo "Finished:   $(date)"
echo "========================================"
