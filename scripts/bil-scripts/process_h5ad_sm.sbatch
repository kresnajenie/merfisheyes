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

# ── Auto-detect gene/cell-id column names from the CSV header ──
# segmented_spot_table.csv (ace-lot/ace-low samples) uses gene_names/cell_ids;
# Xenium's *_selected_transcripts.csv (ace-owl samples) uses feature_name/cell_id.
HEADER=$(head -1 "$CSV_PATH")
if echo "$HEADER" | grep -qw "feature_name"; then
    GENE_COL="feature_name"
else
    GENE_COL="gene_names"
fi
if echo "$HEADER" | grep -qw "cell_id"; then
    CELL_ID_COL="cell_id"
else
    CELL_ID_COL="cell_ids"
fi
echo "Gene col:   $GENE_COL"
echo "Cell col:   $CELL_ID_COL"
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
    --gene-col "$GENE_COL" \
    --cell-id-col "$CELL_ID_COL" \
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
