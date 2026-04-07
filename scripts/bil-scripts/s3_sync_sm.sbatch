#!/bin/bash
#SBATCH --job-name=s3_sync_sm
#SBATCH --output=/bil/users/ijenie/meyes_process_logs/s3_sync_sm_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ijenie@ucsd.edu,eas001@ucsd.edu

# ── Help ─────────────────────────────────────────────────────
if [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ -z "$1" ]; then
    echo "Usage: sbatch s3_sync_sm.sbatch SAMPLE_NAME [--overwrite]"
    echo ""
    echo "Arguments:"
    echo "  SAMPLE_NAME   Sample name (e.g. ace-dip-use)"
    echo "  --overwrite   (optional) Re-upload everything instead of syncing"
    echo ""
    echo "Paths:"
    echo "  Source: /bil/data/meyes/{SAMPLE_NAME}/sm_output/"
    echo "  Dest:   s3://merfisheyes-bil/bil-psc-data2/{SAMPLE_NAME}/sm_output/"
    echo ""
    echo "Examples:"
    echo "  sbatch s3_sync_sm.sbatch ace-dip-use"
    echo "  sbatch s3_sync_sm.sbatch ace-dip-use --overwrite"
    exit 0
fi

# ── Arguments ────────────────────────────────────────────────
SAMPLE_NAME="$1"
SOURCE="/bil/data/meyes/${SAMPLE_NAME}/sm_output/"
DESTINATION="s3://merfisheyes-bil/bil-psc-data2/${SAMPLE_NAME}/sm_output/"
MAX_PARALLEL=8

OVERWRITE=false
if [ "${2:-}" = "--overwrite" ]; then
    OVERWRITE=true
fi

echo "=== s3_sync_sm ==="
echo "Sample:       $SAMPLE_NAME"
echo "Source:       $SOURCE"
echo "Dest:         $DESTINATION"
echo "Parallel:     $MAX_PARALLEL sample syncs"
echo "Mode:         $(if $OVERWRITE; then echo 'overwrite (cp --recursive)'; else echo 'sync (new/changed only)'; fi)"
echo "Started at $(date)"

# ── Environment setup ────────────────────────────────────────
module load aws-cli

# Increase AWS S3 transfer concurrency (default is 10)
export AWS_MAX_CONCURRENT_REQUESTS=50
aws configure set default.s3.max_concurrent_requests 50
aws configure set default.s3.multipart_chunksize 16MB

# ── Preflight check ─────────────────────────────────────────
if [ ! -d "$SOURCE" ]; then
    echo "ERROR: Source directory $SOURCE does not exist"
    exit 1
fi

# ── Sync mapping.json first (fast, single file) ─────────────
if [ -f "$SOURCE/mapping.json" ]; then
    echo ""
    echo "[$(date +%H:%M:%S)] Uploading mapping.json..."
    aws s3 cp "$SOURCE/mapping.json" "${DESTINATION}mapping.json"
fi

# ── Sync sample subdirectories in parallel ──────────────────
SAMPLE_DIRS=()
for d in "$SOURCE"/*/; do
    [ -d "$d" ] && SAMPLE_DIRS+=("$d")
done

TOTAL=${#SAMPLE_DIRS[@]}
echo ""
echo "[$(date +%H:%M:%S)] Found $TOTAL sample directories to sync"

if [ "$TOTAL" -eq 0 ]; then
    echo "No sample subdirectories found. Syncing entire directory..."
    if $OVERWRITE; then
        aws s3 cp "$SOURCE" "$DESTINATION" --recursive
    else
        aws s3 sync "$SOURCE" "$DESTINATION" --size-only
    fi
else
    # Sync samples in parallel batches
    RUNNING=0
    DONE=0
    FAILED=0
    PIDS=()
    NAMES=()

    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        sample_name=$(basename "$sample_dir")
        s3_path="${DESTINATION}${sample_name}/"

        echo "[$(date +%H:%M:%S)] Starting: $sample_name"

        if $OVERWRITE; then
            aws s3 cp "$sample_dir" "$s3_path" --recursive &
        else
            aws s3 sync "$sample_dir" "$s3_path" --size-only &
        fi

        PIDS+=($!)
        NAMES+=("$sample_name")
        RUNNING=$((RUNNING + 1))

        # Wait for a slot if we hit max parallel
        if [ "$RUNNING" -ge "$MAX_PARALLEL" ]; then
            # Wait for any one to finish
            wait -n 2>/dev/null || true
            # Reap finished processes
            NEW_PIDS=()
            NEW_NAMES=()
            for i in "${!PIDS[@]}"; do
                if kill -0 "${PIDS[$i]}" 2>/dev/null; then
                    NEW_PIDS+=("${PIDS[$i]}")
                    NEW_NAMES+=("${NAMES[$i]}")
                else
                    wait "${PIDS[$i]}" 2>/dev/null
                    STATUS=$?
                    DONE=$((DONE + 1))
                    if [ "$STATUS" -eq 0 ]; then
                        echo "[$(date +%H:%M:%S)] Done ($DONE/$TOTAL): ${NAMES[$i]}"
                    else
                        echo "[$(date +%H:%M:%S)] FAILED ($DONE/$TOTAL): ${NAMES[$i]} (exit $STATUS)"
                        FAILED=$((FAILED + 1))
                    fi
                fi
            done
            PIDS=("${NEW_PIDS[@]}")
            NAMES=("${NEW_NAMES[@]}")
            RUNNING=${#PIDS[@]}
        fi
    done

    # Wait for remaining jobs
    for i in "${!PIDS[@]}"; do
        wait "${PIDS[$i]}" 2>/dev/null
        STATUS=$?
        DONE=$((DONE + 1))
        if [ "$STATUS" -eq 0 ]; then
            echo "[$(date +%H:%M:%S)] Done ($DONE/$TOTAL): ${NAMES[$i]}"
        else
            echo "[$(date +%H:%M:%S)] FAILED ($DONE/$TOTAL): ${NAMES[$i]} (exit $STATUS)"
            FAILED=$((FAILED + 1))
        fi
    done

    echo ""
    echo "[$(date +%H:%M:%S)] Sync summary: $((TOTAL - FAILED))/$TOTAL succeeded, $FAILED failed"
fi

echo ""
echo "Sync complete at $(date)"
