# ═══════════════════════════════════════════════════════════════
# env.sh — the one place the SpatialData pipeline's settings live
#
# Sourced by every *_store.sbatch and by launch_store_pipeline.sh, so a
# path, a bucket or the interpreter is changed here and nowhere else.
#
# Override any of these from the environment before submitting:
#   MEYES_S3_BUCKET=somewhere-else ./launch_store_pipeline.sh samples.csv --sync
#
# NOTE: #SBATCH directives are read by Slurm before this file is sourced,
# so the log directory and the mail address in the sbatch headers are
# literal. Change those in the sbatch files themselves.
# ═══════════════════════════════════════════════════════════════

# Where the checkout is, and where output goes
MEYES_REPO="${MEYES_REPO:-$HOME/merfisheyes}"
MEYES_BASE="${MEYES_BASE:-/bil/data/meyes}"
MEYES_LOGS="${MEYES_LOGS:-/bil/users/ijenie/meyes_process_logs}"

# The new bucket: one SpatialData store per dataset, read by URL
MEYES_S3_BUCKET="${MEYES_S3_BUCKET:-bil-psc-xarr}"
MEYES_S3_PREFIX="${MEYES_S3_PREFIX:-}"          # optional folder inside the bucket
MEYES_S3_REGION="${MEYES_S3_REGION:-us-west-2}"

# MapMyCells references (14 GB, not in the image — bind-mounted / read in place)
MEYES_REFERENCE="${MEYES_REFERENCE:-/bil/data/meyes/mapmycells-reference}"

# ── The interpreter ──────────────────────────────────────────
# This pipeline needs python 3.12 with scripts/requirements.txt pinned
# (spatialdata, zarr 3, anndata). That is NOT the old ~/merfisheyes_env,
# which the legacy process_*.py scripts still use — keep both.
#
#   python3.12 -m venv ~/merfisheyes_xarr_env
#   source ~/merfisheyes_xarr_env/bin/activate
#   pip install -r $MEYES_REPO/scripts/requirements.txt
#
# When the container lands (BIL-CLUSTER.md step 1) set MEYES_IMAGE and
# nothing else changes:
#   MEYES_IMAGE=/bil/data/meyes/images/merfisheyes.sif
MEYES_PYTHON="${MEYES_PYTHON:-$HOME/merfisheyes_xarr_env/bin/python}"
MEYES_IMAGE="${MEYES_IMAGE:-}"

# meyes_run <script.py> [args...] — runs one of scripts/ either way
meyes_run () {
    local script="$1"; shift
    if [ -n "$MEYES_IMAGE" ]; then
        apptainer exec \
            --bind "${MEYES_BASE}:${MEYES_BASE}" \
            --bind "${MEYES_REFERENCE}:${MEYES_REFERENCE}" \
            "$MEYES_IMAGE" \
            python -I "/opt/merfisheyes/scripts/${script}" "$@"
    else
        "$MEYES_PYTHON" -I "${MEYES_REPO}/scripts/${script}" "$@"
    fi
}

# One thread per process: these scripts parallelise by process, and MKL
# threads on top of that fight for the same cores.
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1

meyes_banner () {
    echo "========================================"
    echo "  $1 — SLURM job ${SLURM_JOB_ID:-local}"
    echo "========================================"
    echo "Node:    $(hostname)"
    echo "CPUs:    ${SLURM_CPUS_PER_TASK:-?}"
    echo "Python:  $(if [ -n "$MEYES_IMAGE" ]; then echo "apptainer $MEYES_IMAGE"; else echo "$MEYES_PYTHON"; fi)"
    echo "Started: $(date)"
    echo "========================================"
}

meyes_done () {
    local code=$1
    echo "========================================"
    echo "Exit code: ${code}   Finished: $(date)"
    echo "========================================"
    exit "$code"
}
