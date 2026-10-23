# Server-side ingestion worker.
#
# Runs the processing stages against a dataset's raw S3 prefix and produces the
# chunked layout the viewer reads — or, for an upload made with
# NEXT_PUBLIC_INGEST_PIPELINE=prepare, runs scripts/prepare.py and produces one
# SpatialData store. Build from the REPO ROOT so the scripts/ copy resolves:
#   docker build -f worker/Dockerfile -t merfisheyes-worker .
#
# Cell-type annotation is NOT done here. Users run MapMyCells themselves and
# upload the resulting CSV with their dataset; the chunk stage merges it. That
# keeps the reference data (1.3 GB mouse / 8.2 GB human) and the Allen Institute
# packages out of the image entirely.
#
# ARCH: local builds on Apple Silicon produce linux/arm64, which is fine for
# hand-testing. AWS Fargate runs linux/amd64 — Phase 3 / CI must build with
#   docker build --platform linux/amd64 -f worker/Dockerfile -t ... .
FROM python:3.12-slim

WORKDIR /app

# Install Python deps first for layer caching.
COPY worker/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# scripts/prepare.py (the one pipeline, docs/plans/ONE-PIPELINE.md) runs in
# its own environment with exactly the packages a person installs to run it
# themselves; the legacy processors above keep theirs during the rollout.
# entrypoint.py finds it through PREPARE_PYTHON.
COPY scripts/requirements.txt /opt/prepare-requirements.txt
RUN python -m venv /opt/prepare \
    && /opt/prepare/bin/pip install --no-cache-dir -r /opt/prepare-requirements.txt
ENV PREPARE_PYTHON=/opt/prepare/bin/python
# scanpy's numba kernels cache what they compile; site-packages is not
# writable for a non-root user
ENV NUMBA_CACHE_DIR=/tmp/numba-cache
# The processors' lines reach the log (and the progress callback) as they are
# printed, not when a pipe's buffer fills
ENV PYTHONUNBUFFERED=1

# Processors are reused verbatim; entrypoint orchestrates S3 + the run.
COPY scripts/*.py ./scripts/
# The schemas scripts/validate_store.py checks every store against
COPY docs/specs/schemas ./docs/specs/schemas
COPY worker/entrypoint.py .

# prepare.py's version in each stats.json ("pipeline"); build with
#   --build-arg PIPELINE_VERSION=$(git rev-parse --short HEAD)
ARG PIPELINE_VERSION=unknown
ENV MERFISHEYES_PIPELINE_VERSION=$PIPELINE_VERSION

ENTRYPOINT ["python", "entrypoint.py"]
