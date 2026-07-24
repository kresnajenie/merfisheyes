# Server-side ingestion worker.
#
# Runs the processing stages against a dataset's raw S3 prefix and produces the
# chunked layout the viewer reads. Build from the REPO ROOT so the scripts/ copy
# resolves:
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

# Processors are reused verbatim; entrypoint orchestrates S3 + the run.
COPY scripts/process_spatial_data.py scripts/process_single_molecule.py ./scripts/
COPY worker/entrypoint.py .

ENTRYPOINT ["python", "entrypoint.py"]
