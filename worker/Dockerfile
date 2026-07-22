# Server-side ingestion worker.
#
# Runs the processing stages against a dataset's raw S3 prefix and produces the
# chunked layout the viewer reads. Build from the REPO ROOT so the scripts/ copy
# resolves:
#   docker build -f worker/Dockerfile -t merfisheyes-worker .
#
# Intentionally lean: no MapMyCells reference data baked in (1.3 GB for mouse
# alone) — it is synced from S3 at job start.
#
# ARCH: local builds on Apple Silicon produce linux/arm64, which is fine for
# hand-testing. AWS Fargate runs linux/amd64 — Phase 3 / CI must build with
#   docker build --platform linux/amd64 -f worker/Dockerfile -t ... .
#
# Python 3.12 to match the validated `cellmapper` environment.
FROM python:3.12-slim

WORKDIR /app

# git is needed to install the MapMyCells packages from pinned GitHub commits.
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

# Install Python deps first for layer caching.
COPY worker/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# MapMyCells, pinned to the exact commits from the validated `cellmapper` env.
# Installed bottom-up in separate steps on purpose: resolving all three at once
# fails, because mmc_gene_mapper declares its own abc_atlas_access requirement
# that pip cannot reconcile against an explicit pin in the same pass.
# abc_atlas_access goes LAST: installing it first lets mmc_gene_mapper's own
# declared requirement overwrite it, drifting off the validated commit.
RUN pip install --no-cache-dir \
      "mmc_gene_mapper @ git+https://github.com/AllenInstitute/mmc_gene_mapper@e37a15194b1c7101d4a9df1f506eacd7dae5a4ef" \
 && pip install --no-cache-dir \
      "cell_type_mapper @ git+https://github.com/AllenInstitute/cell_type_mapper@15f9b0b1767d4937f5fd0652d303b9590eb75249" \
 && pip install --no-cache-dir \
      "abc_atlas_access @ git+https://github.com/AllenInstitute/abc_atlas_access@fbaeb667913ceb6ff72e5299c4c13907d037249e"

# Processors are reused verbatim; entrypoint orchestrates S3 + the run.
COPY scripts/process_spatial_data.py scripts/process_single_molecule.py scripts/map_my_cell.py ./scripts/
COPY worker/entrypoint.py .

ENTRYPOINT ["python", "entrypoint.py"]
