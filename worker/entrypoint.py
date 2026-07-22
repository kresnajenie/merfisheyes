#!/usr/bin/env python3
"""Server-side ingestion worker entrypoint.

Downloads a dataset's raw S3 prefix, runs the chunking processor, uploads the
chunked output to datasets/{id}/, and (in Phase 3) reports progress/status back
to the app via a signed callback. In Phase 2 the callbacks are STUBBED — printed
as CALLBACK lines — and the canonical fingerprint is a placeholder (README §6).

v1 scope: single-cell (process_spatial_data.py). Single-molecule is stubbed.

Env:
  DATASET_ID          required — e.g. ds_xxx
  DATASET_KIND        single_cell (default) | single_molecule
  PROCESSING_PARAMS   optional JSON, e.g. {"kind":"single_cell","stages":{"chunk":{"chunkSize":1}}}
  AWS_S3_BUCKET       required (or S3_BUCKET)
  AWS_REGION          default us-east-1
  AWS_S3_ENDPOINT     optional — S3-compatible endpoint (e.g. MinIO), path-style
  DELETE_RAW          "true" to delete raw/{id}/ after success (default: false)
  SCRATCH_DIR         default /scratch
  CALLBACK_URL/SECRET Phase 3 — unused here beyond being echoed
"""
import hashlib
import hmac
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3


def env(name, default=None, required=False):
    val = os.environ.get(name, default)
    if required and not val:
        fail(f"Missing required env var: {name}")
    return val


def post_callback(payload, attempts=3):
    """POST the callback to the app, signed with an HMAC of the raw body.

    No-ops when CALLBACK_URL/CALLBACK_SECRET are unset (hand-testing mode), so
    the worker still runs standalone. Failures are logged, never fatal — except
    that a lost COMPLETE means the app never learns, hence the retries.
    """
    url = os.environ.get("CALLBACK_URL")
    secret = os.environ.get("CALLBACK_SECRET")

    if not url or not secret:
        return False

    body = json.dumps(payload).encode()
    sig = hmac.new(secret.encode(), body, hashlib.sha256).hexdigest()
    headers = {
        "Content-Type": "application/json",
        "x-ingest-signature": f"sha256={sig}",
    }

    # Vercel Deployment Protection blocks unauthenticated requests to protected
    # deployments (the callback would 401 with "Protected deployment"). The
    # sanctioned automation path is Protection Bypass, sent as a header.
    bypass = os.environ.get("VERCEL_BYPASS_SECRET")
    if bypass:
        headers["x-vercel-protection-bypass"] = bypass

    req = urllib.request.Request(url, data=body, method="POST", headers=headers)

    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                resp.read()
            return True
        except Exception as e:  # noqa: BLE001 — network errors are expected
            detail = ""
            if isinstance(e, urllib.error.HTTPError):
                try:
                    detail = f" body={e.read().decode()[:200]}"
                except Exception:
                    pass
            print(
                f"  !! callback POST failed (attempt {attempt}/{attempts}): {e}{detail}",
                flush=True,
            )
            if attempt < attempts:
                time.sleep(2 * attempt)

    return False


def callback(status, **extra):
    """Print the callback payload and POST it to the app when configured."""
    payload = {"status": status, **extra}
    print(f"CALLBACK {json.dumps(payload)}", flush=True)
    post_callback(payload)


_last_progress = [0.0]


def progress(stage, percent=None, min_interval=0.0):
    """Report staged progress (design §5.5). Throttled via min_interval."""
    now = time.time()
    if min_interval and (now - _last_progress[0]) < min_interval:
        return
    _last_progress[0] = now

    prog = {"stage": stage}
    if percent is not None:
        prog["percent"] = percent
    print(f"PROGRESS {json.dumps(prog)}", flush=True)
    post_callback({"status": "PROCESSING", "progress": prog})


def fail(msg):
    callback("FAILED", error=msg)
    sys.exit(1)


def content_type_for(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".gz"):
        return "application/gzip"
    if name.endswith(".json"):
        return "application/json"
    return "application/octet-stream"


def make_s3():
    endpoint = os.environ.get("AWS_S3_ENDPOINT") or os.environ.get("S3_ENDPOINT")
    kwargs = {"region_name": env("AWS_REGION", "us-east-1")}
    if endpoint:
        kwargs.update(endpoint_url=endpoint)
        return boto3.client(
            "s3",
            config=boto3.session.Config(s3={"addressing_style": "path"}),
            **kwargs,
        )
    return boto3.client("s3", **kwargs)


def download_prefix(s3, bucket, prefix, dest: Path):
    """Download every object under `prefix` into `dest`, preserving structure."""
    paginator = s3.get_paginator("list_objects_v2")
    count = 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            rel = key[len(prefix):]
            if not rel:  # the prefix "directory" placeholder, if any
                continue
            target = dest / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            s3.download_file(bucket, key, str(target))
            count += 1
            print(f"  downloaded {key} -> {target} ({obj['Size']} bytes)", flush=True)
    return count


def upload_dir(s3, bucket, src: Path, prefix, concurrency=16):
    """Upload everything under `src` to `prefix`, preserving structure.

    Uploaded concurrently: with chunkSize=1 the output is thousands of small
    per-gene objects, so the cost is per-request round-trip latency, not
    bandwidth. Overlapping the waits is the whole win. (boto3 clients are
    thread-safe; large individual files still get boto3's automatic multipart.)
    """
    files = [p for p in sorted(src.rglob("*")) if p.is_file()]
    total = len(files)
    started = time.time()
    done = 0
    errors = []

    def put(path: Path):
        rel = path.relative_to(src).as_posix()
        key = f"{prefix}{rel}"
        s3.upload_file(
            str(path), bucket, key,
            ExtraArgs={"ContentType": content_type_for(path)},
        )
        return key

    with ThreadPoolExecutor(max_workers=concurrency) as pool:
        futures = {pool.submit(put, p): p for p in files}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:  # noqa: BLE001 — collected and re-raised below
                errors.append(f"{futures[fut]}: {type(e).__name__}: {e}")
            done += 1
            if done % 200 == 0 or done == total:
                rate = done / max(time.time() - started, 1e-6)
                print(f"  uploaded {done}/{total} ({rate:.0f} obj/s)", flush=True)

    if errors:
        raise RuntimeError(f"{len(errors)} upload(s) failed; first: {errors[0]}")

    elapsed = time.time() - started
    print(
        f"  upload finished: {total} objects in {elapsed:.1f}s "
        f"({total / max(elapsed, 1e-6):.0f} obj/s, concurrency={concurrency})",
        flush=True,
    )
    return total


def delete_prefix(s3, bucket, prefix):
    paginator = s3.get_paginator("list_objects_v2")
    to_delete = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            to_delete.append({"Key": obj["Key"]})
    for i in range(0, len(to_delete), 1000):
        s3.delete_objects(Bucket=bucket, Delete={"Objects": to_delete[i:i + 1000], "Quiet": True})
    return len(to_delete)


# process_spatial_data.py logs stages like:
#   [10:16:21] [7.3s] === STEP 4: Processing expression matrix (815 genes, ...) ===
# and per-chunk lines like "  [12/815] chunk_00011 (...)". Map both to progress.
_STEP_RE = re.compile(r"===\s*STEP\s+([\w.]+):\s*(.+?)\s*===")
_ITEM_RE = re.compile(r"\[(\d+)/(\d+)\]")


def run_processor(cmd) -> int:
    """Run the processor, echoing its output and reporting staged progress."""
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    stage = "Processing"

    for line in proc.stdout:  # type: ignore[union-attr]
        sys.stdout.write(line)
        sys.stdout.flush()

        step = _STEP_RE.search(line)
        if step:
            # Keep the stage label short for the UI.
            stage = step.group(2).split("(")[0].strip()
            progress(stage)
            continue

        item = _ITEM_RE.search(line)
        if item:
            done, total = int(item.group(1)), int(item.group(2))
            if total:
                # Throttle: with chunkSize=1 this fires hundreds of times.
                progress(stage, percent=round(done * 100 / total), min_interval=3.0)

    proc.wait()
    return proc.returncode


def sync_reference(s3, bucket, prefix, dest: Path, species: str) -> Path:
    """Download only the requested species' reference data.

    Mouse is ~1.3 GB and human ~8.2 GB, so we never pull both. Fargate tasks are
    ephemeral, so this re-downloads per job; in-region that is fast enough not to
    warrant EFS yet.
    """
    species_prefix = f"{prefix.rstrip('/')}/{species}/"
    dest.mkdir(parents=True, exist_ok=True)

    started = time.time()
    n = download_prefix(s3, bucket, species_prefix, dest / species)
    if n == 0:
        raise RuntimeError(
            f"No reference objects under s3://{bucket}/{species_prefix} — "
            "upload the MapMyCells reference data before enabling the annotate stage."
        )
    print(
        f"  reference: {n} file(s) in {time.time() - started:.1f}s",
        flush=True,
    )

    # map_my_cell.py resolves paths as <reference_dir>/<species>/<file>, so the
    # reference dir is the PARENT of the species folder.
    return dest


def run_annotate(input_path: Path, out_dir: Path, reference_dir: Path, cfg: dict) -> Path:
    """Run MapMyCells, returning the mapping CSV for the chunk stage to consume."""
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "scripts/map_my_cell.py",
        str(input_path), str(out_dir),
        "--species", cfg.get("species", "mouse"),
        "--reference_dir", str(reference_dir),
    ]
    # Flattened mapping is the default; hierarchical is far slower and opt-in.
    if cfg.get("hierarchical"):
        cmd.append("--hierarchical")

    print(f"== Running: {' '.join(cmd)} ==", flush=True)
    returncode = run_processor(cmd)
    if returncode != 0:
        raise RuntimeError(f"map_my_cell.py exited {returncode}")

    csv_path = out_dir / "mapping_output.csv"
    if not csv_path.exists():
        raise RuntimeError("map_my_cell.py finished but mapping_output.csv is missing")

    return csv_path


def find_single_cell_input(raw_dir: Path) -> Path:
    """A single .h5ad → that file; otherwise the folder (Xenium/MERSCOPE)."""
    h5ads = sorted(raw_dir.rglob("*.h5ad"))
    if len(h5ads) == 1:
        return h5ads[0]
    if len(h5ads) == 0:
        return raw_dir  # folder format; process_spatial_data.py auto-detects
    raise RuntimeError(
        f"Expected a single .h5ad, found {len(h5ads)}: {[str(p) for p in h5ads]}"
    )


def read_stats(out_dir: Path) -> dict:
    """num_cells / num_genes from the chunked manifest.

    process_spatial_data.py nests these under `statistics` as
    total_cells / total_genes; the fallbacks cover other spellings.
    """
    manifest = out_dir / "manifest.json"
    if not manifest.exists():
        return {}
    try:
        data = json.loads(manifest.read_text())
    except Exception:
        return {}

    stats = data.get("statistics") or {}

    def pick(*keys):
        for src in (stats, data):
            for k in keys:
                if isinstance(src.get(k), int):
                    return src[k]
        return None

    return {
        "numCells": pick("total_cells", "num_cells", "numCells", "n_obs"),
        "numGenes": pick("total_genes", "num_genes", "numGenes", "n_vars"),
        "spatialDimensions": pick("spatial_dimensions"),
    }


def main():
    dataset_id = env("DATASET_ID", required=True)
    kind = env("DATASET_KIND", "single_cell")
    bucket = env("AWS_S3_BUCKET") or env("S3_BUCKET", required=True)
    delete_raw = env("DELETE_RAW", "false").lower() == "true"
    scratch = Path(env("SCRATCH_DIR", "/scratch"))
    concurrency = int(env("UPLOAD_CONCURRENCY", "16"))

    try:
        params = json.loads(env("PROCESSING_PARAMS", "") or "{}")
    except json.JSONDecodeError as e:
        fail(f"PROCESSING_PARAMS is not valid JSON: {e}")

    if kind != "single_cell":
        fail(f"DATASET_KIND '{kind}' not supported yet (single_cell only in v1)")

    callback("PROCESSING", datasetId=dataset_id)

    s3 = make_s3()
    scratch.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f"{dataset_id}_", dir=scratch))
    raw_dir = work / "raw"
    out_dir = work / "out"
    raw_dir.mkdir(parents=True, exist_ok=True)

    raw_prefix = f"raw/{dataset_id}/"
    print(f"== Downloading {raw_prefix} from {bucket} ==", flush=True)
    progress("Downloading raw data")
    n = download_prefix(s3, bucket, raw_prefix, raw_dir)
    if n == 0:
        fail(f"No raw objects under {raw_prefix}")

    input_path = find_single_cell_input(raw_dir)
    print(f"== Input: {input_path} ==", flush=True)

    stages = params.get("stages", {}) or {}

    # ── annotate (MapMyCells) — opt-in, runs BEFORE chunking ──────────────
    # Sequential on purpose: process_spatial_data.py computes DE stats for the
    # mapped columns, which it can only do if the labels are present while it
    # still has the expression matrix.
    mmc_csv = None
    annotate = stages.get("annotate") or None
    if annotate:
        species = annotate.get("species", "mouse")
        ref_prefix = env("MMC_REFERENCE_S3_PREFIX", "reference/mapmycells")
        print(f"== Annotate stage: MapMyCells ({species}) ==", flush=True)
        progress(f"Fetching {species} reference data")
        reference_dir = sync_reference(
            s3, bucket, ref_prefix, work / "reference", species
        )
        progress("Mapping cell types")
        mmc_csv = run_annotate(input_path, work / "mmc", reference_dir, annotate)
        print(f"== Mapping CSV: {mmc_csv} ==", flush=True)

    cmd = [sys.executable, "scripts/process_spatial_data.py", str(input_path), str(out_dir)]
    chunk = stages.get("chunk", {}) or {}
    if isinstance(chunk.get("chunkSize"), int):
        cmd += ["--chunk-size", str(chunk["chunkSize"])]
    if mmc_csv is not None:
        cmd += ["--mmc-csv", str(mmc_csv)]

    print(f"== Running: {' '.join(cmd)} ==", flush=True)
    returncode = run_processor(cmd)
    if returncode != 0:
        fail(f"process_spatial_data.py exited {returncode}")

    if not (out_dir / "manifest.json").exists():
        fail("processor finished but out/manifest.json is missing")

    out_prefix = f"datasets/{dataset_id}/"

    # The worker is the authoritative producer of datasets/{id}/ — clear it so a
    # re-run can't leave stale objects behind (e.g. old 200-genes-per-chunk
    # files sitting alongside new per-gene ones).
    stale = delete_prefix(s3, bucket, out_prefix)
    if stale:
        print(f"== Cleared {stale} stale object(s) under {out_prefix} ==", flush=True)

    print(f"== Uploading out/ to {out_prefix} ==", flush=True)
    progress("Uploading results")
    uploaded = upload_dir(s3, bucket, out_dir, out_prefix, concurrency=concurrency)

    stats = read_stats(out_dir)

    # Phase 3 computes the canonical fingerprint here (README §6) and dedups.
    fingerprint = f"STUB_{dataset_id}"
    print(f"== Fingerprint (STUB, Phase 3 computes canonical): {fingerprint} ==", flush=True)

    if delete_raw:
        removed = delete_prefix(s3, bucket, raw_prefix)
        print(f"== Deleted {removed} raw objects under {raw_prefix} ==", flush=True)
    else:
        print("== DELETE_RAW not set; leaving raw/ in place ==", flush=True)

    # stats go under `stats` — that is the shape /api/ingest/[id]/callback reads.
    callback(
        "COMPLETE",
        datasetId=dataset_id,
        manifestKey=f"{out_prefix}manifest.json",
        uploadedFiles=uploaded,
        fingerprint=fingerprint,
        stats=stats,
    )


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 — top-level guard reports via callback
        fail(f"{type(e).__name__}: {e}")
