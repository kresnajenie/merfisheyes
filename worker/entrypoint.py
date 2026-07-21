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
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3


def env(name, default=None, required=False):
    val = os.environ.get(name, default)
    if required and not val:
        fail(f"Missing required env var: {name}")
    return val


def callback(status, **extra):
    """Phase 2 stub: print the callback payload. Phase 3 POSTs it (signed)."""
    payload = {"status": status, **extra}
    print(f"CALLBACK {json.dumps(payload)}", flush=True)


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
    n = download_prefix(s3, bucket, raw_prefix, raw_dir)
    if n == 0:
        fail(f"No raw objects under {raw_prefix}")

    input_path = find_single_cell_input(raw_dir)
    print(f"== Input: {input_path} ==", flush=True)

    cmd = [sys.executable, "scripts/process_spatial_data.py", str(input_path), str(out_dir)]
    chunk = (params.get("stages", {}) or {}).get("chunk", {}) or {}
    if isinstance(chunk.get("chunkSize"), int):
        cmd += ["--chunk-size", str(chunk["chunkSize"])]

    print(f"== Running: {' '.join(cmd)} ==", flush=True)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        fail(f"process_spatial_data.py exited {result.returncode}")

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

    callback(
        "COMPLETE",
        datasetId=dataset_id,
        manifestKey=f"{out_prefix}manifest.json",
        uploadedFiles=uploaded,
        fingerprint=fingerprint,
        **stats,
    )


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 — top-level guard reports via callback
        fail(f"{type(e).__name__}: {e}")
