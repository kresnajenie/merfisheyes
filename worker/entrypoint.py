#!/usr/bin/env python3
"""Server-side ingestion worker entrypoint.

Downloads a dataset's raw S3 prefix, runs the chunking processor, uploads the
chunked output to datasets/{id}/, and (in Phase 3) reports progress/status back
to the app via a signed callback. In Phase 2 the callbacks are STUBBED — printed
as CALLBACK lines — and the canonical fingerprint is a placeholder (README §6).

Handles single-cell (process_spatial_data.py → manifest.json) and single-molecule
(process_single_molecule.py → manifest.json.gz + genes/*.bin.gz).

Env:
  DATASET_ID          required — e.g. ds_xxx
  DATASET_KIND        single_cell (default) | single_molecule
  PROCESSING_PARAMS   optional JSON, e.g.
                      {"kind":"single_cell","stages":{"chunk":{"chunkSize":1,
                       "mmcCsv":"__annotations__/mapping.csv"}}}
                      stages.chunk.mmcCsv names a user-supplied per-cell label
                      CSV uploaded with the dataset (see extract_annotation_csv).
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
from collections import deque
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import boto3


def env(name, default=None, required=False):
    val = os.environ.get(name, default)
    if required and not val:
        fail(f"Missing required env var: {name}")
    return val


# Consecutive callback failures. Once the endpoint is clearly unreachable
# (misconfigured URL, deployment protection, app down) we stop attempting
# best-effort progress posts entirely: without this, every progress update pays
# the full retry cost and a 2-minute job stretches past 10 minutes.
_callback_failures = [0]
_CALLBACK_GIVE_UP_AFTER = 3


def post_callback(payload, attempts=3, best_effort=False):
    """POST the callback to the app, signed with an HMAC of the raw body.

    No-ops when CALLBACK_URL/CALLBACK_SECRET are unset (hand-testing mode), so
    the worker still runs standalone. Failures are logged, never fatal — except
    that a lost COMPLETE means the app never learns, hence the retries there.

    best_effort=True (progress updates) means: one attempt, no backoff, and skip
    entirely once the endpoint has repeatedly failed. Terminal status posts stay
    persistent because losing them strands the dataset.
    """
    url = os.environ.get("CALLBACK_URL")
    secret = os.environ.get("CALLBACK_SECRET")

    if not url or not secret:
        return False

    if best_effort and _callback_failures[0] >= _CALLBACK_GIVE_UP_AFTER:
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
            _callback_failures[0] = 0

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

    _callback_failures[0] += 1
    if best_effort and _callback_failures[0] == _CALLBACK_GIVE_UP_AFTER:
        print(
            "  !! giving up on progress callbacks after "
            f"{_CALLBACK_GIVE_UP_AFTER} consecutive failures — processing continues, "
            "terminal status will still be attempted",
            flush=True,
        )

    return False


def callback(status, **extra):
    """Print the callback payload and POST it to the app when configured."""
    payload = {"status": status, **extra}
    print(f"CALLBACK {json.dumps(payload)}", flush=True)
    post_callback(payload)


_last_progress = [0.0]


def progress(stage, percent=None, min_interval=0.0):
    """Report staged progress (design §5.5). Throttled via min_interval.

    Progress is best-effort: a single attempt, no retry/backoff. Retrying here
    was actively harmful — each failed post cost ~6s of backoff, which is longer
    than the throttle window, so every update slipped through the throttle and a
    2-minute job ran for over 10 minutes when the callback URL was unreachable.
    """
    now = time.time()
    if min_interval and (now - _last_progress[0]) < min_interval:
        return
    _last_progress[0] = now

    prog = {"stage": stage}
    if percent is not None:
        prog["percent"] = percent
    print(f"PROGRESS {json.dumps(prog)}", flush=True)
    post_callback({"status": "PROCESSING", "progress": prog}, attempts=1, best_effort=True)


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


def run_processor(cmd):
    """Run the processor, echoing its output and reporting staged progress.

    Returns (returncode, reason). `reason` is the processor's own explanation of
    a failure, so the dataset's error message says something actionable — "CSV
    has 99 rows but dataset has 57,489 cells" rather than "exited 1", which is
    all the user would otherwise see in the upload overlay.
    """
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    stage = "Processing"
    # process_spatial_data.py reports failures as "❌ Error during processing: …"
    # before its traceback, so prefer that line; otherwise fall back to the tail.
    error_line = None
    recent = deque(maxlen=5)

    for line in proc.stdout:  # type: ignore[union-attr]
        sys.stdout.write(line)
        sys.stdout.flush()

        stripped = line.strip()
        if stripped:
            recent.append(stripped)
        if "Error during processing:" in line and error_line is None:
            error_line = stripped.split("Error during processing:", 1)[1].strip()

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

    reason = error_line or (" | ".join(recent) if proc.returncode else None)

    return proc.returncode, reason


def extract_annotation_csv(raw_dir: Path, dest_dir: Path, rel_key):
    """Move the user's cell-type CSV out of the raw input tree.

    The CSV is uploaded under a reserved key alongside the dataset. It must be
    moved *out* before the processor runs: for MERSCOPE the raw directory itself
    is the input, and format autodetection scans it for metadata / cell-by-gene
    CSVs — an extra CSV sitting there is asking for a misdetection.

    Identified by an explicit key rather than by sniffing for "a .csv", because
    a MERSCOPE upload is nothing but CSVs.

    Returns the moved path, or None when no CSV was supplied.
    """
    if not rel_key:
        return None

    source = raw_dir / rel_key
    if not source.is_file():
        raise RuntimeError(
            f"processingParams names an annotation CSV at '{rel_key}' but it was "
            f"not uploaded (looked in {raw_dir})"
        )

    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / source.name
    source.replace(dest)

    # An empty parent left behind still counts as a directory entry that
    # autodetection would list.
    for parent in source.parents:
        if parent == raw_dir:
            break
        try:
            parent.rmdir()
        except OSError:
            break

    print(f"== Annotation CSV: {dest} (moved out of raw input) ==", flush=True)

    return dest


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


def compute_raw_fingerprint(raw_dir: Path) -> str:
    """SHA-256 over the raw uploaded bytes.

    Streamed (1 MB chunks) so a 5 GB input never loads into memory. Files are
    hashed in relative-path order for determinism, and each file's path + size
    are folded in so a rename or truncation changes the result. Re-uploading the
    identical file(s) reproduces the same digest, which the app dedups on.
    """
    files = sorted(
        (p for p in raw_dir.rglob("*") if p.is_file()),
        key=lambda p: p.relative_to(raw_dir).as_posix(),
    )
    h = hashlib.sha256()
    for p in files:
        rel = p.relative_to(raw_dir).as_posix()
        h.update(f"{rel}:{p.stat().st_size}\n".encode())
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                h.update(chunk)

    # Bare hex digest — 64 chars, which fits the fingerprint column's
    # VarChar(64) exactly. No prefix (a "sha256:" prefix pushed it to 71 and the
    # DB rejected it). Placeholders are "raw_"/"STUB_"-prefixed, so a 64-hex
    # string is unambiguously a real content fingerprint downstream.
    return h.hexdigest()


def find_single_molecule_input(raw_dir: Path) -> Path:
    """Locate the transcripts file. The client uploads exactly one .csv/.parquet
    (the classifier selects a single detected_transcripts file), so there is
    normally one candidate; if several, prefer a detected_transcripts-style name.
    """
    candidates = sorted(
        p for p in raw_dir.rglob("*") if p.suffix.lower() in (".csv", ".parquet")
    )
    if not candidates:
        raise RuntimeError(f"No .csv/.parquet transcripts file under {raw_dir}")
    if len(candidates) == 1:
        return candidates[0]

    def score(p: Path) -> int:
        name = p.stem.lower()
        if "detected" in name and "transcript" in name:
            return 2
        if "transcript" in name:
            return 1
        return 0

    return max(candidates, key=score)


def detect_molecule_type(input_path: Path) -> str:
    """Pick the process_single_molecule.py --dataset-type from the column names.

    MERSCOPE uses global_x/global_y; Xenium uses x_location/y_location. Peeking
    the header (a single line for CSV, the schema for parquet) keeps type
    detection server-side — no schema hint has to be threaded from the client.
    """
    if input_path.suffix.lower() == ".parquet":
        import pyarrow.parquet as pq

        columns = set(pq.read_schema(str(input_path)).names)
    else:
        with open(input_path, "r", encoding="utf-8", errors="replace") as f:
            header = f.readline()
        columns = {c.strip().strip('"') for c in header.split(",")}

    if "global_x" in columns:
        return "merscope"
    if "x_location" in columns:
        return "xenium"
    raise RuntimeError(
        "Could not detect molecule schema: expected 'global_x' (MERSCOPE) or "
        f"'x_location' (Xenium) in the header. Columns: {sorted(columns)[:12]}"
    )


def read_sm_stats(out_dir: Path) -> dict:
    """Molecule/gene counts from the single-molecule manifest.

    process_single_molecule.py writes a GZIPPED manifest and nests the counts as
    statistics.total_molecules / unique_genes. numCells carries the molecule
    count — the SM API surfaces Dataset.numCells as `numMolecules`.
    """
    manifest = out_dir / "manifest.json.gz"
    if not manifest.exists():
        return {}
    try:
        import gzip

        with gzip.open(manifest, "rt", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return {}

    stats = data.get("statistics") or {}

    return {
        "numCells": stats.get("total_molecules"),
        "numGenes": stats.get("unique_genes"),
        "spatialDimensions": stats.get("spatial_dimensions"),
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

    if kind not in ("single_cell", "single_molecule"):
        fail(f"DATASET_KIND '{kind}' not supported (single_cell | single_molecule)")

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

    # Fingerprint the raw upload before anything mutates raw_dir (single-cell
    # moves the annotation CSV out). Re-uploading the identical file(s) yields
    # the same fingerprint; the app dedups on it (@unique column).
    fingerprint = compute_raw_fingerprint(raw_dir)
    print(f"== Raw fingerprint: {fingerprint} ==", flush=True)

    stages = params.get("stages", {}) or {}
    chunk = stages.get("chunk", {}) or {}

    if kind == "single_molecule":
        # SM writes a GZIPPED manifest and no genes-per-chunk concept — it's a
        # different processor and stats shape from single-cell.
        input_path = find_single_molecule_input(raw_dir)
        dataset_type = detect_molecule_type(input_path)
        print(f"== Input: {input_path} (type={dataset_type}) ==", flush=True)

        cmd = [
            sys.executable, "scripts/process_single_molecule.py",
            str(input_path), str(out_dir),
            "--dataset-type", dataset_type,
        ]
        processor_name = "process_single_molecule.py"
        manifest_name = "manifest.json.gz"
        read_output_stats = read_sm_stats
    else:
        # Before resolving the input: for a folder dataset the input IS raw_dir,
        # so the annotation CSV has to be gone before anything inspects the dir.
        mmc_csv = extract_annotation_csv(
            raw_dir, work / "annotations", chunk.get("mmcCsv")
        )
        input_path = find_single_cell_input(raw_dir)
        print(f"== Input: {input_path} ==", flush=True)

        cmd = [
            sys.executable, "scripts/process_spatial_data.py",
            str(input_path), str(out_dir),
        ]
        if isinstance(chunk.get("chunkSize"), int):
            cmd += ["--chunk-size", str(chunk["chunkSize"])]
        if mmc_csv is not None:
            cmd += ["--mmc-csv", str(mmc_csv)]
        processor_name = "process_spatial_data.py"
        manifest_name = "manifest.json"
        read_output_stats = read_stats

    print(f"== Running: {' '.join(cmd)} ==", flush=True)
    returncode, reason = run_processor(cmd)
    if returncode != 0:
        fail(reason or f"{processor_name} exited {returncode}")

    if not (out_dir / manifest_name).exists():
        fail(f"processor finished but out/{manifest_name} is missing")

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

    stats = read_output_stats(out_dir)

    if delete_raw:
        removed = delete_prefix(s3, bucket, raw_prefix)
        print(f"== Deleted {removed} raw objects under {raw_prefix} ==", flush=True)
    else:
        print("== DELETE_RAW not set; leaving raw/ in place ==", flush=True)

    # stats go under `stats` — that is the shape /api/ingest/[id]/callback reads.
    callback(
        "COMPLETE",
        datasetId=dataset_id,
        manifestKey=f"{out_prefix}{manifest_name}",
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
