#!/usr/bin/env python3
"""Download one BIL sample's legacy output from S3, without listing the bucket.

  python download_legacy_sample.py SAMPLE OUT_DIR [--workers 32] [--cells-only]

The bucket's objects are publicly readable but the bucket cannot be listed
(and the data-sync IAM user may not list it either), so `aws s3 sync` fails.
The legacy format names every file in its manifests, though, so this reads
them and fetches exactly what reconstruct_inputs.py needs, over plain HTTPS:

  meyes_output/manifest.json, coords/, obs/ (from obs/metadata.json),
                expr/index.json + every chunk it counts
  sm_output/    the single folder, or — for a multi-slice sample — one folder
                per slice named by the cells' `_sample_id`; each manifest and
                every gene file its molecule counts say exists

palettes/ and de/ are not needed and are skipped. Writes OUT_DIR/SAMPLE/.

A re-run only fetches what is missing: each file is written to `.part` and
renamed when complete, so a file that exists is whole. Then:

  python reconstruct_inputs.py OUT_DIR/SAMPLE REBUILT_DIR

Exit codes: 0 every file downloaded; 2 the sample has no meyes_output;
3 some files failed (re-run to retry them).
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from reconstruct_inputs import UNASSIGNED_SUFFIX, read_obs_column, sanitize_gene_name  # noqa: E402

BASE = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2"
SAFE = re.compile(r"^[A-Za-z0-9._-]+$")
T0 = time.perf_counter()


def log(msg: str) -> None:
    print(f"[{time.perf_counter() - T0:7.1f}s] {msg}", flush=True)


def fetch(url: str, dest: Path, retries: int = 4) -> str | None:
    """Download url to dest (atomically). None on success or when dest is
    already there; 'missing' when S3 answers 403/404; else the error."""
    if dest.exists():
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(part, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            part.rename(dest)
            return None
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):  # S3 answers 403 for a missing key when listing is denied
                part.unlink(missing_ok=True)
                return "missing"
            err = f"HTTP {e.code}"
        except Exception as e:  # network hiccup: retry
            err = str(e)
        time.sleep(2 ** attempt)
    part.unlink(missing_ok=True)
    return err


def get_all(jobs: list[tuple[str, Path]], workers: int, label: str) -> list[str]:
    """Fetch many files in parallel; returns 'path: reason' for each failure."""
    todo = [(u, p) for u, p in jobs if not p.exists()]
    log(f"{label}: {len(jobs):,} files, {len(jobs) - len(todo):,} already here")
    failed, done = [], 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch, u, p): p for u, p in todo}
        for fut in as_completed(futures):
            done += 1
            if (err := fut.result()) is not None:
                failed.append(f"{futures[fut]}: {err}")
            if done % 200 == 0 or done == len(todo):
                log(f"  {done:,}/{len(todo):,}")
    return failed


def must(url: str, dest: Path) -> bool:
    err = fetch(url, dest)
    if err:
        log(f"  {dest.name}: {err}")
    return err is None


def molecule_jobs(url: str, folder: Path) -> list[tuple[str, Path]]:
    manifest = json.loads(gzip.decompress((folder / "manifest.json.gz").read_bytes()))
    counts = manifest["genes"].get("molecule_counts", {})
    jobs = []
    for gene in manifest["genes"]["unique_gene_names"]:
        name = sanitize_gene_name(gene)
        count = counts.get(gene)
        assigned = count.get("assigned") if isinstance(count, dict) else count
        unassigned = count.get("unassigned", 0) if isinstance(count, dict) else 0
        # a file is written only when it has molecules
        if assigned is None or assigned > 0:
            jobs.append((f"{url}/genes/{name}.bin.gz", folder / "genes" / f"{name}.bin.gz"))
        if unassigned:
            jobs.append((f"{url}/genes/{name}{UNASSIGNED_SUFFIX}.bin.gz",
                         folder / "genes" / f"{name}{UNASSIGNED_SUFFIX}.bin.gz"))
    total = int(manifest["statistics"]["total_molecules"])
    log(f"  {folder.name}: {total:,} molecules, {len(jobs)} gene files")
    return jobs


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("sample", help="the sample's folder name in bil-psc-data2, e.g. ace-dip-vat")
    ap.add_argument("out_dir", type=Path, help="downloads go to OUT_DIR/SAMPLE/")
    ap.add_argument("--workers", type=int, default=32, help="parallel downloads (default 32)")
    ap.add_argument("--cells-only", action="store_true", help="skip sm_output (the molecules)")
    args = ap.parse_args(argv)
    if not SAFE.match(args.sample):
        print(f"not a sample name: {args.sample!r}", file=sys.stderr)
        return 2

    url, root = f"{BASE}/{args.sample}", args.out_dir / args.sample
    m_url, meyes = f"{url}/meyes_output", root / "meyes_output"
    log(f"{args.sample} → {root}")

    # ── cells ──
    if not must(f"{m_url}/manifest.json", meyes / "manifest.json"):
        print(f"{args.sample}: no meyes_output/manifest.json in the bucket (check the name)", file=sys.stderr)
        return 2
    for name in ("obs/metadata.json", "expr/index.json"):
        if not must(f"{m_url}/{name}", meyes / name):
            return 3
    manifest = json.loads((meyes / "manifest.json").read_text())
    obs_meta = json.loads((meyes / "obs" / "metadata.json").read_text())
    index = json.loads((meyes / "expr" / "index.json").read_text())

    jobs = [(f"{m_url}/coords/{c}.bin.gz", meyes / "coords" / f"{c}.bin.gz")
            for c in ["spatial"] + list(manifest["statistics"].get("available_embeddings", []))]
    for col, info in obs_meta.items():
        ext = ".bin.gz" if info.get("format") == "bin-v1" else ".json.gz"
        if SAFE.match(col):
            jobs.append((f"{m_url}/obs/{col}{ext}", meyes / "obs" / f"{col}{ext}"))
        else:  # a column name that is not safe as a path: download by URL-quoting it
            q = urllib.request.quote(col + ext)
            jobs.append((f"{m_url}/obs/{q}", meyes / "obs" / f"{col}{ext}".replace("/", "_")))
    jobs += [(f"{m_url}/expr/chunk_{i:05d}.bin.gz", meyes / "expr" / f"chunk_{i:05d}.bin.gz")
             for i in range(index["num_chunks"])]
    failed = get_all(jobs, args.workers, "meyes_output")

    # ── molecules ──
    if not args.cells_only:
        s_url, sm = f"{url}/sm_output", root / "sm_output"
        if fetch(f"{s_url}/manifest.json.gz", sm / "manifest.json.gz") is None:
            folders = [(s_url, sm)]
            log("sm_output: one folder (single slice)")
        elif "_sample_id" in obs_meta:
            n = manifest["statistics"]["total_cells"]
            sids, _ = read_obs_column(meyes / "obs", "_sample_id", obs_meta["_sample_id"], n)
            folders = []
            for sid in map(str, sids.categories):
                if not SAFE.match(sid):
                    log(f"  skipping slice id {sid!r}: not a safe folder name")
                    continue
                if fetch(f"{s_url}/{sid}/manifest.json.gz", sm / sid / "manifest.json.gz") is None:
                    folders.append((f"{s_url}/{sid}", sm / sid))
                else:
                    log(f"  slice {sid}: no molecules in the bucket")
            log(f"sm_output: {len(folders)} slice folders")
        else:
            folders = []
            log("sm_output: none in the bucket (cells only)")
        mjobs = [j for u, f in folders for j in molecule_jobs(u, f)]
        failed += get_all(mjobs, args.workers, "sm_output")

    size = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    log(f"{size / 1e9:.2f} GB in {root}")
    if failed:
        (root / "download_failures.txt").write_text("\n".join(failed) + "\n")
        log(f"⚠ {len(failed)} files failed — see {root}/download_failures.txt; re-run to retry")
        return 3
    print(f"\nNext:\n  python scripts/bil-scripts/reconstruct_inputs.py {root} REBUILT_DIR/{args.sample}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
