#!/usr/bin/env python3
"""Run every dataset under a folder through the processor and report what happened.

Point it at a directory and it finds the datasets inside, processes each one the
way the ingestion worker would, validates the output, and writes a results
table. Adding datasets to test means dropping them in the folder — nothing here
needs editing.

    python3 scripts/test_datasets.py ~/data/test-data-sizes/single-cell

Datasets run smallest-first, so a format bug surfaces on a 140 MB file rather
than after half an hour on a 4 GB one. Outputs are deleted after validation
(``--keep`` to retain them): a full pass over large inputs otherwise needs far
more free disk than a laptop usually has.

A dataset that fails is recorded and the run continues — one bad input should
not cost you the rest of the batch.

Discovery and format detection reuse the real functions from the processing
scripts rather than reimplementing them, so this harness cannot disagree with
the code it is testing.
"""
import argparse
import contextlib
import csv
import io
import json
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from process_spatial_data import detect_input_format  # noqa: E402
from process_single_molecule import match_transcript_file  # noqa: E402


# macOS sprays these through any folder that has been opened in Finder or copied
# from a Mac. `._x.h5ad` in particular is an AppleDouble stub, not a dataset —
# it is a few KB of resource fork that would fail with a confusing HDF5 error.
def is_junk(name: str) -> bool:
    return name.startswith(".") or name == "__MACOSX"


def dir_size(path: Path) -> int:
    if path.is_file():
        return path.stat().st_size

    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def try_detect_format(path: Path):
    """detect_input_format() for folders, or None when this isn't a dataset folder.

    It raises for anything unrecognised and prints a running commentary either
    way; during discovery both are noise.
    """
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            return detect_input_format(path)
    except Exception:
        return None


def find_transcripts(path: Path):
    """Best single-molecule transcripts file in `path`, or None.

    Scored with the processor's own matcher. That matcher only accepts `.csv`
    (its directory scan is CSV-only), but process_single_molecule.py reads
    `.parquet` perfectly well when handed the path explicitly — so a .parquet
    candidate is scored under a .csv name and then returned as itself.
    """
    best, best_score = None, 0
    for entry in sorted(path.iterdir()):
        if not entry.is_file() or is_junk(entry.name):
            continue

        name = entry.name
        if entry.suffix.lower() == ".parquet":
            name = entry.stem + ".csv"

        score, _ = match_transcript_file(name)
        if score > best_score:
            best, best_score = entry, score

    return best


def discover(root: Path):
    """Find every dataset under `root`.

    A dataset is a .h5ad file, a folder detect_input_format() recognises (Xenium
    / MERSCOPE), or a folder holding a transcripts file. A folder that is itself a
    dataset is not descended into — its contents are its parts, not more
    datasets.
    """
    found = []

    if root.is_file():
        if root.suffix == ".h5ad" and not is_junk(root.name):
            found.append({"path": root, "format": "h5ad", "kind": "single_cell"})

        return found

    if not root.is_dir():
        return found

    entries = sorted(p for p in root.iterdir() if not is_junk(p.name))

    for p in entries:
        if p.is_file() and p.suffix == ".h5ad":
            found.append({"path": p, "format": "h5ad", "kind": "single_cell"})

    fmt = try_detect_format(root)
    transcripts = find_transcripts(root)

    if fmt:
        found.append({"path": root, "format": fmt, "kind": "single_cell"})
    if transcripts:
        found.append(
            {"path": transcripts, "format": "transcripts", "kind": "single_molecule"}
        )

    if fmt or transcripts:
        return found

    for p in entries:
        if p.is_dir():
            found.extend(discover(p))

    return found


def peak_rss_mb(proc) -> float:
    """Poll the child's RSS while it runs.

    getrusage(RUSAGE_CHILDREN) would be cheaper but reports a high-water mark
    across every child the harness has ever reaped, so it cannot attribute a
    peak to one dataset.
    """
    peak = 0

    def poll():
        nonlocal peak
        while proc.poll() is None:
            try:
                out = subprocess.run(
                    ["ps", "-o", "rss=", "-p", str(proc.pid)],
                    capture_output=True, text=True, timeout=5,
                )
                rss = int(out.stdout.strip() or 0)
                peak = max(peak, rss)
            except Exception:
                pass
            time.sleep(0.5)

    t = threading.Thread(target=poll, daemon=True)
    t.start()

    return lambda: peak / 1024  # ps reports KB on macOS and Linux


def validate(out_dir: Path) -> tuple:
    """Check the output is actually loadable, not merely that the run exited 0.

    Returns (ok, problems, stats). Every path the viewer will later request is
    checked here, because a manifest that promises a file the viewer then 404s
    on is exactly the failure this harness exists to catch before an upload.
    """
    problems = []
    stats = {}

    manifest_path = out_dir / "manifest.json"
    if not manifest_path.exists():
        return False, ["manifest.json missing"], stats

    try:
        manifest = json.loads(manifest_path.read_text())
    except Exception as e:
        return False, [f"manifest.json unparseable: {e}"], stats

    m_stats = manifest.get("statistics", {})
    files = manifest.get("files", {})

    stats = {
        "cells": m_stats.get("total_cells"),
        "genes": m_stats.get("total_genes"),
        "spatial_dims": m_stats.get("spatial_dimensions"),
        "embeddings": ",".join(m_stats.get("available_embeddings") or []),
        "cluster_count": m_stats.get("cluster_count"),
        "chunks": manifest.get("processing", {}).get("num_chunks"),
    }

    # A run can produce a structurally perfect but empty output — every file in
    # place, nothing in them. That is a failure, not a pass.
    if not m_stats.get("total_cells"):
        problems.append(f"manifest reports {m_stats.get('total_cells')!r} cells")
    if not m_stats.get("total_genes"):
        problems.append(f"manifest reports {m_stats.get('total_genes')!r} genes")

    if not (out_dir / "coords" / "spatial.bin.gz").exists():
        problems.append("coords/spatial.bin.gz missing")

    for emb in m_stats.get("available_embeddings") or []:
        if not (out_dir / "coords" / f"{emb}.bin.gz").exists():
            problems.append(f"coords/{emb}.bin.gz missing (manifest lists it)")

    if not (out_dir / "expr" / "index.json").exists():
        problems.append("expr/index.json missing")

    declared_chunks = files.get("expression_chunks")
    actual_chunks = len(list((out_dir / "expr").glob("chunk_*.bin.gz")))
    if declared_chunks is not None and declared_chunks != actual_chunks:
        problems.append(
            f"manifest declares {declared_chunks} expression chunks, found {actual_chunks}"
        )
    stats["chunk_files"] = actual_chunks

    obs_meta_path = out_dir / "obs" / "metadata.json"
    obs_meta = {}
    if not obs_meta_path.exists():
        problems.append("obs/metadata.json missing")
    else:
        try:
            obs_meta = json.loads(obs_meta_path.read_text())
        except Exception as e:
            problems.append(f"obs/metadata.json unparseable: {e}")

    categorical = [c for c, m in obs_meta.items() if m.get("type") == "categorical"]
    numerical = [c for c, m in obs_meta.items() if m.get("type") == "numerical"]
    stats["obs_categorical"] = len(categorical)
    stats["obs_numerical"] = len(numerical)

    for col in files.get("observation_columns") or []:
        if col not in obs_meta:
            problems.append(f"manifest lists obs column '{col}' absent from obs/metadata.json")
        if not (out_dir / "obs" / f"{col}.json.gz").exists():
            problems.append(f"obs/{col}.json.gz missing")

    # Palettes exist for categorical columns only — a numerical column has a
    # gradient, not a discrete palette.
    for col in categorical:
        if not (out_dir / "palettes" / f"{col}.json").exists():
            problems.append(f"palettes/{col}.json missing (categorical column)")

    de_cols = files.get("de_stats") or []
    stats["de_columns"] = len(de_cols)
    for col in de_cols:
        if not (out_dir / "de" / f"{col}.bin.gz").exists():
            problems.append(f"de/{col}.bin.gz missing (manifest lists it)")

    return len(problems) == 0, problems, stats


def free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1024 ** 3


def run_one(ds, out_root: Path, args) -> dict:
    path = ds["path"]
    label = path.name
    in_size = dir_size(path)

    row = {
        "dataset": label,
        "format": ds["format"],
        "kind": ds["kind"],
        "input_mb": round(in_size / 1024 ** 2, 1),
        "status": "",
        "seconds": None,
        "peak_rss_mb": None,
        "output_mb": None,
        "objects": None,
        "error": "",
    }

    if ds["kind"] != "single_cell":
        row["status"] = "skipped"
        row["error"] = "single-molecule not wired into this harness yet"

        return row

    out_dir = out_root / f"{label}__out"
    if out_dir.exists():
        shutil.rmtree(out_dir)

    # A run that fills the disk corrupts nothing but wastes the whole batch, and
    # the output can exceed the input. Refuse rather than find out.
    needed = max(args.min_free_gb, in_size / 1024 ** 3)
    if free_gb(out_root) < needed:
        row["status"] = "skipped"
        row["error"] = f"only {free_gb(out_root):.1f} GiB free, want >{needed:.1f}"

        return row

    cmd = [
        sys.executable, str(SCRIPTS_DIR / "process_spatial_data.py"),
        str(path), str(out_dir),
        "--chunk-size", str(args.chunk_size),
    ]

    print(f"\n{'=' * 70}\n== {label}  ({row['input_mb']} MB, {ds['format']})\n{'=' * 70}", flush=True)
    print(f"   {' '.join(cmd)}", flush=True)

    started = time.time()
    log_path = out_root / f"{label}.log"

    with open(log_path, "w") as log_file:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1
        )
        read_peak = peak_rss_mb(proc)
        try:
            for line in proc.stdout:
                log_file.write(line)
                if args.verbose:
                    sys.stdout.write(line)
            proc.wait(timeout=args.timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            row["status"] = "timeout"
            row["error"] = f"exceeded {args.timeout}s"

    row["seconds"] = round(time.time() - started, 1)
    row["peak_rss_mb"] = round(read_peak(), 1)

    if row["status"] == "timeout":
        if not args.keep and out_dir.exists():
            shutil.rmtree(out_dir)

        return row

    if proc.returncode != 0:
        tail = log_path.read_text().strip().splitlines()[-3:]
        row["status"] = "failed"
        row["error"] = f"exit {proc.returncode}: {' | '.join(tail)}"[:400]
        if not args.keep and out_dir.exists():
            shutil.rmtree(out_dir)

        return row

    ok, problems, stats = validate(out_dir)
    row.update(stats)
    row["output_mb"] = round(dir_size(out_dir) / 1024 ** 2, 1)
    row["objects"] = sum(1 for p in out_dir.rglob("*") if p.is_file())
    row["status"] = "passed" if ok else "invalid"
    if not ok:
        row["error"] = "; ".join(problems)[:400]

    print(
        f"   -> {row['status']}  {row['seconds']}s  peak {row['peak_rss_mb']} MB  "
        f"{stats.get('cells')} cells x {stats.get('genes')} genes  "
        f"{row['objects']} objects / {row['output_mb']} MB",
        flush=True,
    )
    if not ok:
        for p in problems:
            print(f"      !! {p}", flush=True)

    if not args.keep:
        shutil.rmtree(out_dir)

    return row


COLUMNS = [
    "dataset", "format", "kind", "input_mb", "status", "seconds", "peak_rss_mb",
    "cells", "genes", "spatial_dims", "embeddings", "cluster_count",
    "obs_categorical", "obs_numerical", "de_columns", "chunks", "chunk_files",
    "objects", "output_mb", "error",
]


def write_results(rows, out_root: Path):
    csv_path = out_root / "results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    md_cols = [
        "dataset", "format", "input_mb", "status", "seconds", "peak_rss_mb",
        "cells", "genes", "de_columns", "objects", "output_mb",
    ]
    lines = [
        "| " + " | ".join(md_cols) + " |",
        "|" + "|".join("---" for _ in md_cols) + "|",
    ]
    for row in rows:
        # Blank only what is genuinely absent — 0 cells / 0 DE columns is a
        # finding, and `or ""` would hide it.
        cells = ["" if row.get(c) is None else str(row.get(c)) for c in md_cols]
        lines.append("| " + " | ".join(cells) + " |")

    failures = [r for r in rows if r["status"] not in ("passed", "skipped")]
    if failures:
        lines.append("")
        lines.append("### Failures")
        for row in failures:
            lines.append(f"- **{row['dataset']}** ({row['status']}): {row['error']}")

    md_path = out_root / "results.md"
    md_path.write_text("\n".join(lines) + "\n")

    return csv_path, md_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process every dataset under a folder and report the results."
    )
    parser.add_argument("root", type=Path, help="Folder (or single file) to search for datasets.")
    parser.add_argument("--out", type=Path, default=None,
                        help="Where to write outputs and results (default: ./dataset-test-results).")
    parser.add_argument("--keep", action="store_true",
                        help="Keep each processed output instead of deleting it after validation.")
    parser.add_argument("--filter", default=None,
                        help="Only run datasets whose name contains this substring.")
    parser.add_argument("--chunk-size", type=int, default=1,
                        help="Genes per chunk (default: 1, matching the ingestion worker).")
    parser.add_argument("--timeout", type=int, default=7200,
                        help="Per-dataset timeout in seconds (default: 7200).")
    parser.add_argument("--min-free-gb", type=float, default=5.0,
                        help="Skip a dataset unless at least this much disk is free (default: 5).")
    parser.add_argument("--list", action="store_true",
                        help="List the datasets that would run, then exit.")
    parser.add_argument("--verbose", action="store_true",
                        help="Echo processor output (it is always written to a per-dataset .log).")

    return parser.parse_args()


def main():
    args = parse_args()

    root = args.root.expanduser()
    if not root.exists():
        raise SystemExit(f"Not found: {root}")

    datasets = discover(root)
    if args.filter:
        datasets = [d for d in datasets if args.filter.lower() in d["path"].name.lower()]

    if not datasets:
        raise SystemExit(f"No datasets found under {root}")

    # Smallest first: a format bug should surface on the 140 MB file, not after
    # half an hour on the 4 GB one.
    datasets.sort(key=lambda d: dir_size(d["path"]))

    print(f"Found {len(datasets)} dataset(s) under {root}:\n")
    for d in datasets:
        size_mb = dir_size(d["path"]) / 1024 ** 2
        print(f"  {size_mb:9.1f} MB  {d['format']:12s} {d['kind']:16s} {d['path']}")

    if args.list:
        return

    out_root = (args.out or Path.cwd() / "dataset-test-results").expanduser()
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"\nResults -> {out_root}   ({free_gb(out_root):.1f} GiB free)\n")

    rows = []
    for ds in datasets:
        rows.append(run_one(ds, out_root, args))
        # Written every iteration so a long run can be inspected while it is
        # still going, and survives an interrupt.
        write_results(rows, out_root)

    csv_path, md_path = write_results(rows, out_root)

    counts = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1

    print(f"\n{'=' * 70}")
    print("  " + "  ".join(f"{k}: {v}" for k, v in sorted(counts.items())))
    print(f"  {csv_path}\n  {md_path}")
    print("=" * 70)

    if any(r["status"] not in ("passed", "skipped") for r in rows):
        sys.exit(1)


if __name__ == "__main__":
    main()
