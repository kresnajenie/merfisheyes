#!/usr/bin/env python3
"""Benchmark Python preprocessing speed.

Tests the Python preprocessing scripts (process_h5ad.py / process_single_molecule.py)
that convert raw data into the chunked format ready for S3 upload.

This measures the server-side / offline processing path for datasets that are too
large to process in the browser.

Usage:
    python run_python_processing.py --data-dir benchmark_data/synthetic
    python run_python_processing.py --data-dir /path/to/real/data --types single_cell
    python run_python_processing.py --filter-csv browser_benchmarks.csv  # only test failures
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    DEFAULT_RESULTS_DIR, get_size_tier, get_git_version,
)
from run_browser_benchmarks import (
    get_file_size_mb, scan_single_cell, scan_single_molecule,
    load_existing_results,
)

SCRIPTS_DIR = Path(__file__).resolve().parent.parent  # scripts/

PROCESSING_CSV_COLUMNS = [
    "platform", "data_type", "size_tier",
    "n_cells", "n_molecules", "n_genes",
    "path", "file_size_mb", "filetype",
    "processing_time_s", "output_size_mb",
    "peak_memory_mb",
    "crashed", "error", "notes",
    "timestamp", "version",
]


def append_processing_csv(csv_path: str, row: dict):
    exists = os.path.exists(csv_path)
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PROCESSING_CSV_COLUMNS)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def measure_peak_memory(pid: int) -> float | None:
    """Get peak RSS of a process in MB."""
    try:
        import psutil
        proc = psutil.Process(pid)
        info = proc.memory_info()
        return round(info.rss / (1024 * 1024), 1)
    except Exception:
        return None


def filter_crashed(filter_csv: str) -> set[str]:
    """Return paths that crashed in browser benchmarks (candidates for Python processing)."""
    paths = set()
    if not os.path.exists(filter_csv):
        return paths
    with open(filter_csv, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("crashed", "False").lower() == "true":
                paths.add(row.get("path", ""))
    return paths


def benchmark_python_single_cell(path: str, filetype: str, n_cells: int, n_genes: int) -> dict:
    """Run process_h5ad.py or process_spatial_data.py and measure time."""
    row = {
        "platform": "Python_preprocessing",
        "data_type": "single_cell",
        "size_tier": get_size_tier(n_cells, "single_cell"),
        "n_cells": n_cells, "n_molecules": None,
        "n_genes": n_genes, "path": path,
        "file_size_mb": round(get_file_size_mb(path), 2),
        "filetype": filetype,
        "processing_time_s": None, "output_size_mb": None,
        "peak_memory_mb": None,
        "crashed": False, "error": None, "notes": "",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "version": get_git_version(),
    }

    # Choose the right script
    if filetype == "h5ad":
        script = SCRIPTS_DIR / "process_h5ad.py"
    else:
        script = SCRIPTS_DIR / "process_spatial_data.py"

    if not script.exists():
        row["error"] = f"Script not found: {script}"
        row["crashed"] = True
        return row

    # Create temp output dir
    with tempfile.TemporaryDirectory(prefix="bench_py_") as tmp_dir:
        output_dir = os.path.join(tmp_dir, "output")
        cmd = [sys.executable, str(script), path, output_dir]

        t0 = time.perf_counter()
        try:
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )

            # Monitor memory periodically
            peak_mem = 0.0
            try:
                import psutil
                ps_proc = psutil.Process(proc.pid)
                while proc.poll() is None:
                    try:
                        mem = ps_proc.memory_info().rss / (1024 * 1024)
                        peak_mem = max(peak_mem, mem)
                    except Exception:
                        pass
                    time.sleep(0.5)
            except ImportError:
                proc.wait()

            stdout, stderr = proc.communicate(timeout=3600)
            elapsed = round(time.perf_counter() - t0, 2)

            if proc.returncode != 0:
                row["crashed"] = True
                row["error"] = stderr.decode()[:500]
            else:
                row["processing_time_s"] = elapsed
                if peak_mem > 0:
                    row["peak_memory_mb"] = round(peak_mem, 1)
                # Measure output size
                if os.path.isdir(output_dir):
                    total = 0
                    for root, _, files in os.walk(output_dir):
                        for f in files:
                            total += os.path.getsize(os.path.join(root, f))
                    row["output_size_mb"] = round(total / (1024 * 1024), 2)

        except subprocess.TimeoutExpired:
            proc.kill()
            row["crashed"] = True
            row["error"] = "Timeout (1h)"
        except Exception as e:
            row["crashed"] = True
            row["error"] = str(e)[:500]

    return row


def benchmark_python_single_molecule(path: str, filetype: str, n_molecules: int, n_genes: int) -> dict:
    """Run process_single_molecule.py and measure time."""
    row = {
        "platform": "Python_preprocessing",
        "data_type": "single_molecule",
        "size_tier": get_size_tier(n_molecules, "single_molecule"),
        "n_cells": None, "n_molecules": n_molecules,
        "n_genes": n_genes, "path": path,
        "file_size_mb": round(get_file_size_mb(path), 2),
        "filetype": filetype,
        "processing_time_s": None, "output_size_mb": None,
        "peak_memory_mb": None,
        "crashed": False, "error": None, "notes": "",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "version": get_git_version(),
    }

    script = SCRIPTS_DIR / "process_single_molecule.py"
    if not script.exists():
        row["error"] = f"Script not found: {script}"
        row["crashed"] = True
        return row

    with tempfile.TemporaryDirectory(prefix="bench_py_sm_") as tmp_dir:
        output_dir = os.path.join(tmp_dir, "output")
        cmd = [sys.executable, str(script), path, output_dir]

        t0 = time.perf_counter()
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            peak_mem = 0.0
            try:
                import psutil
                ps_proc = psutil.Process(proc.pid)
                while proc.poll() is None:
                    try:
                        mem = ps_proc.memory_info().rss / (1024 * 1024)
                        peak_mem = max(peak_mem, mem)
                    except Exception:
                        pass
                    time.sleep(0.5)
            except ImportError:
                proc.wait()

            stdout, stderr = proc.communicate(timeout=3600)
            elapsed = round(time.perf_counter() - t0, 2)

            if proc.returncode != 0:
                row["crashed"] = True
                row["error"] = stderr.decode()[:500]
            else:
                row["processing_time_s"] = elapsed
                if peak_mem > 0:
                    row["peak_memory_mb"] = round(peak_mem, 1)
                if os.path.isdir(output_dir):
                    total = 0
                    for root, _, files in os.walk(output_dir):
                        for f in files:
                            total += os.path.getsize(os.path.join(root, f))
                    row["output_size_mb"] = round(total / (1024 * 1024), 2)

        except subprocess.TimeoutExpired:
            proc.kill()
            row["crashed"] = True
            row["error"] = "Timeout (1h)"
        except Exception as e:
            row["crashed"] = True
            row["error"] = str(e)[:500]

    return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Benchmark Python preprocessing speed")
    parser.add_argument("--data-dir", default=DEFAULT_SYNTH_DIR)
    parser.add_argument("--output", default=os.path.join(DEFAULT_RESULTS_DIR, "python_processing.csv"))
    parser.add_argument("--filter-csv", default=None,
                        help="Only test datasets that CRASHED in browser benchmarks (fallback path)")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--types", nargs="+", default=["single_cell", "single_molecule"],
                        choices=["single_cell", "single_molecule"])
    parser.add_argument("--sc-filetypes", nargs="+", default=["h5ad"])
    parser.add_argument("--sm-filetypes", nargs="+", default=["parquet", "csv"])
    args = parser.parse_args()

    print(f"Python preprocessing benchmark")
    print(f"  Data:   {args.data_dir}")
    print(f"  Output: {args.output}")

    crashed_paths = None
    if args.filter_csv:
        crashed_paths = filter_crashed(args.filter_csv)
        print(f"  Filtering to {len(crashed_paths)} crashed datasets from {args.filter_csv}")

    tested = load_existing_results(args.output) if args.resume else set()
    completed = 0

    if "single_cell" in args.types:
        sc_datasets = scan_single_cell(args.data_dir, args.sc_filetypes)
        print(f"\nSingle cell datasets: {len(sc_datasets)}")
        for i, ds in enumerate(sc_datasets):
            if ds["path"] in tested:
                continue
            if crashed_paths and ds["path"] not in crashed_paths:
                continue
            print(f"  [{i+1}/{len(sc_datasets)}] {ds['filetype']} — "
                  f"{ds['n_cells']:,} cells × {ds['n_genes']:,} genes")

            row = benchmark_python_single_cell(
                ds["path"], ds["filetype"], ds["n_cells"], ds["n_genes"],
            )
            append_processing_csv(args.output, row)
            completed += 1
            if row["crashed"]:
                print(f"    CRASHED: {row['error']}")
            else:
                print(f"    time={row['processing_time_s']}s  "
                      f"output={row['output_size_mb']}MB  "
                      f"mem={row['peak_memory_mb']}MB")

    if "single_molecule" in args.types:
        sm_datasets = scan_single_molecule(args.data_dir, args.sm_filetypes)
        print(f"\nSingle molecule datasets: {len(sm_datasets)}")
        for i, ds in enumerate(sm_datasets):
            if ds["path"] in tested:
                continue
            if crashed_paths and ds["path"] not in crashed_paths:
                continue
            print(f"  [{i+1}/{len(sm_datasets)}] {ds['filetype']} — "
                  f"{ds['n_molecules']:,} molecules × {ds['n_genes']:,} genes")

            row = benchmark_python_single_molecule(
                ds["path"], ds["filetype"], ds["n_molecules"], ds["n_genes"],
            )
            append_processing_csv(args.output, row)
            completed += 1
            if row["crashed"]:
                print(f"    CRASHED: {row['error']}")
            else:
                print(f"    time={row['processing_time_s']}s  "
                      f"output={row['output_size_mb']}MB  "
                      f"mem={row['peak_memory_mb']}MB")

    print(f"\nDone. {completed} tests.")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    DEFAULT_SYNTH_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "benchmark_data", "synthetic",
    )
    DEFAULT_RESULTS_DIR = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "benchmark_results",
    )
    main()
