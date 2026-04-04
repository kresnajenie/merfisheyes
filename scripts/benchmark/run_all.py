#!/usr/bin/env python3
"""Main orchestrator for the MERFISHeyes benchmark pipeline.

Phases:
  1. generate    - Generate synthetic test data (single cell + single molecule)
  2. browser     - Run browser benchmarks on synthetic data
  3. real        - Run browser benchmarks on real data (requires --real-data-dir)
  4. upload      - Test upload speed (for datasets that loaded successfully)
  5. python      - Test Python preprocessing speed (for datasets that crashed)
  6. graphs      - Generate comparison graphs from results

Usage:
    # Full pipeline (generate + browser benchmarks)
    python run_all.py

    # Specific phases
    python run_all.py --phases generate browser
    python run_all.py --phases browser --max-cells 500000 --max-molecules 50000000

    # With real data
    python run_all.py --phases real --real-data-dir /path/to/test-data

    # Upload + Python processing (after browser benchmarks)
    python run_all.py --phases upload python

    # Quick smoke test
    python run_all.py --phases generate browser --sc-cells 1000 10000 --sc-genes 100 500 \\
        --sc-filetypes h5ad --sm-molecules 10000000 --sm-genes 1000 --sm-filetypes parquet

    # Dry run (show what would be generated)
    python run_all.py --phases generate --dry-run
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

DEFAULT_SYNTH_DIR = str(REPO_ROOT / "benchmark_data" / "synthetic")
DEFAULT_RESULTS_DIR = str(REPO_ROOT / "benchmark_results")


def run_cmd(cmd: list[str], label: str) -> bool:
    """Run a subprocess and print output. Returns True on success."""
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"{'='*70}\n")
    t0 = time.time()
    result = subprocess.run(cmd, cwd=str(SCRIPT_DIR))
    elapsed = time.time() - t0
    status = "OK" if result.returncode == 0 else "FAILED"
    print(f"\n  [{status}] {label} ({elapsed:.1f}s)")
    return result.returncode == 0


def main():
    all_phases = ["generate", "browser", "real", "upload", "python", "graphs"]

    parser = argparse.ArgumentParser(
        description="MERFISHeyes benchmark pipeline orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--phases", nargs="+", default=["generate", "browser"],
                        choices=all_phases,
                        help="Pipeline phases to run (default: generate browser)")

    # Directories
    parser.add_argument("--synth-dir", default=DEFAULT_SYNTH_DIR,
                        help="Synthetic data output directory")
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR,
                        help="Benchmark results directory")
    parser.add_argument("--real-data-dir", default=None,
                        help="Real data directory (required for 'real' phase)")

    # Generation options
    parser.add_argument("--sc-cells", type=int, nargs="+", default=None,
                        help="Single cell counts (default: all 12)")
    parser.add_argument("--sc-genes", type=int, nargs="+", default=None,
                        help="Single cell gene counts (default: all 8)")
    parser.add_argument("--sc-filetypes", nargs="+", default=None,
                        help="Single cell file types (default: h5ad xenium merscope)")
    parser.add_argument("--sm-molecules", type=int, nargs="+", default=None,
                        help="Single molecule counts (default: all 12)")
    parser.add_argument("--sm-genes", type=int, nargs="+", default=None,
                        help="Single molecule gene counts (default: 1000 2000)")
    parser.add_argument("--sm-filetypes", nargs="+", default=None,
                        help="Single molecule file types (default: parquet csv)")

    # Benchmark options
    parser.add_argument("--url", default="http://localhost:3000",
                        help="Dev server URL")
    parser.add_argument("--max-cells", type=int, default=None)
    parser.add_argument("--max-molecules", type=int, default=None)
    parser.add_argument("--max-genes", type=int, default=None)
    parser.add_argument("--max-file-size-gb", type=float, default=10.0)
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--resume", action="store_true",
                        help="Skip already-completed tests across all phases")

    parser.add_argument("--dry-run", action="store_true",
                        help="Show plan without executing (generate phase)")

    args = parser.parse_args()

    phases = args.phases
    print(f"MERFISHeyes Benchmark Pipeline")
    print(f"  Phases:     {', '.join(phases)}")
    print(f"  Synth dir:  {args.synth_dir}")
    print(f"  Results:    {args.results_dir}")
    if args.real_data_dir:
        print(f"  Real data:  {args.real_data_dir}")
    print()

    os.makedirs(args.results_dir, exist_ok=True)

    results_ok = True

    # --- Phase 1: Generate synthetic data ---
    if "generate" in phases:
        # Single cell
        cmd = [sys.executable, str(SCRIPT_DIR / "generate_single_cell.py"),
               "--output-dir", os.path.join(args.synth_dir, "single_cell"),
               "--max-file-size-gb", str(args.max_file_size_gb)]
        if args.sc_cells:
            cmd += ["--cells"] + [str(x) for x in args.sc_cells]
        if args.sc_genes:
            cmd += ["--genes"] + [str(x) for x in args.sc_genes]
        if args.sc_filetypes:
            cmd += ["--filetypes"] + args.sc_filetypes
        if args.dry_run:
            cmd += ["--dry-run"]
        results_ok &= run_cmd(cmd, "Generate synthetic single cell data")

        # Single molecule
        cmd = [sys.executable, str(SCRIPT_DIR / "generate_single_molecule.py"),
               "--output-dir", os.path.join(args.synth_dir, "single_molecule"),
               "--max-file-size-gb", str(args.max_file_size_gb)]
        if args.sm_molecules:
            cmd += ["--molecules"] + [str(x) for x in args.sm_molecules]
        if args.sm_genes:
            cmd += ["--genes"] + [str(x) for x in args.sm_genes]
        if args.sm_filetypes:
            cmd += ["--filetypes"] + args.sm_filetypes
        if args.dry_run:
            cmd += ["--dry-run"]
        results_ok &= run_cmd(cmd, "Generate synthetic single molecule data")

    browser_csv = os.path.join(args.results_dir, "browser_benchmarks.csv")

    # --- Phase 2: Browser benchmarks on synthetic data ---
    if "browser" in phases:
        cmd = [sys.executable, str(SCRIPT_DIR / "run_browser_benchmarks.py"),
               "--data-dir", args.synth_dir,
               "--output", browser_csv,
               "--url", args.url]
        if args.max_cells:
            cmd += ["--max-cells", str(args.max_cells)]
        if args.max_molecules:
            cmd += ["--max-molecules", str(args.max_molecules)]
        if args.max_genes:
            cmd += ["--max-genes", str(args.max_genes)]
        if args.headless:
            cmd += ["--headless"]
        if args.resume:
            cmd += ["--resume"]
        if args.sc_filetypes:
            cmd += ["--sc-filetypes"] + args.sc_filetypes
        if args.sm_filetypes:
            cmd += ["--sm-filetypes"] + args.sm_filetypes
        results_ok &= run_cmd(cmd, "Browser benchmarks on synthetic data")

    # --- Phase 3: Real data benchmarks ---
    if "real" in phases:
        if not args.real_data_dir:
            print("\nERROR: --real-data-dir required for 'real' phase")
            sys.exit(1)
        cmd = [sys.executable, str(SCRIPT_DIR / "run_real_data.py"),
               "--data-dir", args.real_data_dir,
               "--output", os.path.join(args.results_dir, "real_data_benchmarks.csv"),
               "--url", args.url]
        if args.headless:
            cmd += ["--headless"]
        if args.resume:
            cmd += ["--resume"]
        results_ok &= run_cmd(cmd, "Browser benchmarks on real data")

    # --- Phase 4: Upload speed ---
    if "upload" in phases:
        cmd = [sys.executable, str(SCRIPT_DIR / "run_upload_benchmarks.py"),
               "--data-dir", args.synth_dir,
               "--output", os.path.join(args.results_dir, "upload_benchmarks.csv"),
               "--url", args.url]
        if os.path.exists(browser_csv):
            cmd += ["--filter-csv", browser_csv]
        if args.headless:
            cmd += ["--headless"]
        if args.resume:
            cmd += ["--resume"]
        results_ok &= run_cmd(cmd, "Upload speed benchmarks")

    # --- Phase 5: Python preprocessing ---
    if "python" in phases:
        cmd = [sys.executable, str(SCRIPT_DIR / "run_python_processing.py"),
               "--data-dir", args.synth_dir,
               "--output", os.path.join(args.results_dir, "python_processing.csv")]
        if os.path.exists(browser_csv):
            cmd += ["--filter-csv", browser_csv]
        if args.resume:
            cmd += ["--resume"]
        results_ok &= run_cmd(cmd, "Python preprocessing benchmarks")

    # --- Phase 6: Generate graphs ---
    if "graphs" in phases:
        cmd = [sys.executable, str(SCRIPT_DIR / "make_graphs.py"),
               "--results-dir", args.results_dir,
               "--output-dir", str(REPO_ROOT / "benchmark_graphs")]
        results_ok &= run_cmd(cmd, "Generate comparison graphs")

    # Summary
    print(f"\n{'='*70}")
    print(f"  Pipeline complete {'(all passed)' if results_ok else '(some failures)'}")
    print(f"{'='*70}")
    print(f"\nResults directory: {args.results_dir}")
    for f in sorted(os.listdir(args.results_dir)):
        if f.endswith(".csv"):
            path = os.path.join(args.results_dir, f)
            size = os.path.getsize(path)
            # Count lines (rows)
            with open(path) as fh:
                lines = sum(1 for _ in fh) - 1  # minus header
            print(f"  {f}: {lines} results ({size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
