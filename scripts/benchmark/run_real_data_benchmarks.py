#!/usr/bin/env python3
"""
Benchmark MERFISHeyes with real test data from /data/merfisheyes-test-data/test-data-sizes/

Tests:
  - Single Molecule MERSCOPE CSVs (800MB, 2GB, 5GB, 6GB, 11GB)
  - Single Cell MERSCOPE folders (150MB, 800MB, 1.2GB)
  - Single Cell Xenium folders (10GB, 50GB, 65GB)

Usage:
    python run_real_data_benchmarks.py
    python run_real_data_benchmarks.py --only sm          # single molecule only
    python run_real_data_benchmarks.py --only sc          # single cell only
    python run_real_data_benchmarks.py --max-size 2000    # skip files > 2GB
    python run_real_data_benchmarks.py --headless         # headless browser
"""

import argparse
import asyncio
import json
import os
import shutil
import tempfile
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional


TEST_DATA_ROOT = "/data/merfisheyes-test-data/test-data-sizes"

# Define all test datasets
SINGLE_MOLECULE_DATASETS = [
    {
        "name": "SM MERSCOPE 800MB (10M molecules)",
        "path": f"{TEST_DATA_ROOT}/single-molecule/merscope/800m-ace-dry-dog/detected_transcripts.csv",
        "size_tier": "small",
        "n_molecules": 10_000_000,
        "upload_type": "merscope",
    },
    {
        "name": "SM MERSCOPE 2GB (32M molecules)",
        "path": f"{TEST_DATA_ROOT}/single-molecule/merscope/2g-ace-ear-nap-202204011403_20220401135320220401TREM5xADF12mo4hemip_VMSC00101-Region_0/detected_transcripts.csv",
        "size_tier": "medium",
        "n_molecules": 31_886_321,
        "upload_type": "merscope",
    },
    {
        "name": "SM MERSCOPE 5GB (67M molecules)",
        "path": f"{TEST_DATA_ROOT}/single-molecule/merscope/5g-ace-dry-dog/detected_transcripts.csv",
        "size_tier": "large",
        "n_molecules": 67_494_355,
        "upload_type": "merscope",
    },
    {
        "name": "SM MERSCOPE 6GB (66M molecules)",
        "path": f"{TEST_DATA_ROOT}/single-molecule/merscope/6g-ace-dip-use-1231122498/detected_transcripts.csv",
        "size_tier": "large",
        "n_molecules": 66_407_574,
        "upload_type": "merscope",
    },
    {
        "name": "SM MERSCOPE 11GB (121M molecules)",
        "path": f"{TEST_DATA_ROOT}/single-molecule/merscope/11g-ace-dip-use-1231122504/detected_transcripts.csv",
        "size_tier": "xl",
        "n_molecules": 120_707_002,
        "upload_type": "merscope",
    },
]

SINGLE_CELL_MERSCOPE_DATASETS = [
    {
        "name": "SC MERSCOPE 150MB (Ovarian Cancer)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/merscope/Human_Ovarian_Cancer_150mb",
        "size_tier": "small",
        "n_cells": None,  # will be determined from data
    },
    {
        "name": "SC MERSCOPE 800MB (Uterine Cancer)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/merscope/Human_Uterine_Cancer_800mb",
        "size_tier": "medium",
        "n_cells": None,
    },
    {
        "name": "SC MERSCOPE 1.2GB (Colon Cancer)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/merscope/Human_Colon_Cancer_1.2gb",
        "size_tier": "large",
        "n_cells": None,
    },
]

SINGLE_CELL_XENIUM_DATASETS = [
    {
        "name": "SC Xenium 10GB (Melanoma)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/xenium/Human_Melanoma_10gb",
        "size_tier": "large",
        "n_cells": None,
    },
    {
        "name": "SC Xenium 50GB (Renal Carcinoma)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/xenium/Human_Renal_Carcinoma_50gb",
        "size_tier": "xl",
        "n_cells": None,
    },
    {
        "name": "SC Xenium 65GB (Whole Mouse Pup)",
        "folder": f"{TEST_DATA_ROOT}/single-cell/xenium/Whole_mouse_pup_65gb",
        "size_tier": "xl",
        "n_cells": None,
    },
]


@dataclass
class BenchmarkResult:
    platform: str
    dataset_name: str
    data_type: str
    size_tier: str
    n_items: Optional[int]
    file_size_mb: float
    load_time_s: Optional[float] = None
    gene_query_time_s: Optional[float] = None
    memory_peak_mb: Optional[float] = None
    ui_responsive: Optional[bool] = None
    fps_after_load: Optional[float] = None
    crashed: bool = False
    error: Optional[str] = None
    notes: str = ""
    timestamp: str = field(default_factory=lambda: time.strftime("%Y-%m-%d %H:%M:%S"))


def get_file_size_mb(path: str) -> float:
    """Get file size in MB, following symlinks."""
    real_path = os.path.realpath(path)
    if os.path.isfile(real_path):
        return round(os.path.getsize(real_path) / (1024 * 1024), 1)
    elif os.path.isdir(real_path):
        total = 0
        for dirpath, _, filenames in os.walk(real_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                if os.path.isfile(fp):
                    total += os.path.getsize(fp)
        return round(total / (1024 * 1024), 1)
    return 0


def get_folder_size_mb(folder: str) -> float:
    """Get total size of relevant files in folder (csv only, skip images/zips)."""
    total = 0
    for f in os.listdir(folder):
        fp = os.path.realpath(os.path.join(folder, f))
        if os.path.isfile(fp) and f.endswith(".csv"):
            total += os.path.getsize(fp)
    return round(total / (1024 * 1024), 1)


def resolve_folder_symlinks(folder: str) -> tuple[str, bool]:
    """If folder contains symlinks, copy to temp dir with resolved files.
    Returns (path_to_use, is_temp)."""
    has_symlinks = any(
        os.path.islink(os.path.join(folder, f)) for f in os.listdir(folder)
    )
    if not has_symlinks:
        return folder, False

    # Copy with symlink resolution
    tmp = tempfile.mkdtemp(prefix="benchmark_")
    for f in os.listdir(folder):
        src = os.path.realpath(os.path.join(folder, f))
        dst = os.path.join(tmp, f)
        if os.path.isfile(src):
            shutil.copy2(src, dst)
    return tmp, True


async def measure_memory(cdp) -> Optional[float]:
    try:
        metrics = await cdp.send("Performance.getMetrics")
        for m in metrics.get("metrics", []):
            if m["name"] == "JSHeapUsedSize":
                return round(m["value"] / (1024 * 1024), 1)
    except Exception:
        pass
    return None


async def measure_fps(page) -> Optional[float]:
    try:
        return await page.evaluate("""() => new Promise(resolve => {
            let frames = 0;
            const start = performance.now();
            function tick() {
                frames++;
                if (performance.now() - start < 2000)
                    requestAnimationFrame(tick);
                else
                    resolve(Math.round(frames / ((performance.now() - start) / 1000)));
            }
            requestAnimationFrame(tick);
        })""")
    except Exception:
        return None


async def benchmark_single_molecule(
    page, cdp, ds: dict, base_url: str
) -> BenchmarkResult:
    """Benchmark a single molecule CSV file upload through the browser UI."""
    csv_path = os.path.realpath(ds["path"])
    size_mb = get_file_size_mb(csv_path)

    result = BenchmarkResult(
        platform="MERFISHeyes",
        dataset_name=ds["name"],
        data_type="single_molecule",
        size_tier=ds["size_tier"],
        n_items=ds["n_molecules"],
        file_size_mb=size_mb,
    )

    try:
        # Navigate to homepage
        await page.goto(base_url, wait_until="networkidle", timeout=30_000)

        # Navigate directly to single molecule mode via URL param
        await page.goto(f"{base_url}?mode=sm", wait_until="networkidle", timeout=30_000)
        await page.wait_for_timeout(500)

        # Find the MERSCOPE upload input
        input_id = f"sm-file-input-{ds['upload_type']}"
        file_input = page.locator(f"#{input_id}")

        # Start timing from file selection
        t0 = time.perf_counter()

        # Set the file on the input element
        await file_input.set_input_files(csv_path)

        print(f"      File set, waiting for processing and navigation...")

        # Wait for navigation to sm-viewer (can take very long for large files)
        timeout_ms = max(600_000, int(size_mb * 1000))  # at least 10 min, ~1s per MB
        try:
            await page.wait_for_url("**/sm-viewer/**", timeout=timeout_ms)
            await page.wait_for_selector("canvas", timeout=120_000)
            await page.wait_for_timeout(5000)  # let rendering stabilize

            result.load_time_s = round(time.perf_counter() - t0, 2)
            result.ui_responsive = True

            # Measure memory and FPS
            result.memory_peak_mb = await measure_memory(cdp)
            result.fps_after_load = await measure_fps(page)

        except Exception as e:
            elapsed = time.perf_counter() - t0
            result.load_time_s = round(elapsed, 2)
            if "timeout" in str(e).lower():
                result.crashed = True
                result.error = f"Timed out after {elapsed:.0f}s"
            else:
                result.crashed = True
                result.error = str(e)[:500]

    except Exception as e:
        result.crashed = True
        result.error = str(e)[:500]

    return result


async def benchmark_single_cell_merscope(
    page, cdp, ds: dict, base_url: str
) -> BenchmarkResult:
    """Benchmark a MERSCOPE folder upload through the browser UI."""
    folder = ds["folder"]
    size_mb = get_folder_size_mb(folder)

    result = BenchmarkResult(
        platform="MERFISHeyes",
        dataset_name=ds["name"],
        data_type="single_cell",
        size_tier=ds["size_tier"],
        n_items=ds.get("n_cells"),
        file_size_mb=size_mb,
    )

    # Resolve symlinks - Playwright can't handle symlinks in webkitdirectory inputs
    upload_folder, is_temp = resolve_folder_symlinks(folder)

    try:
        await page.goto(base_url, wait_until="networkidle", timeout=30_000)

        # Single cell mode is default. Find the MERSCOPE folder input.
        file_input = page.locator("#sc-file-input-merscope")

        t0 = time.perf_counter()

        # For webkitdirectory inputs, Playwright requires a directory path
        await file_input.set_input_files(upload_folder, timeout=60_000)

        print(f"      Folder set, waiting for processing...")

        # Wait for navigation to viewer
        timeout_ms = max(600_000, int(size_mb * 2000))
        try:
            await page.wait_for_url("**/viewer/**", timeout=timeout_ms)
            await page.wait_for_selector("canvas", timeout=120_000)
            await page.wait_for_timeout(5000)

            result.load_time_s = round(time.perf_counter() - t0, 2)
            result.ui_responsive = True

            result.memory_peak_mb = await measure_memory(cdp)
            result.fps_after_load = await measure_fps(page)

        except Exception as e:
            elapsed = time.perf_counter() - t0
            result.load_time_s = round(elapsed, 2)
            result.crashed = True
            result.error = str(e)[:500]

    except Exception as e:
        result.crashed = True
        result.error = str(e)[:500]
    finally:
        if is_temp and os.path.isdir(upload_folder):
            shutil.rmtree(upload_folder, ignore_errors=True)

    return result


async def benchmark_single_cell_xenium(
    page, cdp, ds: dict, base_url: str
) -> BenchmarkResult:
    """Benchmark a Xenium folder upload through the browser UI."""
    folder = ds["folder"]

    # Only count the relevant files (cells.csv.gz, cell_feature_matrix.h5, etc.)
    relevant_exts = {".csv", ".csv.gz", ".h5", ".parquet", ".json", ".xenium"}
    size_mb = 0
    files_to_upload = []
    for f in sorted(os.listdir(folder)):
        fp = os.path.join(folder, f)
        real_fp = os.path.realpath(fp)
        if os.path.isfile(real_fp):
            # Skip huge morphology/zip files that the adapter doesn't need
            if any(f.endswith(ext) for ext in [".ome.tif", ".zip", ".tar.gz", ".tar",
                                                 ".zarr.zip", ".html"]):
                continue
            files_to_upload.append(real_fp)
            size_mb += os.path.getsize(real_fp) / (1024 * 1024)

    size_mb = round(size_mb, 1)

    result = BenchmarkResult(
        platform="MERFISHeyes",
        dataset_name=ds["name"],
        data_type="single_cell",
        size_tier=ds["size_tier"],
        n_items=ds.get("n_cells"),
        file_size_mb=size_mb,
        notes=f"Relevant files only ({len(files_to_upload)} files, excluding morphology/zips)",
    )

    try:
        await page.goto(base_url, wait_until="networkidle", timeout=30_000)

        file_input = page.locator("#sc-file-input-xenium")

        real_folder = os.path.realpath(folder)

        t0 = time.perf_counter()

        # For webkitdirectory inputs, Playwright requires a directory path
        await file_input.set_input_files(real_folder, timeout=60_000)

        print(f"      Folder set ({size_mb:.0f}MB relevant files), waiting...")

        timeout_ms = max(600_000, int(size_mb * 3000))
        try:
            await page.wait_for_url("**/viewer/**", timeout=timeout_ms)
            await page.wait_for_selector("canvas", timeout=120_000)
            await page.wait_for_timeout(5000)

            result.load_time_s = round(time.perf_counter() - t0, 2)
            result.ui_responsive = True

            result.memory_peak_mb = await measure_memory(cdp)
            result.fps_after_load = await measure_fps(page)

        except Exception as e:
            elapsed = time.perf_counter() - t0
            result.load_time_s = round(elapsed, 2)
            result.crashed = True
            result.error = str(e)[:500]

    except Exception as e:
        result.crashed = True
        result.error = str(e)[:500]

    return result


def print_results(results: list[BenchmarkResult]):
    """Print formatted results table."""
    print("\n" + "=" * 100)
    print("  BENCHMARK RESULTS")
    print("=" * 100)

    for dtype in ["single_molecule", "single_cell"]:
        group = [r for r in results if r.data_type == dtype]
        if not group:
            continue

        label = "SINGLE MOLECULE" if dtype == "single_molecule" else "SINGLE CELL"
        print(f"\n  {label}")
        print("-" * 100)
        print(f"  {'Dataset':<45} {'Size':>8} {'Load':>10} {'Memory':>10} {'FPS':>6} {'Status':>10}")
        print("-" * 100)

        for r in group:
            load = f"{r.load_time_s:.1f}s" if r.load_time_s else "--"
            mem = f"{r.memory_peak_mb:.0f}MB" if r.memory_peak_mb else "--"
            fps = f"{r.fps_after_load:.0f}" if r.fps_after_load else "--"
            status = "CRASH" if r.crashed else "OK"
            size = f"{r.file_size_mb:.0f}MB"
            items = ""
            if r.n_items:
                if r.n_items >= 1_000_000:
                    items = f" ({r.n_items/1_000_000:.0f}M)"
                else:
                    items = f" ({r.n_items/1_000:.0f}K)"

            name = r.dataset_name[:45]
            print(f"  {name:<45} {size:>8} {load:>10} {mem:>10} {fps:>6} {status:>10}")
            if r.error:
                print(f"    ERROR: {r.error[:80]}")
            if r.notes:
                print(f"    NOTE: {r.notes[:80]}")

    print()


def save_results(results: list[BenchmarkResult], output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    path = os.path.join(output_dir, f"real_data_benchmark_{stamp}.json")
    with open(path, "w") as f:
        json.dump(
            {
                "timestamp": stamp,
                "test_data_root": TEST_DATA_ROOT,
                "results": [asdict(r) for r in results],
            },
            f,
            indent=2,
        )
    return path


async def async_main():
    parser = argparse.ArgumentParser(description="Benchmark MERFISHeyes with real test data")
    parser.add_argument("--base-url", default="http://localhost:3000")
    parser.add_argument("--output-dir", default="./benchmark_results")
    parser.add_argument("--only", choices=["sm", "sc", "sc-merscope", "sc-xenium"],
                        help="Only run specific benchmark type")
    parser.add_argument("--max-size", type=float, default=float("inf"),
                        help="Skip datasets larger than this many MB")
    parser.add_argument("--headless", action="store_true", help="Run in headless mode")
    parser.add_argument("--skip-large", action="store_true",
                        help="Skip datasets > 5GB (useful for quick tests)")
    args = parser.parse_args()

    from playwright.async_api import async_playwright

    print("=" * 70)
    print("  MERFISHeyes Real Data Benchmark")
    print("=" * 70)
    print(f"  Server:    {args.base_url}")
    print(f"  Data root: {TEST_DATA_ROOT}")
    print(f"  Max size:  {args.max_size}MB" if args.max_size < float("inf") else "  Max size:  unlimited")
    if args.only:
        print(f"  Filter:    {args.only}")
    print()

    results: list[BenchmarkResult] = []

    async with async_playwright() as pw:
        browser = await pw.chromium.launch(
            headless=args.headless,
            args=[
                "--disable-web-security",
                "--no-sandbox",
                "--js-flags=--max-old-space-size=8192",
                "--disable-dev-shm-usage",
            ]
        )

        # --- Single Molecule Benchmarks ---
        if args.only is None or args.only == "sm":
            print("\n--- SINGLE MOLECULE BENCHMARKS ---\n")
            for ds in SINGLE_MOLECULE_DATASETS:
                csv_path = os.path.realpath(ds["path"])
                if not os.path.exists(csv_path):
                    print(f"  SKIP (missing): {ds['name']}")
                    continue

                size_mb = get_file_size_mb(csv_path)
                if size_mb > args.max_size:
                    print(f"  SKIP (too large: {size_mb}MB): {ds['name']}")
                    continue
                if args.skip_large and size_mb > 5000:
                    print(f"  SKIP (--skip-large): {ds['name']}")
                    continue

                print(f"  Testing: {ds['name']} ({size_mb:.0f}MB)")

                # Fresh page for each test to avoid memory accumulation
                page = await browser.new_page()
                cdp = await page.context.new_cdp_session(page)
                await cdp.send("Performance.enable")

                try:
                    r = await benchmark_single_molecule(page, cdp, ds, args.base_url)
                    results.append(r)
                    status = "CRASH" if r.crashed else "OK"
                    load = f"{r.load_time_s:.1f}s" if r.load_time_s else "--"
                    mem = f"{r.memory_peak_mb:.0f}MB" if r.memory_peak_mb else "--"
                    print(f"    → {status}  load={load}  mem={mem}")
                    if r.error:
                        print(f"    → ERROR: {r.error[:120]}")
                finally:
                    await page.close()

        # --- Single Cell MERSCOPE Benchmarks ---
        if args.only is None or args.only in ("sc", "sc-merscope"):
            print("\n--- SINGLE CELL MERSCOPE BENCHMARKS ---\n")
            for ds in SINGLE_CELL_MERSCOPE_DATASETS:
                if not os.path.isdir(ds["folder"]):
                    print(f"  SKIP (missing): {ds['name']}")
                    continue

                size_mb = get_folder_size_mb(ds["folder"])
                if size_mb > args.max_size:
                    print(f"  SKIP (too large: {size_mb}MB): {ds['name']}")
                    continue

                print(f"  Testing: {ds['name']} ({size_mb:.0f}MB)")

                page = await browser.new_page()
                cdp = await page.context.new_cdp_session(page)
                await cdp.send("Performance.enable")

                try:
                    r = await benchmark_single_cell_merscope(page, cdp, ds, args.base_url)
                    results.append(r)
                    status = "CRASH" if r.crashed else "OK"
                    load = f"{r.load_time_s:.1f}s" if r.load_time_s else "--"
                    mem = f"{r.memory_peak_mb:.0f}MB" if r.memory_peak_mb else "--"
                    print(f"    → {status}  load={load}  mem={mem}")
                    if r.error:
                        print(f"    → ERROR: {r.error[:120]}")
                finally:
                    await page.close()

        # --- Single Cell Xenium Benchmarks ---
        if args.only is None or args.only in ("sc", "sc-xenium"):
            print("\n--- SINGLE CELL XENIUM BENCHMARKS ---\n")
            for ds in SINGLE_CELL_XENIUM_DATASETS:
                if not os.path.isdir(ds["folder"]):
                    print(f"  SKIP (missing): {ds['name']}")
                    continue

                # Rough size estimate of relevant files
                size_mb = 0
                for f in os.listdir(ds["folder"]):
                    fp = os.path.realpath(os.path.join(ds["folder"], f))
                    if os.path.isfile(fp) and not any(f.endswith(ext) for ext in [".ome.tif", ".zip", ".tar.gz", ".tar", ".zarr.zip", ".html"]):
                        size_mb += os.path.getsize(fp) / (1024 * 1024)

                if size_mb > args.max_size:
                    print(f"  SKIP (too large: {size_mb:.0f}MB): {ds['name']}")
                    continue
                if args.skip_large and size_mb > 5000:
                    print(f"  SKIP (--skip-large): {ds['name']}")
                    continue

                print(f"  Testing: {ds['name']} (~{size_mb:.0f}MB relevant files)")

                page = await browser.new_page()
                cdp = await page.context.new_cdp_session(page)
                await cdp.send("Performance.enable")

                try:
                    r = await benchmark_single_cell_xenium(page, cdp, ds, args.base_url)
                    results.append(r)
                    status = "CRASH" if r.crashed else "OK"
                    load = f"{r.load_time_s:.1f}s" if r.load_time_s else "--"
                    mem = f"{r.memory_peak_mb:.0f}MB" if r.memory_peak_mb else "--"
                    print(f"    → {status}  load={load}  mem={mem}")
                    if r.error:
                        print(f"    → ERROR: {r.error[:120]}")
                finally:
                    await page.close()

        await browser.close()

    if results:
        out_path = save_results(results, args.output_dir)
        print(f"\nResults saved to {out_path}")
        print_results(results)
    else:
        print("\nNo results collected.")


def main():
    asyncio.run(async_main())


if __name__ == "__main__":
    main()
