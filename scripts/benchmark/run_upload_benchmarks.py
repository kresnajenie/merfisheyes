#!/usr/bin/env python3
"""Benchmark upload speed for datasets that load successfully in the browser.

This measures:
  1. The time to complete the "Upload & Save" flow (presigned URLs → S3 upload → complete)
  2. Overall upload throughput (MB/s)

Prerequisites:
  - MERFISHeyes dev server running with valid S3 config (.env with AWS credentials)
  - Playwright installed
  - Datasets that successfully load in the browser (use run_browser_benchmarks.py first)

Usage:
    python run_upload_benchmarks.py --data-dir benchmark_data/synthetic
    python run_upload_benchmarks.py --data-dir benchmark_data/synthetic --filter-csv browser_results.csv
    python run_upload_benchmarks.py --types single_cell --filetypes h5ad
"""

import argparse
import asyncio
import csv
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (
    BenchmarkResult, DEV_SERVER_URL, BENCHMARK_TIMEOUT_S, DEFAULT_RESULTS_DIR,
    get_size_tier, get_git_version,
)
from run_browser_benchmarks import (
    BrowserBenchmarkRunner, append_result_csv, get_file_size_mb,
    load_existing_results, scan_single_cell, scan_single_molecule,
    list_folder_files,
)

# Upload-specific CSV columns
UPLOAD_CSV_COLUMNS = [
    "platform", "data_type", "size_tier",
    "n_cells", "n_molecules", "n_genes",
    "path", "file_size_mb", "filetype",
    "upload_time_s", "upload_throughput_mbps",
    "crashed", "error", "notes",
    "timestamp", "version",
]


def append_upload_csv(csv_path: str, row: dict):
    exists = os.path.exists(csv_path)
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=UPLOAD_CSV_COLUMNS)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def filter_successful(filter_csv: str) -> set[str]:
    """Return set of paths that loaded successfully (crashed=False)."""
    paths = set()
    if not os.path.exists(filter_csv):
        return paths
    with open(filter_csv, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("crashed", "True").lower() == "false":
                paths.add(row.get("path", ""))
    return paths


class UploadBenchmarkRunner:
    """Measures upload speed by triggering the Upload & Save flow in the browser."""

    def __init__(self, base_url: str, timeout_s: int, headless: bool = False):
        self.base_url = base_url
        self.timeout_ms = timeout_s * 1000
        self.headless = headless

    async def benchmark_upload_sc(
        self, path: str, filetype: str, n_cells: int, n_genes: int,
    ) -> dict:
        from playwright.async_api import async_playwright

        row = {
            "platform": "MERFISHeyes",
            "data_type": "single_cell",
            "size_tier": get_size_tier(n_cells, "single_cell"),
            "n_cells": n_cells, "n_molecules": None,
            "n_genes": n_genes, "path": path,
            "file_size_mb": round(get_file_size_mb(path), 2),
            "filetype": filetype,
            "upload_time_s": None, "upload_throughput_mbps": None,
            "crashed": False, "error": None, "notes": "",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": get_git_version(),
        }

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=self.headless)
            page = await browser.new_page()

            try:
                await page.goto(self.base_url, wait_until="networkidle", timeout=30_000)

                # Upload the file first (load into browser)
                if filetype == "h5ad":
                    async with page.expect_file_chooser() as fc_info:
                        await page.get_by_text("H5AD FileSingle .h5ad file").click()
                    fc = await fc_info.value
                    await fc.set_files(path)
                elif filetype == "xenium":
                    files = list_folder_files(path)
                    async with page.expect_file_chooser() as fc_info:
                        await page.get_by_text("Xenium FolderSelect Xenium output folder").click()
                    fc = await fc_info.value
                    await fc.set_files(files)
                elif filetype == "merscope":
                    files = list_folder_files(path)
                    async with page.expect_file_chooser() as fc_info:
                        await page.get_by_text("Merscope FolderSelect Merscope output folder").click()
                    fc = await fc_info.value
                    await fc.set_files(files)

                # Wait for viewer to load
                await page.wait_for_url("**/viewer/**", timeout=self.timeout_ms)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(2000)

                # Click "Upload & Save" button
                upload_btn = page.get_by_role("button", name="Upload & Save")
                if await upload_btn.count() == 0:
                    upload_btn = page.get_by_text("Upload & Save")
                if await upload_btn.count() == 0:
                    row["notes"] = "Upload button not found"
                else:
                    t0 = time.perf_counter()
                    await upload_btn.click()

                    # Wait for upload modal / confirmation
                    # The modal has a "Upload" or "Confirm" button
                    confirm = page.get_by_role("button", name="Upload")
                    if await confirm.count() > 0:
                        await confirm.click()

                    # Wait for upload to complete (URL changes to /viewer/{id})
                    await page.wait_for_url("**/viewer/*", timeout=self.timeout_ms)
                    # Wait for success toast or completion indicator
                    await page.wait_for_timeout(2000)

                    upload_time = round(time.perf_counter() - t0, 2)
                    row["upload_time_s"] = upload_time
                    if row["file_size_mb"] > 0 and upload_time > 0:
                        row["upload_throughput_mbps"] = round(row["file_size_mb"] / upload_time, 2)

            except Exception as e:
                row["crashed"] = True
                row["error"] = str(e)[:500]
            finally:
                await browser.close()

        return row

    async def benchmark_upload_sm(
        self, path: str, filetype: str, n_molecules: int, n_genes: int,
    ) -> dict:
        from playwright.async_api import async_playwright

        row = {
            "platform": "MERFISHeyes",
            "data_type": "single_molecule",
            "size_tier": get_size_tier(n_molecules, "single_molecule"),
            "n_cells": None, "n_molecules": n_molecules,
            "n_genes": n_genes, "path": path,
            "file_size_mb": round(get_file_size_mb(path), 2),
            "filetype": filetype,
            "upload_time_s": None, "upload_throughput_mbps": None,
            "crashed": False, "error": None, "notes": "",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": get_git_version(),
        }

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=self.headless)
            page = await browser.new_page()

            try:
                await page.goto(self.base_url, wait_until="networkidle", timeout=30_000)
                await page.get_by_role("button", name="single molecule", exact=True).click()
                await page.wait_for_timeout(500)

                async with page.expect_file_chooser() as fc_info:
                    await page.get_by_text("Xenium Parquet/CSVSelect .parquet or .csv file").click()
                fc = await fc_info.value
                await fc.set_files(path)

                await page.wait_for_url("**/sm-viewer/**", timeout=self.timeout_ms)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(2000)

                upload_btn = page.get_by_role("button", name="Upload & Save")
                if await upload_btn.count() == 0:
                    upload_btn = page.get_by_text("Upload & Save")
                if await upload_btn.count() > 0:
                    t0 = time.perf_counter()
                    await upload_btn.click()

                    confirm = page.get_by_role("button", name="Upload")
                    if await confirm.count() > 0:
                        await confirm.click()

                    await page.wait_for_url("**/sm-viewer/*", timeout=self.timeout_ms)
                    await page.wait_for_timeout(2000)

                    upload_time = round(time.perf_counter() - t0, 2)
                    row["upload_time_s"] = upload_time
                    if row["file_size_mb"] > 0 and upload_time > 0:
                        row["upload_throughput_mbps"] = round(row["file_size_mb"] / upload_time, 2)
                else:
                    row["notes"] = "Upload button not found"

            except Exception as e:
                row["crashed"] = True
                row["error"] = str(e)[:500]
            finally:
                await browser.close()

        return row


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
async def run_upload_benchmarks(args):
    import urllib.request
    try:
        urllib.request.urlopen(args.url, timeout=5)
    except Exception:
        print(f"ERROR: Dev server not running at {args.url}")
        sys.exit(1)

    runner = UploadBenchmarkRunner(args.url, args.timeout, args.headless)

    # Optionally filter to only datasets that loaded successfully
    allowed_paths = None
    if args.filter_csv:
        allowed_paths = filter_successful(args.filter_csv)
        print(f"Filtering to {len(allowed_paths)} successful datasets from {args.filter_csv}")

    tested = load_existing_results(args.output) if args.resume else set()
    completed = 0

    if "single_cell" in args.types:
        sc_datasets = scan_single_cell(args.data_dir, args.sc_filetypes)
        for ds in sc_datasets:
            if ds["path"] in tested:
                continue
            if allowed_paths and ds["path"] not in allowed_paths:
                continue
            print(f"  Upload SC: {ds['filetype']} {ds['n_cells']:,}c × {ds['n_genes']:,}g")
            row = await runner.benchmark_upload_sc(
                ds["path"], ds["filetype"], ds["n_cells"], ds["n_genes"],
            )
            append_upload_csv(args.output, row)
            completed += 1
            if row["crashed"]:
                print(f"    CRASHED: {row['error']}")
            else:
                print(f"    upload={row['upload_time_s']}s  throughput={row['upload_throughput_mbps']} MB/s")

    if "single_molecule" in args.types:
        sm_datasets = scan_single_molecule(args.data_dir, args.sm_filetypes)
        for ds in sm_datasets:
            if ds["path"] in tested:
                continue
            if allowed_paths and ds["path"] not in allowed_paths:
                continue
            print(f"  Upload SM: {ds['filetype']} {ds['n_molecules']:,}m × {ds['n_genes']:,}g")
            row = await runner.benchmark_upload_sm(
                ds["path"], ds["filetype"], ds["n_molecules"], ds["n_genes"],
            )
            append_upload_csv(args.output, row)
            completed += 1
            if row["crashed"]:
                print(f"    CRASHED: {row['error']}")
            else:
                print(f"    upload={row['upload_time_s']}s  throughput={row['upload_throughput_mbps']} MB/s")

    print(f"\nDone. {completed} upload tests.")
    print(f"Results: {args.output}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark upload speed")
    parser.add_argument("--data-dir", default=DEFAULT_SYNTH_DIR)
    parser.add_argument("--output", default=os.path.join(DEFAULT_RESULTS_DIR, "upload_benchmarks.csv"))
    parser.add_argument("--filter-csv", default=None,
                        help="Only test datasets that passed in this browser benchmark CSV")
    parser.add_argument("--url", default=DEV_SERVER_URL)
    parser.add_argument("--timeout", type=int, default=BENCHMARK_TIMEOUT_S)
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--types", nargs="+", default=["single_cell", "single_molecule"],
                        choices=["single_cell", "single_molecule"])
    parser.add_argument("--sc-filetypes", nargs="+", default=["h5ad", "xenium", "merscope"])
    parser.add_argument("--sm-filetypes", nargs="+", default=["parquet", "csv"])
    args = parser.parse_args()
    asyncio.run(run_upload_benchmarks(args))


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
