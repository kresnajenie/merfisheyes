#!/usr/bin/env python3
"""Run browser-based benchmarks against MERFISHeyes using Playwright.

Loads each synthetic (or real) dataset into the browser and measures:
  load_time_s, gene_query_time_s, memory_peak_mb, fps_after_load, ui_responsive, crashed

Prerequisites:
  1. Generate data:     python generate_single_cell.py / generate_single_molecule.py
  2. Install Playwright: pip install playwright && python -m playwright install chromium
  3. Start dev server:   npm run dev  (http://localhost:3000)

Usage:
    python run_browser_benchmarks.py
    python run_browser_benchmarks.py --data-dir benchmark_data/synthetic --output results.csv
    python run_browser_benchmarks.py --types single_cell --filetypes h5ad --max-cells 500000
    python run_browser_benchmarks.py --resume   # skip already-tested combinations
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
    BenchmarkResult, DEV_SERVER_URL, BENCHMARK_TIMEOUT_S,
    GENE_TO_QUERY, SM_GENE_TO_QUERY, DEFAULT_RESULTS_DIR,
    get_size_tier, get_git_version,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def get_file_size_mb(path: str) -> float:
    if os.path.isfile(path):
        return os.path.getsize(path) / (1024 * 1024)
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            total += os.path.getsize(os.path.join(root, f))
    return total / (1024 * 1024)


def list_folder_files(folder: str) -> list[str]:
    """Recursively list all files in a folder (for directory upload)."""
    files = []
    for root, _, names in os.walk(folder):
        for name in names:
            files.append(os.path.join(root, name))
    return sorted(files)


def load_existing_results(csv_path: str) -> set[str]:
    """Load already-tested paths from an existing CSV for resume support."""
    tested = set()
    if not os.path.exists(csv_path):
        return tested
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tested.add(row.get("path", ""))
    return tested


def append_result_csv(csv_path: str, result: BenchmarkResult):
    """Append one result row to CSV. Creates file with header if needed."""
    exists = os.path.exists(csv_path)
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=BenchmarkResult.CSV_COLUMNS)
        if not exists:
            writer.writeheader()
        writer.writerow(result.to_dict())


# ---------------------------------------------------------------------------
# Browser benchmark runner
# ---------------------------------------------------------------------------
class BrowserBenchmarkRunner:
    def __init__(self, base_url: str, timeout_s: int, headless: bool = False):
        self.base_url = base_url
        self.timeout_ms = timeout_s * 1000
        self.headless = headless

    def is_server_up(self) -> bool:
        import urllib.request
        try:
            urllib.request.urlopen(self.base_url, timeout=5)
            return True
        except Exception:
            return False

    # -- single cell -------------------------------------------------------

    async def run_single_cell(
        self, path: str, filetype: str, n_cells: int, n_genes: int,
    ) -> BenchmarkResult:
        result = BenchmarkResult(
            data_type="single_cell",
            size_tier=get_size_tier(n_cells, "single_cell"),
            n_cells=n_cells,
            n_genes=n_genes,
            path=path,
            file_size_mb=round(get_file_size_mb(path), 2),
            filetype=filetype,
        )

        from playwright.async_api import async_playwright

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=self.headless)
            page = await browser.new_page()
            cdp = await page.context.new_cdp_session(page)
            await cdp.send("Performance.enable")

            try:
                await page.goto(self.base_url, wait_until="networkidle", timeout=30_000)

                t0 = time.perf_counter()

                if filetype == "h5ad":
                    await self._upload_h5ad(page, path)
                elif filetype == "xenium":
                    await self._upload_xenium_folder(page, path)
                elif filetype == "merscope":
                    await self._upload_merscope_folder(page, path)

                # Wait for viewer
                await page.wait_for_url("**/viewer/**", timeout=self.timeout_ms)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(3000)

                result.load_time_s = round(time.perf_counter() - t0, 2)
                result.ui_responsive = True
                result.memory_peak_mb = await self._measure_memory(cdp)
                result.fps_after_load = await self._measure_fps(page)
                result.gene_query_time_s = await self._measure_gene_query_sc(page, GENE_TO_QUERY)

            except Exception as e:
                result.crashed = True
                result.error = str(e)[:500]
                if result.load_time_s is None:
                    result.load_time_s = round(time.perf_counter() - t0, 2)
            finally:
                await browser.close()

        return result

    # -- single molecule ---------------------------------------------------

    async def run_single_molecule(
        self, path: str, filetype: str, n_molecules: int, n_genes: int,
    ) -> BenchmarkResult:
        result = BenchmarkResult(
            data_type="single_molecule",
            size_tier=get_size_tier(n_molecules, "single_molecule"),
            n_molecules=n_molecules,
            n_genes=n_genes,
            path=path,
            file_size_mb=round(get_file_size_mb(path), 2),
            filetype=filetype,
        )

        from playwright.async_api import async_playwright

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=self.headless)
            page = await browser.new_page()
            cdp = await page.context.new_cdp_session(page)
            await cdp.send("Performance.enable")

            try:
                t0 = time.perf_counter()
                await self._upload_sm_file(page, path)

                await page.wait_for_url("**/sm-viewer/**", timeout=self.timeout_ms)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(3000)

                result.load_time_s = round(time.perf_counter() - t0, 2)
                result.ui_responsive = True
                result.memory_peak_mb = await self._measure_memory(cdp)
                result.fps_after_load = await self._measure_fps(page)
                result.gene_query_time_s = await self._measure_gene_query_sm(page, SM_GENE_TO_QUERY)

            except Exception as e:
                result.crashed = True
                result.error = str(e)[:500]
                if result.load_time_s is None:
                    result.load_time_s = round(time.perf_counter() - t0, 2)
            finally:
                await browser.close()

        return result

    # -- upload helpers ----------------------------------------------------

    async def _upload_h5ad(self, page, path: str):
        # Use .first to avoid strict mode violation when duplicate IDs exist
        await page.locator("#sc-file-input-h5ad").first.set_input_files(path)

    async def _upload_xenium_folder(self, page, folder: str):
        await page.locator("#sc-file-input-xenium").first.set_input_files(folder, timeout=60_000)

    async def _upload_merscope_folder(self, page, folder: str):
        await page.locator("#sc-file-input-merscope").first.set_input_files(folder, timeout=60_000)

    async def _upload_sm_file(self, page, path: str):
        await page.goto(f"{self.base_url}?mode=sm", wait_until="networkidle", timeout=30_000)
        await page.wait_for_timeout(500)
        await page.locator("#sm-file-input-xenium").first.set_input_files(path)

    # -- measurement helpers -----------------------------------------------

    @staticmethod
    async def _measure_memory(cdp) -> float | None:
        try:
            metrics = await cdp.send("Performance.getMetrics")
            for m in metrics.get("metrics", []):
                if m["name"] == "JSHeapUsedSize":
                    return round(m["value"] / (1024 * 1024), 1)
        except Exception:
            pass
        return None

    @staticmethod
    async def _measure_fps(page) -> float | None:
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

    @staticmethod
    async def _measure_gene_query_sc(page, gene_name: str) -> float | None:
        """Measure gene query time for single cell viewer."""
        try:
            btn = page.get_by_role("button", name="Gene")
            if await btn.count() == 0:
                return None
            await btn.click()
            await page.wait_for_timeout(400)

            search = page.get_by_role("textbox", name="Search gene")
            if await search.count() == 0:
                return None
            await search.fill(gene_name)
            await page.wait_for_timeout(300)

            radio = page.get_by_role("radio", name=gene_name)
            if await radio.count() == 0:
                # Try case-insensitive / partial match
                radio = page.locator(f'input[type="radio"]').first
                if await radio.count() == 0:
                    return None

            t0 = time.perf_counter()
            await radio.click()
            await page.wait_for_timeout(2000)
            return round(time.perf_counter() - t0, 2)
        except Exception:
            return None

    @staticmethod
    async def _measure_gene_query_sm(page, gene_name: str) -> float | None:
        """Measure gene query time for single molecule viewer (select biggest gene)."""
        try:
            # Single molecule viewer has checkbox-based gene selection
            # Try to find the gene in the control panel
            checkbox = page.get_by_label(gene_name, exact=True)
            if await checkbox.count() == 0:
                # Try finding via text
                checkbox = page.locator(f'text="{gene_name}"').first
                if await checkbox.count() == 0:
                    return None

            t0 = time.perf_counter()
            await checkbox.click()
            await page.wait_for_timeout(3000)  # gene loading from local data
            return round(time.perf_counter() - t0, 2)
        except Exception:
            return None


# ---------------------------------------------------------------------------
# Scan data directory for test files
# ---------------------------------------------------------------------------
def scan_single_cell(data_dir: str, filetypes: list[str]) -> list[dict]:
    """Scan for single cell datasets, either from manifest or directory structure."""
    datasets = []

    # Try manifest first
    manifest_path = os.path.join(data_dir, "single_cell", "manifest_single_cell.json")
    if not os.path.exists(manifest_path):
        manifest_path = os.path.join(data_dir, "manifest_single_cell.json")
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
        for ds in manifest.get("datasets", []):
            if ds["filetype"] in filetypes and os.path.exists(ds["path"]):
                datasets.append(ds)
        return datasets

    # Fallback: scan directory structure
    sc_dir = os.path.join(data_dir, "single_cell")
    if not os.path.exists(sc_dir):
        sc_dir = data_dir

    for ft in filetypes:
        ft_dir = os.path.join(sc_dir, ft)
        if not os.path.isdir(ft_dir):
            continue
        for entry in sorted(os.listdir(ft_dir)):
            path = os.path.join(ft_dir, entry)
            # Parse n_cells and n_genes from filename: sc_{cells}c_{genes}g
            parts = entry.replace(".h5ad", "").replace(".parquet", "").replace(".csv", "")
            try:
                tokens = parts.split("_")
                n_cells = int(tokens[1].rstrip("c"))
                n_genes = int(tokens[2].rstrip("g"))
            except (IndexError, ValueError):
                continue
            datasets.append({
                "n_cells": n_cells,
                "n_genes": n_genes,
                "filetype": ft,
                "path": path,
                "file_size_mb": round(get_file_size_mb(path), 2),
            })

    return datasets


def scan_single_molecule(data_dir: str, filetypes: list[str]) -> list[dict]:
    """Scan for single molecule datasets."""
    datasets = []

    manifest_path = os.path.join(data_dir, "single_molecule", "manifest_single_molecule.json")
    if not os.path.exists(manifest_path):
        manifest_path = os.path.join(data_dir, "manifest_single_molecule.json")
    if os.path.exists(manifest_path):
        with open(manifest_path) as f:
            manifest = json.load(f)
        for ds in manifest.get("datasets", []):
            if ds["filetype"] in filetypes and os.path.exists(ds["path"]):
                datasets.append(ds)
        return datasets

    sm_dir = os.path.join(data_dir, "single_molecule")
    if not os.path.exists(sm_dir):
        sm_dir = data_dir

    for ft in filetypes:
        ft_dir = os.path.join(sm_dir, ft)
        if not os.path.isdir(ft_dir):
            continue
        for entry in sorted(os.listdir(ft_dir)):
            path = os.path.join(ft_dir, entry)
            parts = entry.replace(f".{ft}", "")
            try:
                tokens = parts.split("_")
                n_mol = int(tokens[1].rstrip("m"))
                n_genes = int(tokens[2].rstrip("g"))
            except (IndexError, ValueError):
                continue
            datasets.append({
                "n_molecules": n_mol,
                "n_genes": n_genes,
                "filetype": ft,
                "path": path,
                "file_size_mb": round(get_file_size_mb(path), 2),
            })

    return datasets


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
async def run_benchmarks(args):
    runner = BrowserBenchmarkRunner(
        base_url=args.url,
        timeout_s=args.timeout,
        headless=args.headless,
    )

    if not runner.is_server_up():
        print(f"ERROR: Dev server not running at {args.url}")
        print("Start it with:  npm run dev")
        sys.exit(1)

    print(f"MERFISHeyes browser benchmark")
    print(f"  Server:  {args.url}")
    print(f"  Data:    {args.data_dir}")
    print(f"  Output:  {args.output}")
    print(f"  Timeout: {args.timeout}s per test")
    print()

    tested = load_existing_results(args.output) if args.resume else set()
    if tested:
        print(f"  Resuming: {len(tested)} tests already complete\n")

    total_tests = 0
    completed = 0
    errors = 0

    # Single cell tests
    if "single_cell" in args.types:
        sc_datasets = scan_single_cell(args.data_dir, args.sc_filetypes)
        # Apply filters
        if args.max_cells:
            sc_datasets = [d for d in sc_datasets if d["n_cells"] <= args.max_cells]
        if args.max_genes:
            sc_datasets = [d for d in sc_datasets if d["n_genes"] <= args.max_genes]

        print(f"Single cell datasets: {len(sc_datasets)}")
        for i, ds in enumerate(sc_datasets):
            if ds["path"] in tested:
                print(f"  [{i+1}/{len(sc_datasets)}] SKIP (already tested): {ds['path']}")
                continue

            print(f"  [{i+1}/{len(sc_datasets)}] {ds['filetype']} — "
                  f"{ds['n_cells']:,} cells × {ds['n_genes']:,} genes "
                  f"({ds['file_size_mb']:.1f} MB)")

            result = await runner.run_single_cell(
                path=ds["path"],
                filetype=ds["filetype"],
                n_cells=ds["n_cells"],
                n_genes=ds["n_genes"],
            )
            append_result_csv(args.output, result)
            total_tests += 1

            if result.crashed:
                errors += 1
                print(f"    CRASHED: {result.error}")
            else:
                completed += 1
                print(f"    load={result.load_time_s}s  gene_query={result.gene_query_time_s}s  "
                      f"mem={result.memory_peak_mb}MB  fps={result.fps_after_load}")

    # Single molecule tests
    if "single_molecule" in args.types:
        sm_datasets = scan_single_molecule(args.data_dir, args.sm_filetypes)
        if args.max_molecules:
            sm_datasets = [d for d in sm_datasets if d["n_molecules"] <= args.max_molecules]
        if args.max_genes:
            sm_datasets = [d for d in sm_datasets if d["n_genes"] <= args.max_genes]

        print(f"\nSingle molecule datasets: {len(sm_datasets)}")
        for i, ds in enumerate(sm_datasets):
            if ds["path"] in tested:
                print(f"  [{i+1}/{len(sm_datasets)}] SKIP (already tested): {ds['path']}")
                continue

            print(f"  [{i+1}/{len(sm_datasets)}] {ds['filetype']} — "
                  f"{ds['n_molecules']:,} molecules × {ds['n_genes']:,} genes "
                  f"({ds['file_size_mb']:.1f} MB)")

            result = await runner.run_single_molecule(
                path=ds["path"],
                filetype=ds["filetype"],
                n_molecules=ds["n_molecules"],
                n_genes=ds["n_genes"],
            )
            append_result_csv(args.output, result)
            total_tests += 1

            if result.crashed:
                errors += 1
                print(f"    CRASHED: {result.error}")
            else:
                completed += 1
                print(f"    load={result.load_time_s}s  gene_query={result.gene_query_time_s}s  "
                      f"mem={result.memory_peak_mb}MB  fps={result.fps_after_load}")

    print(f"\nDone. {completed} passed, {errors} crashed out of {total_tests} tests.")
    print(f"Results: {args.output}")


def main():
    parser = argparse.ArgumentParser(description="Run browser benchmarks against MERFISHeyes")
    parser.add_argument("--data-dir", default=DEFAULT_SYNTH_DIR,
                        help="Directory containing test data")
    parser.add_argument("--output", default=os.path.join(DEFAULT_RESULTS_DIR, "browser_benchmarks.csv"),
                        help="Output CSV path")
    parser.add_argument("--url", default=DEV_SERVER_URL,
                        help="MERFISHeyes dev server URL")
    parser.add_argument("--timeout", type=int, default=BENCHMARK_TIMEOUT_S,
                        help="Timeout per test in seconds")
    parser.add_argument("--headless", action="store_true",
                        help="Run browser in headless mode")
    parser.add_argument("--resume", action="store_true",
                        help="Skip tests already present in output CSV")
    parser.add_argument("--types", nargs="+", default=["single_cell", "single_molecule"],
                        choices=["single_cell", "single_molecule"],
                        help="Data types to test")
    parser.add_argument("--sc-filetypes", nargs="+", default=["h5ad", "xenium", "merscope"],
                        help="Single cell file types to test")
    parser.add_argument("--sm-filetypes", nargs="+", default=["parquet", "csv"],
                        help="Single molecule file types to test")
    parser.add_argument("--max-cells", type=int, default=None,
                        help="Skip single cell tests with more cells than this")
    parser.add_argument("--max-molecules", type=int, default=None,
                        help="Skip single molecule tests with more molecules than this")
    parser.add_argument("--max-genes", type=int, default=None,
                        help="Skip tests with more genes than this")
    args = parser.parse_args()

    asyncio.run(run_benchmarks(args))


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
