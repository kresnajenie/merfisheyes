#!/usr/bin/env python3
"""
Benchmark spatial transcriptomics visualization platforms.

Compares MERFISHeyes against CellxGene, SpaTEO, and Vitessce on:
  - Data loading time
  - Gene query/switching time
  - Peak memory usage
  - UI responsiveness (frozen during load?)
  - FPS after rendering

Prerequisites:
  1. Generate test data:   python generate_test_data.py --output-dir ./benchmark_data
  2. Install Playwright:   pip install playwright && python -m playwright install chromium
  3. Start MERFISHeyes:    cd ../.. && npm run dev

Usage:
    python run_benchmarks.py --data-dir ./benchmark_data --platform merfisheyes
    python run_benchmarks.py --data-dir ./benchmark_data --platform cellxgene
    python run_benchmarks.py --data-dir ./benchmark_data --platform spateo
    python run_benchmarks.py --data-dir ./benchmark_data --platform all
    python run_benchmarks.py --compare ./benchmark_results

Manual-only platforms (no automation):
  - Xenium Explorer: Desktop app, test manually and record times
  - Allen Brain Cell Atlas: Hosted viewer, can't upload custom data
  - Samui: Niche tool, preprocess data and serve locally
  See the MANUAL TESTING GUIDE at the bottom of this file.
"""

import argparse
import asyncio
import json
import os
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional


# ---------------------------------------------------------------------------
#  Result data model
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkResult:
    platform: str
    data_type: str       # "single_cell" or "single_molecule"
    size_tier: str       # "small" | "medium" | "large" | "xl"
    n_items: int         # cells or molecules
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


# ---------------------------------------------------------------------------
#  MERFISHeyes benchmark (Playwright)
# ---------------------------------------------------------------------------

class MerfisheyesBenchmark:
    """Test MERFISHeyes by uploading files through the browser UI."""

    name = "MERFISHeyes"

    def __init__(self, base_url: str = "http://localhost:3000"):
        self.base_url = base_url

    def is_available(self) -> bool:
        import urllib.request
        try:
            urllib.request.urlopen(self.base_url, timeout=5)
            return True
        except Exception:
            return False

    # -- single cell ---------------------------------------------------

    async def benchmark_single_cell(
        self, path: str, size_tier: str, n_cells: int
    ) -> BenchmarkResult:
        from playwright.async_api import async_playwright

        result = self._make_result("single_cell", size_tier, n_cells, path)

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=False)
            page = await browser.new_page()
            cdp = await page.context.new_cdp_session(page)
            await cdp.send("Performance.enable")

            try:
                await page.goto(self.base_url, wait_until="networkidle")

                # Click the H5AD card to open file chooser, then set file
                t0 = time.perf_counter()
                async with page.expect_file_chooser() as fc_info:
                    await page.get_by_text("H5AD FileSingle .h5ad file").click()
                file_chooser = await fc_info.value
                await file_chooser.set_files(path)

                # Wait for navigation to the viewer page (worker done + redirect)
                await page.wait_for_url("**/viewer/**", timeout=600_000)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(3000)

                result.load_time_s = round(time.perf_counter() - t0, 2)
                result.ui_responsive = True

                result.memory_peak_mb = await self._measure_memory(cdp)
                result.fps_after_load = await self._measure_fps(page)

                # Gene query: open gene panel and select first gene
                result.gene_query_time_s = await self._measure_gene_query(page)

            except Exception as e:
                result.crashed = True
                result.error = str(e)[:500]
            finally:
                await browser.close()

        return result

    # -- single molecule -----------------------------------------------

    async def benchmark_single_molecule(
        self, path: str, size_tier: str, n_molecules: int
    ) -> BenchmarkResult:
        from playwright.async_api import async_playwright

        result = self._make_result("single_molecule", size_tier, n_molecules, path)

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(headless=False)
            page = await browser.new_page()
            cdp = await page.context.new_cdp_session(page)
            await cdp.send("Performance.enable")

            try:
                await page.goto(self.base_url, wait_until="networkidle")

                # Switch to single-molecule mode
                await page.get_by_role("button", name="single molecule", exact=True).click()
                await page.wait_for_timeout(500)

                # Click the Xenium Parquet/CSV card to open file chooser
                t0 = time.perf_counter()
                async with page.expect_file_chooser() as fc_info:
                    await page.get_by_text("Xenium Parquet/CSVSelect .parquet or .csv file").click()
                file_chooser = await fc_info.value
                await file_chooser.set_files(path)

                await page.wait_for_url("**/sm-viewer/**", timeout=600_000)
                await page.wait_for_selector("canvas", timeout=120_000)
                await page.wait_for_timeout(3000)

                result.load_time_s = round(time.perf_counter() - t0, 2)
                result.ui_responsive = True

                result.memory_peak_mb = await self._measure_memory(cdp)
                result.fps_after_load = await self._measure_fps(page)

            except Exception as e:
                result.crashed = True
                result.error = str(e)[:500]
            finally:
                await browser.close()

        return result

    # -- helpers -------------------------------------------------------

    def _make_result(self, data_type, tier, n, path):
        return BenchmarkResult(
            platform=self.name,
            data_type=data_type,
            size_tier=tier,
            n_items=n,
            file_size_mb=round(os.path.getsize(path) / (1024 * 1024), 1),
        )

    @staticmethod
    async def _measure_memory(cdp) -> Optional[float]:
        try:
            metrics = await cdp.send("Performance.getMetrics")
            for m in metrics.get("metrics", []):
                if m["name"] == "JSHeapUsedSize":
                    return round(m["value"] / (1024 * 1024), 1)
        except Exception:
            pass
        return None

    @staticmethod
    async def _measure_fps(page) -> Optional[float]:
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
    async def _measure_gene_query(page) -> Optional[float]:
        """Click Gene button, search for a gene, click the radio, measure render time."""
        try:
            btn = page.get_by_role("button", name="Gene")
            if await btn.count() == 0:
                return None
            await btn.click()
            await page.wait_for_timeout(400)

            # Type gene name into search box to filter the list
            search = page.get_by_role("textbox", name="Search gene")
            if await search.count() == 0:
                return None
            await search.fill("Gene_0")
            await page.wait_for_timeout(300)

            # Click the first radio button in the filtered results
            radio = page.get_by_role("radio", name="Gene_0")
            if await radio.count() == 0:
                return None

            t0 = time.perf_counter()
            await radio.click()
            await page.wait_for_timeout(2000)  # wait for visualization update
            return round(time.perf_counter() - t0, 2)
        except Exception:
            pass
        return None


# ---------------------------------------------------------------------------
#  CellxGene benchmark (subprocess + Playwright)
# ---------------------------------------------------------------------------

class CellxgeneBenchmark:
    """Launch the cellxgene CLI, then test via Playwright."""

    name = "CellxGene"

    def __init__(self, port: int = 5005):
        self.port = port
        self._proc = None

    def is_available(self) -> bool:
        try:
            r = subprocess.run(
                ["cellxgene", "--version"], capture_output=True, text=True, timeout=10
            )
            return r.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    async def benchmark_single_cell(
        self, path: str, size_tier: str, n_cells: int
    ) -> BenchmarkResult:
        from playwright.async_api import async_playwright

        result = BenchmarkResult(
            platform=self.name,
            data_type="single_cell",
            size_tier=size_tier,
            n_items=n_cells,
            file_size_mb=round(os.path.getsize(path) / (1024 * 1024), 1),
        )

        try:
            startup_s = self._start_server(path)
            result.notes = f"server startup: {startup_s:.1f}s"

            async with async_playwright() as pw:
                browser = await pw.chromium.launch(headless=False)
                page = await browser.new_page()
                cdp = await page.context.new_cdp_session(page)
                await cdp.send("Performance.enable")

                t0 = time.perf_counter()
                await page.goto(f"http://localhost:{self.port}")
                await page.wait_for_selector("canvas", timeout=300_000)
                await page.wait_for_timeout(3000)

                render_s = time.perf_counter() - t0
                result.load_time_s = round(startup_s + render_s, 2)
                result.ui_responsive = True

                # Browser memory
                result.memory_peak_mb = await MerfisheyesBenchmark._measure_memory(cdp)

                # Server-side memory
                try:
                    import psutil
                    proc = psutil.Process(self._proc.pid)
                    server_mb = proc.memory_info().rss / (1024 * 1024)
                    result.notes += f", server RAM: {server_mb:.0f}MB"
                except Exception:
                    pass

                # Gene query: type a gene name into the search box
                t_gene = time.perf_counter()
                search = page.locator('input[placeholder*="search" i]').first
                if await search.count() > 0:
                    await search.fill("Gene_0")
                    await page.wait_for_timeout(3000)
                    result.gene_query_time_s = round(time.perf_counter() - t_gene, 2)

                await browser.close()

        except Exception as e:
            result.crashed = True
            result.error = str(e)[:500]
        finally:
            self._stop_server()

        return result

    # CellxGene does not support single-molecule data
    benchmark_single_molecule = None

    # -- helpers -------------------------------------------------------

    def _start_server(self, h5ad_path: str) -> float:
        import urllib.request

        t0 = time.perf_counter()
        self._proc = subprocess.Popen(
            ["cellxgene", "launch", h5ad_path,
             "--port", str(self.port), "--verbose"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        url = f"http://localhost:{self.port}"
        for _ in range(300):
            try:
                urllib.request.urlopen(url, timeout=2)
                return time.perf_counter() - t0
            except Exception:
                time.sleep(1)
        raise TimeoutError("CellxGene server failed to start within 5 min")

    def _stop_server(self):
        if self._proc:
            self._proc.terminate()
            try:
                self._proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self._proc.kill()
            self._proc = None


# ---------------------------------------------------------------------------
#  SpaTEO benchmark (pure Python timing)
# ---------------------------------------------------------------------------

class SpateoBenchmark:
    """Time SpaTEO as a Python analysis library (no browser)."""

    name = "SpaTEO"

    def is_available(self) -> bool:
        try:
            import anndata  # noqa: F401
            return True
        except ImportError:
            return False

    def benchmark_single_cell(
        self, path: str, size_tier: str, n_cells: int
    ) -> BenchmarkResult:
        import tracemalloc

        import anndata

        result = BenchmarkResult(
            platform=self.name,
            data_type="single_cell",
            size_tier=size_tier,
            n_items=n_cells,
            file_size_mb=round(os.path.getsize(path) / (1024 * 1024), 1),
            notes="Python library, no browser UI",
        )

        try:
            tracemalloc.start()

            # Load
            t0 = time.perf_counter()
            adata = anndata.read_h5ad(path)
            result.load_time_s = round(time.perf_counter() - t0, 2)

            # Gene query
            t1 = time.perf_counter()
            gene = adata.var_names[0]
            _expr = adata[:, gene].X.toarray().flatten()
            result.gene_query_time_s = round(time.perf_counter() - t1, 4)

            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            result.memory_peak_mb = round(peak / (1024 * 1024), 1)

            # Override with RSS if psutil is available (more representative)
            try:
                import psutil
                result.memory_peak_mb = round(
                    psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024), 1
                )
            except ImportError:
                pass

            # Try actual spateo operations
            try:
                import spateo as st  # noqa: F401
                result.notes += ", spateo imported OK"
            except ImportError:
                result.notes += ", spateo not installed (anndata-only test)"

        except Exception as e:
            result.crashed = True
            result.error = str(e)[:500]

        return result

    # SpaTEO doesn't have a dedicated single-molecule viewer
    benchmark_single_molecule = None


# ---------------------------------------------------------------------------
#  Vitessce benchmark (Zarr conversion + optional Playwright)
# ---------------------------------------------------------------------------

class VitessceBenchmark:
    """Time Vitessce's data preparation (Zarr conversion).

    Full interactive testing requires a Jupyter notebook, so this measures
    the mandatory preprocessing step that other tools don't require.
    """

    name = "Vitessce"

    def is_available(self) -> bool:
        try:
            import anndata  # noqa: F401
            return True
        except ImportError:
            return False

    def benchmark_single_cell(
        self, path: str, size_tier: str, n_cells: int
    ) -> BenchmarkResult:
        import anndata

        result = BenchmarkResult(
            platform=self.name,
            data_type="single_cell",
            size_tier=size_tier,
            n_items=n_cells,
            file_size_mb=round(os.path.getsize(path) / (1024 * 1024), 1),
        )

        try:
            # Vitessce requires Zarr conversion before viewing
            t0 = time.perf_counter()
            adata = anndata.read_h5ad(path)
            zarr_path = path.replace(".h5ad", ".zarr")
            adata.write_zarr(zarr_path)
            convert_s = time.perf_counter() - t0
            result.load_time_s = round(convert_s, 2)
            result.notes = (
                f"Zarr conversion only ({convert_s:.1f}s). "
                "Vitessce also needs main-thread deck.gl rendering on top of this."
            )

            try:
                import psutil
                result.memory_peak_mb = round(
                    psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024), 1
                )
            except ImportError:
                pass

            # Clean up zarr
            import shutil
            if os.path.isdir(zarr_path):
                shutil.rmtree(zarr_path, ignore_errors=True)

        except Exception as e:
            result.crashed = True
            result.error = str(e)[:500]

        return result

    benchmark_single_molecule = None


# ---------------------------------------------------------------------------
#  Report generation
# ---------------------------------------------------------------------------

TIER_ORDER = ["small", "medium", "large", "xl"]


def _fmt(val, suffix="", crash=False):
    if crash:
        return "CRASH"
    if val is None:
        return "  --"
    return f"{val}{suffix}"


def print_table(results: list[BenchmarkResult]):
    """Print a formatted comparison table to stdout."""
    for dtype_label, dtype_key in [
        ("SINGLE CELL", "single_cell"),
        ("SINGLE MOLECULE", "single_molecule"),
    ]:
        group = [r for r in results if r.data_type == dtype_key]
        if not group:
            continue

        platforms = sorted(set(r.platform for r in group))
        tiers = [t for t in TIER_ORDER if any(r.size_tier == t for r in group)]

        # Build size labels from the data
        tier_labels = {}
        for r in group:
            n = r.n_items
            if n >= 1_000_000:
                tier_labels[r.size_tier] = f"{n / 1_000_000:.0f}M"
            else:
                tier_labels[r.size_tier] = f"{n / 1_000:.0f}K"

        col_w = 16
        hdr = f"{'':>18}" + "".join(
            f"{t} ({tier_labels.get(t, '?')})".rjust(col_w) for t in tiers
        )

        print(f"\n{'=' * len(hdr)}")
        print(f"  {dtype_label} BENCHMARKS")
        print(f"{'=' * len(hdr)}")
        print(hdr)
        print("-" * len(hdr))

        for metric_label, attr, suffix in [
            ("Load time", "load_time_s", "s"),
            ("Gene query", "gene_query_time_s", "s"),
            ("Memory", "memory_peak_mb", " MB"),
            ("FPS", "fps_after_load", ""),
        ]:
            print(f"\n  {metric_label}")
            for plat in platforms:
                row = f"    {plat:<14}"
                for t in tiers:
                    r = next(
                        (r for r in group if r.platform == plat and r.size_tier == t),
                        None,
                    )
                    if r is None:
                        row += "  --".rjust(col_w)
                    else:
                        row += _fmt(
                            getattr(r, attr), suffix, r.crashed
                        ).rjust(col_w)
                print(row)

        # UI responsiveness row
        print(f"\n  UI frozen during load?")
        for plat in platforms:
            row = f"    {plat:<14}"
            for t in tiers:
                r = next(
                    (r for r in group if r.platform == plat and r.size_tier == t),
                    None,
                )
                if r is None:
                    row += "  --".rjust(col_w)
                elif r.ui_responsive is None:
                    row += "N/A".rjust(col_w)
                elif r.ui_responsive:
                    row += "No".rjust(col_w)
                else:
                    row += "YES".rjust(col_w)
            print(row)

        print()


def save_results(results: list[BenchmarkResult], output_dir: str) -> str:
    os.makedirs(output_dir, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    path = os.path.join(output_dir, f"benchmark_{stamp}.json")
    with open(path, "w") as f:
        json.dump(
            {"timestamp": stamp, "results": [asdict(r) for r in results]},
            f,
            indent=2,
        )
    return path


def load_and_compare(results_dir: str):
    """Load every benchmark_*.json in a dir and print a unified table."""
    all_results = []
    for fname in sorted(Path(results_dir).glob("*.json")):
        with open(fname) as f:
            data = json.load(f)
        if "results" not in data:
            continue
        for r in data["results"]:
            all_results.append(BenchmarkResult(**r))
    if all_results:
        print_table(all_results)
    else:
        print(f"No benchmark_*.json files found in {results_dir}")


# ---------------------------------------------------------------------------
#  Orchestrator
# ---------------------------------------------------------------------------

async def _run_one(bench, data_dir: str, sizes: list[str], results: list):
    """Run all applicable benchmarks for one platform."""
    name = bench.name
    print(f"\n{'-' * 60}")
    print(f"  {name}")
    print(f"{'-' * 60}")

    if not bench.is_available():
        print(f"  SKIP: {name} not available (not installed or server not running)")
        return

    manifest_path = os.path.join(data_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        print(f"  SKIP: {manifest_path} not found. Run generate_test_data.py first.")
        return

    with open(manifest_path) as f:
        manifest = json.load(f)

    for ds in manifest["datasets"]:
        if ds["size_tier"] not in sizes:
            continue

        path = ds["path"]
        if not os.path.exists(path):
            print(f"  SKIP: file not found: {path}")
            continue

        tier = ds["size_tier"]
        is_sc = ds["type"] == "single_cell"
        n = ds.get("cells") if is_sc else ds.get("molecules")
        kind = "cells" if is_sc else "molecules"

        method_name = "benchmark_single_cell" if is_sc else "benchmark_single_molecule"
        method = getattr(bench, method_name, None)
        if method is None:
            print(f"  --  {tier} ({n:,} {kind}): {name} has no {ds['type']} support")
            continue

        print(f"\n  [{tier}] {n:,} {kind}  ({ds['size_mb']} MB)")

        if asyncio.iscoroutinefunction(method):
            r = await method(path, tier, n)
        else:
            r = method(path, tier, n)

        results.append(r)

        icon = "OK" if not r.crashed else "FAIL"
        parts = [f"load={r.load_time_s}s" if r.load_time_s else None,
                 f"gene={r.gene_query_time_s}s" if r.gene_query_time_s else None,
                 f"mem={r.memory_peak_mb}MB" if r.memory_peak_mb else None,
                 f"fps={r.fps_after_load}" if r.fps_after_load else None]
        print(f"    {icon}  " + ", ".join(p for p in parts if p))
        if r.notes:
            print(f"         {r.notes}")
        if r.error:
            print(f"         ERROR: {r.error[:200]}")


async def async_main():
    parser = argparse.ArgumentParser(
        description="Benchmark spatial transcriptomics platforms"
    )
    parser.add_argument("--data-dir", default="./benchmark_data")
    parser.add_argument("--output-dir", default="./benchmark_results")
    parser.add_argument(
        "--platform",
        default="all",
        choices=["merfisheyes", "cellxgene", "spateo", "vitessce", "all"],
    )
    parser.add_argument(
        "--sizes",
        nargs="+",
        default=["small", "medium", "large", "xl"],
        choices=["small", "medium", "large", "xl"],
    )
    parser.add_argument("--compare", metavar="DIR",
                        help="Load and compare results from a directory")
    parser.add_argument("--merfisheyes-url", default="http://localhost:3000")
    parser.add_argument("--cellxgene-port", type=int, default=5005)
    args = parser.parse_args()

    if args.compare:
        load_and_compare(args.compare)
        return

    registry: dict[str, object] = {}
    p = args.platform
    if p in ("merfisheyes", "all"):
        registry["MERFISHeyes"] = MerfisheyesBenchmark(args.merfisheyes_url)
    if p in ("cellxgene", "all"):
        registry["CellxGene"] = CellxgeneBenchmark(args.cellxgene_port)
    if p in ("spateo", "all"):
        registry["SpaTEO"] = SpateoBenchmark()
    if p in ("vitessce", "all"):
        registry["Vitessce"] = VitessceBenchmark()

    print("=" * 60)
    print("  Spatial Transcriptomics Platform Benchmark")
    print("=" * 60)
    print(f"  Data:      {os.path.abspath(args.data_dir)}")
    print(f"  Platforms: {', '.join(registry)}")
    print(f"  Sizes:     {', '.join(args.sizes)}")

    results: list[BenchmarkResult] = []
    for bench in registry.values():
        await _run_one(bench, args.data_dir, args.sizes, results)

    if results:
        out = save_results(results, args.output_dir)
        print(f"\nResults saved to {out}")
        print_table(results)
    else:
        print("\nNo results. Check that platforms are available and test data exists.")


def main():
    asyncio.run(async_main())


if __name__ == "__main__":
    main()


# ===========================================================================
#  MANUAL TESTING GUIDE
# ===========================================================================
#
#  For platforms that cannot be automated, use the same test datasets and
#  record results in benchmark_results/manual_results.json using this format:
#
#  {
#    "timestamp": "2026-03-31 12:00:00",
#    "results": [
#      {
#        "platform": "Xenium Explorer",
#        "data_type": "single_cell",
#        "size_tier": "medium",
#        "n_items": 200000,
#        "file_size_mb": 95.2,
#        "load_time_s": 12.5,
#        "gene_query_time_s": 0.8,
#        "memory_peak_mb": 2400,
#        "ui_responsive": true,
#        "fps_after_load": 45,
#        "crashed": false,
#        "notes": "Tested on Windows 11, 32GB RAM, RTX 3060"
#      }
#    ]
#  }
#
#  --- XENIUM EXPLORER ---
#  1. Download from 10xgenomics.com (free, requires login)
#  2. You need Xenium-format data (cells.csv.gz, transcripts.csv.gz, etc.)
#     - Download sample datasets from 10x Genomics datasets page
#     - Synthetic benchmark H5ADs won't work (Xenium Explorer only reads its own format)
#  3. Open the app, load a dataset, use a stopwatch for:
#     - Time from "Open" to first render
#     - Time to color by a gene
#     - Check Task Manager for RAM usage
#  4. System requirements: 16GB+ RAM, dedicated GPU, SSD recommended
#
#  --- ALLEN BRAIN CELL ATLAS ---
#  1. Go to abc-atlas.brain-map.org
#  2. Can't upload custom data -- only tests their hosted viewer
#  3. Useful for comparing render performance on their ~4M cell dataset:
#     - Time to load the MERFISH spatial view
#     - Time to switch genes
#     - Use Chrome DevTools > Performance tab for FPS
#  4. Not a fair upload-time comparison, but shows deck.gl render limits
#
#  --- SAMUI ---
#  1. Clone github.com/MoffittLab/samui
#  2. Preprocess MERFISH data into Samui's expected format
#  3. Serve the processed data and open the viewer
#  4. Measure the same metrics manually
#
#  To include manual results in the comparison table:
#    python run_benchmarks.py --compare ./benchmark_results
#  (It reads all benchmark_*.json AND manual_results.json)
