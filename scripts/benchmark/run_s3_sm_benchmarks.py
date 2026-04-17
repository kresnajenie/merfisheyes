#!/usr/bin/env python3
"""
Benchmark single molecule S3 lazy-loading via the "Load from S3" feature.

Navigates to the homepage, enters each S3 URL in the Load from S3 modal,
and measures manifest download + initial gene loading + rendering time.

Usage:
    python run_s3_sm_benchmarks.py
    python run_s3_sm_benchmarks.py --headless
"""

import argparse
import asyncio
import csv
import json
import os
import time

BASE_URL = "http://localhost:3000"

S3_DATASETS = [
    {
        "name": "infancyVC_SM",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/infancyVC_SM/",
        "molecules": 11_526_235,
        "genes": 300,
    },
    {
        "name": "ace-dud-vex-sm (monkey)",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-vex-sm/",
        "molecules": 20_678_792,
        "genes": 300,
    },
    {
        "name": "ace-dry-dry-sm",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/ace-dry-dry-sm/",
        "molecules": 29_238_426,
        "genes": 300,
    },
    {
        "name": "ace-dud-wag-sm (monkey)",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-wag/ace-dud-wag-sm/",
        "molecules": 53_242_389,
        "genes": 300,
    },
    {
        "name": "ace-dip-dog-sm",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/ace-dip-dog-sm/",
        "molecules": 67_494_354,
        "genes": 300,
    },
    {
        "name": "ace-dud-was-sm (monkey)",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/monkeys/ace-dud-was/ace-dud-was-sm/",
        "molecules": 102_297_212,
        "genes": 300,
    },
    {
        "name": "12_04_2025_sm (combined)",
        "url": "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/12_04_2025_sm/",
        "molecules": 1_223_821_804,
        "genes": 500,
    },
]

OUTPUT_CSV = "benchmark_results/s3_sm_benchmark.csv"


async def measure_memory(cdp):
    try:
        metrics = await cdp.send("Performance.getMetrics")
        for m in metrics.get("metrics", []):
            if m["name"] == "JSHeapUsedSize":
                return round(m["value"] / (1024 * 1024), 1)
    except Exception:
        pass
    return None


async def measure_fps(page):
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


async def benchmark_s3_dataset(ds, base_url, headless):
    from playwright.async_api import async_playwright

    result = {
        "name": ds["name"],
        "url": ds["url"],
        "molecules": ds["molecules"],
        "genes": ds["genes"],
        "load_time_s": None,
        "memory_peak_mb": None,
        "fps_after_load": None,
        "crashed": False,
        "error": None,
    }

    async with async_playwright() as pw:
        browser = await pw.chromium.launch(
            headless=headless,
            args=["--no-sandbox", "--disable-dev-shm-usage"],
        )
        page = await browser.new_page()
        cdp = await page.context.new_cdp_session(page)
        await cdp.send("Performance.enable")

        try:
            # Navigate to SM homepage
            await page.goto(f"{base_url}?mode=sm", wait_until="networkidle", timeout=30_000)
            await page.wait_for_timeout(2000)

            # Click "Load from S3" — an animated overlay can intercept clicks,
            # so use JavaScript to trigger the modal directly
            await page.evaluate("""() => {
                // Find all buttons with "Load from S3" text and click the visible one
                const buttons = document.querySelectorAll('button');
                for (const btn of buttons) {
                    if (btn.textContent.includes('Load from S3') && btn.offsetParent !== null) {
                        btn.click();
                        return true;
                    }
                }
                return false;
            }""")
            await page.wait_for_timeout(1000)

            # Fill the URL input in the modal
            url_input = page.get_by_label("Dataset Folder URL")
            await url_input.fill(ds["url"])

            # Start timing from the Load button click
            t0 = time.perf_counter()

            load_btn = page.get_by_role("button", name="Load Dataset")
            await load_btn.click()

            # Wait for navigation to sm-viewer/from-s3
            timeout_ms = 300_000  # 5 min max
            await page.wait_for_url("**/sm-viewer/from-s3**", timeout=timeout_ms)
            await page.wait_for_selector("canvas", timeout=120_000)
            await page.wait_for_timeout(5000)  # let rendering stabilize

            result["load_time_s"] = round(time.perf_counter() - t0, 2)
            result["memory_peak_mb"] = await measure_memory(cdp)
            result["fps_after_load"] = await measure_fps(page)

        except Exception as e:
            result["crashed"] = True
            result["error"] = str(e)[:300]
            elapsed = time.perf_counter() - t0 if 't0' in dir() else 0
            result["load_time_s"] = round(elapsed, 2) if elapsed else None

        finally:
            await browser.close()

    return result


def fmt_mol(n):
    if n >= 1_000_000_000:
        return f"{n/1e9:.1f}B"
    if n >= 1_000_000:
        return f"{n/1e6:.1f}M"
    return f"{n/1e3:.0f}K"


async def main():
    parser = argparse.ArgumentParser(description="Benchmark SM S3 lazy-loading")
    parser.add_argument("--headless", action="store_true")
    parser.add_argument("--url", default=BASE_URL)
    args = parser.parse_args()

    print("=" * 70)
    print("  MERFISHeyes S3 Single Molecule Lazy-Loading Benchmark")
    print("=" * 70)
    print(f"  Server: {args.url}")
    print(f"  Datasets: {len(S3_DATASETS)}")
    print()

    results = []

    for i, ds in enumerate(S3_DATASETS):
        print(f"  [{i+1}/{len(S3_DATASETS)}] {ds['name']} ({fmt_mol(ds['molecules'])} molecules, {ds['genes']} genes)")

        r = await benchmark_s3_dataset(ds, args.url, args.headless)
        results.append(r)

        if r["crashed"]:
            print(f"    CRASHED: {r['error'][:100]}")
        else:
            load = f"{r['load_time_s']:.1f}s" if r['load_time_s'] else "--"
            mem = f"{r['memory_peak_mb']:.0f}MB" if r['memory_peak_mb'] else "--"
            fps = f"{r['fps_after_load']}" if r['fps_after_load'] else "--"
            print(f"    load={load}  mem={mem}  fps={fps}")

    # Save results
    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "name", "url", "molecules", "genes",
            "load_time_s", "memory_peak_mb", "fps_after_load",
            "crashed", "error",
        ])
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\nResults saved to {OUTPUT_CSV}")

    # Print summary table
    print("\n" + "=" * 80)
    print(f"  {'Dataset':<30} {'Molecules':>12} {'Load':>8} {'Memory':>8} {'FPS':>6}")
    print("-" * 80)
    for r in results:
        name = r["name"][:30]
        mol = fmt_mol(r["molecules"])
        load = f"{r['load_time_s']:.1f}s" if r["load_time_s"] and not r["crashed"] else "CRASH"
        mem = f"{r['memory_peak_mb']:.0f}MB" if r["memory_peak_mb"] else "--"
        fps = str(r["fps_after_load"]) if r["fps_after_load"] else "--"
        print(f"  {name:<30} {mol:>12} {load:>8} {mem:>8} {fps:>6}")
    print("=" * 80)


if __name__ == "__main__":
    asyncio.run(main())
