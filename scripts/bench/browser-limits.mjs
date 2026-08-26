// Probe the browser's hard allocation limits — the ceilings that decide what
// the local upload path can and cannot open, independent of machine RAM.
//
//   node scripts/bench/browser-limits.mjs [--json bench/results/<env>/limits.json]
//
// Run this on every machine in the comparison. If a 16 GB laptop and a 500 GB
// workstation report the same numbers, that is the finding: the limit is
// per-renderer, so buying a bigger machine does not raise it.
//
// Two limits matter, and they fail differently:
//   * ArrayBuffer — File.arrayBuffer() must materialise the WHOLE file as one
//     buffer before h5wasm/hyparquet sees a byte. Over the cap the read is
//     refused instantly (RangeError), so a big file never starts processing.
//   * JS heap — row-oriented parsing (MERSCOPE / CSV → JS objects) grows the
//     V8 heap until the renderer is killed, which is a tab crash mid-parse.
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { chromium } from "@playwright/test";

const args = process.argv.slice(2);
const jsonAt = (() => {
  const i = args.indexOf("--json");
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : null;
})();

const browser = await chromium.launch({ args: ["--no-sandbox"] });
const page = await browser.newPage();
await page.goto("about:blank");

const probe = await page.evaluate(async () => {
  const out = { arrayBufferGB: {}, wasmMemoryGB: {} };

  // Bisect the largest single ArrayBuffer this renderer will hand out.
  let lo = 0, hi = 8192; // MB
  const canAlloc = (mb) => {
    try {
      const a = new ArrayBuffer(mb * 1024 * 1024);
      return a.byteLength === mb * 1024 * 1024;
    } catch {
      return false;
    }
  };
  while (hi - lo > 16) {
    const mid = Math.floor((lo + hi) / 2);
    if (canAlloc(mid)) lo = mid; else hi = mid;
  }
  out.maxArrayBufferGB = +(lo / 1024).toFixed(2);

  for (const gb of [0.5, 1, 1.5, 2, 3, 4]) out.arrayBufferGB[gb] = canAlloc(gb * 1024);
  for (const gb of [1, 2, 3, 4]) {
    try {
      new WebAssembly.Memory({ initial: Math.floor((gb * 1024 ** 3) / 65536) });
      out.wasmMemoryGB[gb] = true;
    } catch {
      out.wasmMemoryGB[gb] = false;
    }
  }

  out.jsHeapLimitGB = performance.memory
    ? +(performance.memory.jsHeapSizeLimit / 1024 ** 3).toFixed(2)
    : null;
  out.deviceMemoryGB = navigator.deviceMemory ?? null;
  out.userAgent = navigator.userAgent;
  return out;
});

const machine = {
  env: process.env.PERF_ENV || `${os.platform()}-${os.arch()}`,
  platform: `${os.platform()} ${os.release()}`,
  arch: os.arch(),
  cpu: os.cpus()[0]?.model ?? "unknown",
  cores: os.cpus().length,
  totalRamGB: +(os.totalmem() / 1024 ** 3).toFixed(2),
  chromium: browser.version(),
};
await browser.close();

const report = { machine, ...probe, measuredAt: new Date().toISOString() };
console.log(`${machine.cpu} · ${machine.cores} cores · ${machine.totalRamGB} GB RAM · Chromium ${machine.chromium}`);
console.log(`  max single ArrayBuffer : ${report.maxArrayBufferGB} GB   <- the largest FILE the browser can open`);
console.log(`  JS heap limit          : ${report.jsHeapLimitGB ?? "n/a"} GB   <- where CSV/row parsing crashes the tab`);
console.log(`  WebAssembly.Memory 4GB : ${report.wasmMemoryGB[4] ? "ok" : "refused"}`);
console.log(`  machine RAM            : ${machine.totalRamGB} GB (unused above the caps)`);

if (jsonAt) {
  fs.mkdirSync(path.dirname(jsonAt), { recursive: true });
  fs.writeFileSync(jsonAt, JSON.stringify(report, null, 1) + "\n");
  console.log(`\nwrote ${jsonAt}`);
}
