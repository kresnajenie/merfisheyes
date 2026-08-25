// Benchmark the LOCAL (in-browser) path: hand a dataset to the app in a real
// Chromium tab, time the processing, sample memory, and record where it dies.
//
//   node scripts/bench/local-run.mjs --datasets bench/datasets.json \
//        [--base http://localhost:3100] [--only <id> ...] [--format h5ad ...] [--max-gb 5]
//        [--timeout 1800] [--heap-mb N] [--headed]
//
// The browser tab crashing IS the measurement — that's the ceiling the server
// path exists to remove. Timing excludes byte transfer (nothing is uploaded
// here): what's measured is parse + build + first render, the browser's
// equivalent of the worker's process_ms.
//
// Needs a server running the app built with NEXT_PUBLIC_E2E=1 (the run records
// whether the hooks were present; without them timing falls back to wall clock).
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { execSync } from "node:child_process";
import { chromium } from "@playwright/test";
import { classifyDir, stageDir } from "./classify.mjs";

const args = process.argv.slice(2);
const flag = (name, fallback = null) => {
  const i = args.indexOf(`--${name}`);
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : fallback;
};
const has = (name) => args.includes(`--${name}`);
const multi = (name) => {
  const out = [];
  for (let i = 0; i < args.length; i++) {
    if (args[i] === `--${name}`) {
      for (let j = i + 1; j < args.length && !args[j].startsWith("--"); j++) out.push(args[j]);
    }
  }
  return out;
};

const DATASETS = flag("datasets", "bench/datasets.json");
const BASE = flag("base", "http://localhost:3100");
const ONLY = multi("only");
const FORMATS = multi("format");
const MAX_GB = Number(flag("max-gb", "0")) || Infinity;
const TIMEOUT_MS = Number(flag("timeout", "1800")) * 1000;
const HEAP_MB = flag("heap-mb");
const ENV_NAME = process.env.PERF_ENV || `${os.platform()}-${os.arch()}`;
const OUT_DIR = path.join("bench", "results", ENV_NAME);

const FILE_INPUT = {
  h5ad: "#sc-file-input-h5ad",
  xenium: "#sc-file-input-folder",
  merscope: "#sc-file-input-folder",
  "single-molecule": "#sm-file-input-file",
};
const CANVAS = {
  "single-molecule": '[data-testid="sm-scene-canvas"] canvas',
  default: '[data-testid="sc-scene-canvas"] canvas',
};

/** Chromium processes visible right now, as {pid, ppid, isBrowser}. */
function chromiumProcs() {
  const out = [];
  let entries;
  try {
    entries = fs.readdirSync("/proc").filter((n) => /^\d+$/.test(n));
  } catch {
    return out; // not Linux
  }
  for (const pid of entries) {
    try {
      const cmd = fs.readFileSync(`/proc/${pid}/cmdline`, "utf8");
      if (!cmd.includes("ms-playwright")) continue;
      const stat = fs.readFileSync(`/proc/${pid}/stat`, "utf8");
      const ppid = Number(stat.slice(stat.lastIndexOf(")") + 2).split(" ")[1]);
      out.push({ pid: Number(pid), ppid, userDataDir: /user-data-dir=([^\0 ]+)/.exec(cmd)?.[1] ?? null });
    } catch {
      /* process gone */
    }
  }
  return out;
}

/**
 * Total RSS of ONE browser's process tree, in bytes.
 *
 * Scoped to the tree rooted at `rootPid`, not to every Playwright process:
 * two benchmark runs at once would otherwise sum into each other's peaks (a
 * concurrent run once reported an h5ad peak identical to the other ladder's
 * single-molecule peak). Chromium's zygote/renderer/gpu children are all
 * descendants of the browser process, so the walk is exact — and it captures
 * the h5wasm heap, which the JS heap number alone would miss.
 */
function treeRss(rootPid) {
  if (!rootPid) return 0;
  const procs = chromiumProcs();
  const byParent = new Map();
  for (const p of procs) {
    if (!byParent.has(p.ppid)) byParent.set(p.ppid, []);
    byParent.get(p.ppid).push(p.pid);
  }
  const tree = [];
  const stack = [rootPid];
  while (stack.length) {
    const pid = stack.pop();
    tree.push(pid);
    for (const kid of byParent.get(pid) ?? []) stack.push(kid);
  }
  let total = 0;
  for (const pid of tree) {
    try {
      const m = fs.readFileSync(`/proc/${pid}/status`, "utf8").match(/VmRSS:\s+(\d+) kB/);
      if (m) total += Number(m[1]) * 1024;
    } catch {
      /* gone */
    }
  }
  return total;
}

const gb = (b) => +(b / 1024 ** 3).toFixed(2);

async function runOne(dataset, browser, rootPid) {
  const fmt = dataset.declaredFormat;
  const source = path.join(dataset.root, dataset.datasetRoot ?? dataset.path);
  const isSm = fmt === "single-molecule";
  // Folder formats: hand the input a staged directory holding only the files
  // the app would consume (see classify.mjs) — Playwright would otherwise
  // stream the entire export over CDP, which a real browser never does.
  let staged = null;
  let abs = source;
  let uploadFiles = null;
  if (fmt === "xenium" || fmt === "merscope") {
    const picked = classifyDir(source, fmt);
    // Stage beside the repo, not in /tmp: same filesystem as the databank, so
    // the files are hardlinked (instant, no extra disk) instead of copied.
    staged = path.join(process.cwd(), ".bench-stage", dataset.id);
    abs = stageDir(picked, staged);
    uploadFiles = { count: picked.length, gb: gb(picked.reduce((a, f) => a + f.size, 0)) };
  }
  const context = await browser.newContext({ viewport: { width: 1280, height: 800 } });
  const page = await context.newPage();

  let crashed = false;
  page.on("crash", () => { crashed = true; });
  const consoleErrors = [];
  page.on("console", (m) => { if (m.type() === "error") consoleErrors.push(m.text().slice(0, 200)); });

  // Sample the whole Chromium tree: JS heap alone misses the h5wasm heap.
  let peakRss = 0;
  const sampler = setInterval(() => { peakRss = Math.max(peakRss, treeRss(rootPid)); }, 500);

  const started = Date.now();
  const result = {
    datasetId: dataset.id,
    format: fmt,
    path: dataset.path,
    sizeGB: dataset.sizeGB,
    relevantGB: dataset.relevantGB ?? null,
    shape: dataset.shape,
    outcome: "ok",
    processMs: null,
    sceneMs: null,
    peakRssGB: null,
    jsHeapGB: null,
    hooks: false,
    error: null,
  };

  try {
    await page.goto(`${BASE}${isSm ? "/?mode=sm" : "/"}`, {
      waitUntil: "domcontentloaded", timeout: 60_000,
    });
    const input = page.locator(FILE_INPUT[fmt]);
    await input.waitFor({ state: "attached", timeout: 30_000 });

    const t0 = await page.evaluate(() => performance.now());
    // Explicit timeout: a folder upload streams every file's contents over CDP
    // and blows past Playwright's 30 s default long before processing starts.
    await input.setInputFiles(abs, { timeout: TIMEOUT_MS });

    if (isSm) {
      // Both the wait AND the click need the explicit timeout: the confirm
      // dialog re-renders while the file is still being read, so Playwright's
      // 30 s default can expire mid-actionability-check on a large file.
      const confirm = page.getByRole("button", { name: /Looks good/i });
      await confirm.waitFor({ state: "visible", timeout: TIMEOUT_MS });
      await confirm.click({ timeout: TIMEOUT_MS });
    }

    // Race the viewer against the app's own error toast: a dataset the browser
    // cannot allocate fails in seconds, and without this the run would sit out
    // the whole timeout before recording an outcome it already knew.
    const errorToast = page.locator(".Toastify__toast--error");
    const navigated = page.waitForURL(isSm ? "**/sm-viewer**" : "**/viewer**", { timeout: TIMEOUT_MS })
      .then(() => "viewer");
    const failed = errorToast.first().waitFor({ state: "visible", timeout: TIMEOUT_MS })
      .then(() => "toast");
    const first = await Promise.race([navigated, failed]).catch((e) => { throw e; });
    if (first === "toast") {
      const text = (await errorToast.first().innerText().catch(() => "")).replace(/\s+/g, " ").trim();
      throw new Error(`app reported: ${text.slice(0, 200)}`);
    }
    await page.locator(CANVAS[isSm ? "single-molecule" : "default"]).first()
      .waitFor({ state: "visible", timeout: TIMEOUT_MS });

    const hooks = await page.evaluate(() => window.__merfish ?? null);
    if (hooks) {
      result.hooks = true;
      await page.waitForFunction(() => window.__merfish?.renderComplete != null, undefined, { timeout: TIMEOUT_MS });
      const h = await page.evaluate(() => window.__merfish);
      result.processMs = Math.round((h.renderComplete ?? 0) - t0);
      result.sceneMs = h.sceneReady != null ? Math.round(h.sceneReady - t0) : null;
      result.stats = h.stats ?? null;
    } else {
      const end = await page.evaluate(() => performance.now());
      result.processMs = Math.round(end - t0);
    }

    try {
      const cdp = await context.newCDPSession(page);
      const { metrics } = await cdp.send("Performance.getMetrics");
      const used = metrics.find((m) => m.name === "JSHeapUsedSize");
      if (used) result.jsHeapGB = gb(used.value);
    } catch { /* tab may be gone */ }
  } catch (e) {
    const msg = String(e?.message ?? e);
    result.outcome = crashed
      ? "oom"
      : /^app reported:/.test(msg) ? "failed"
      : /Timeout|exceeded/i.test(msg) ? "timeout"
      : /Target (page|browser).*closed|Target closed/i.test(msg) ? "aborted"
      : "error";
    result.error = msg.split("\n")[0].slice(0, 300);
  } finally {
    clearInterval(sampler);
    peakRss = Math.max(peakRss, treeRss(rootPid));
    result.peakRssGB = gb(peakRss);
    result.wallMs = Date.now() - started;
    if (uploadFiles) { result.inputFiles = uploadFiles.count; result.inputGB = uploadFiles.gb; }
    if (staged) fs.rmSync(staged, { recursive: true, force: true });
    if (consoleErrors.length) result.consoleErrors = consoleErrors.slice(0, 5);
    await context.close().catch(() => {});
  }
  return result;
}

async function main() {
  const db = JSON.parse(fs.readFileSync(DATASETS, "utf8"));
  let list = db.datasets.map((d) => ({ ...d, root: db.root }));
  if (ONLY.length) list = list.filter((d) => ONLY.includes(d.id));
  if (FORMATS.length) list = list.filter((d) => FORMATS.includes(d.declaredFormat));
  // A dataset with no cells and no molecules is one the processor would reject
  // (e.g. a raw imaging folder). Feeding a 60 GB tree of images to the file
  // input just makes the browser enumerate it for minutes; skip it up front.
  const usable = (d) => (d.shape?.cells ?? 0) > 0 || (d.shape?.molecules ?? 0) > 0;
  const skipped = list.filter((d) => !usable(d));
  list = list.filter(usable);
  for (const d of skipped) console.log(`skipped ${d.id} — no cells/molecules detected`);

  // Rung size: what the pipeline reads, but never 0 (that would sort an
  // unusable folder to the front of the ladder).
  const rung = (d) => (d.relevantGB ? d.relevantGB : d.sizeGB);
  list = list.filter((d) => FILE_INPUT[d.declaredFormat] && rung(d) <= MAX_GB);
  // Smallest first: an ascending ladder finds the ceiling without wasting a
  // long run on a dataset that was never going to fit.
  list.sort((a, b) => rung(a) - rung(b));

  if (!list.length) {
    console.error("No datasets selected."); process.exit(2);
  }

  const launchArgs = ["--no-sandbox"];
  if (HEAP_MB) launchArgs.push(`--js-flags=--max-old-space-size=${HEAP_MB}`);
  // Note which browsers already exist so this run's own tree can be identified
  // even when another benchmark is running.
  const before = new Set(chromiumProcs().map((p) => p.pid));
  const browser = await chromium.launch({ headless: !has("headed"), args: launchArgs });
  const rootPid = chromiumProcs().find((p) => p.userDataDir && !before.has(p.pid))?.pid ?? null;
  if (!rootPid) console.log("warning: could not identify this run's browser process; peak RSS will read 0");

  const machine = {
    env: ENV_NAME,
    platform: `${os.platform()} ${os.release()}`,
    arch: os.arch(),
    cpu: os.cpus()[0]?.model ?? "unknown",
    cores: os.cpus().length,
    totalRamGB: gb(os.totalmem()),
    node: process.version,
    chromium: browser.version(),
    base: BASE,
    heapCapMB: HEAP_MB ? Number(HEAP_MB) : null,
  };
  console.log(`machine: ${machine.cpu} · ${machine.cores} cores · ${machine.totalRamGB} GB RAM · Chromium ${machine.chromium}`);
  console.log(`${"dataset".padEnd(52)} ${"GB".padStart(6)} ${"outcome".padEnd(8)} ${"process".padStart(9)} ${"peakRSS".padStart(8)}`);
  console.log("-".repeat(92));

  fs.mkdirSync(OUT_DIR, { recursive: true });
  const results = [];
  for (const d of list) {
    const r = await runOne(d, browser, rootPid);
    results.push(r);
    const secs = r.processMs != null ? `${(r.processMs / 1000).toFixed(1)}s` : "-";
    console.log(
      `${r.datasetId.slice(0, 52).padEnd(52)} ${String(r.relevantGB || r.sizeGB).padStart(6)} ` +
      `${r.outcome.padEnd(8)} ${secs.padStart(9)} ${(r.peakRssGB + " GB").padStart(8)}` +
      (r.error ? `\n${" ".repeat(52)}   ! ${r.error}` : ""),
    );
    fs.writeFileSync(
      path.join(OUT_DIR, `local-${r.datasetId}.json`),
      JSON.stringify({ machine, path: "local", ...r }, null, 1) + "\n",
    );
  }
  await browser.close();

  // A single-cell run that loaded zero genes is degraded (Xenium h5/tar.gz
  // expression is unreadable in the browser) — it must not count as the
  // largest local success.
  const degraded = (r) => r.stats && r.stats.dataType === "single_cell" && (r.stats.exprGeneCount ?? r.stats.geneCount) === 0;
  const ok = results.filter((r) => r.outcome === "ok" && !degraded(r));
  const failed = results.filter((r) => r.outcome !== "ok");
  console.log("-".repeat(92));
  console.log(`${ok.length}/${results.length} succeeded → ${OUT_DIR}/`);
  if (ok.length) {
    const biggest = ok.reduce((a, b) => ((a.relevantGB || a.sizeGB) > (b.relevantGB || b.sizeGB) ? a : b));
    console.log(`largest processed locally: ${biggest.datasetId} (${biggest.relevantGB || biggest.sizeGB} GB, peak ${biggest.peakRssGB} GB)`);
  }
  if (failed.length) {
    console.log(`ceiling hit at: ${failed.map((f) => `${f.datasetId} [${f.outcome}]`).join(", ")}`);
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
