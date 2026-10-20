/**
 * The budgets (docs/specs/BUDGETS.md) measured on the REAL stores, through
 * headless Chromium against a local dev server started with the E2E hooks:
 *
 *   NEXT_PUBLIC_E2E=1 npm run dev -- --port 3002        # in one terminal
 *   node scripts/perf/cors-proxy.cjs 3977               # in another
 *   npx tsx scripts/perf/real-store.ts [--app http://localhost:3002] [--proxy http://127.0.0.1:3977] [--store all_sets|atera|skin] [--gene Slc17a7]
 *
 * The proxy is needed because the bucket's CORS only allows the standard
 * dev ports, and it logs bytes per request, which the browser cannot see
 * for a cross-origin bucket. Results go to perf/real/<store>.json (committed, unlike perf/results/)
 * (median of the samples, plus every sample), so the "today" column of
 * BUDGETS.md and the before/after of each step come from the same script.
 *
 * Not part of the CI gate: it depends on the network and on the bucket.
 */
import fs from "node:fs";
import path from "node:path";

import { chromium, type Page } from "playwright";

const arg = (flag: string, fallback: string) => {
  const i = process.argv.indexOf(flag);

  return i >= 0 ? process.argv[i + 1] : fallback;
};

const APP = arg("--app", "http://localhost:3002");
const PROXY = arg("--proxy", "http://127.0.0.1:3977");
const WHICH = arg("--store", "all_sets");
const GENE = arg("--gene", "");
const SAMPLES = Number(arg("--samples", "3"));

const STORES: Record<string, { folder: string; gene: string }> = {
  all_sets: { folder: "spatialdata/12_16_2025_all_sets/", gene: "Slc17a7" },
  atera: { folder: "spatialdata/atera_breast_cancer/", gene: "EPCAM" },
  skin: { folder: "spatialdata/skin_adapt_v2/", gene: "EEF1G" },
  set1: { folder: "spatialdata/12_16_2025_sets/set1/", gene: "Slc17a7" },
};

interface Hooks {
  layers: Record<string, { firstContentAt: number | null; settledAt: number | null; changedAt: number | null; drawn: number; wanted: number; resident: number; bytes: number; requests: number }>;
  frames: { count: number; over16ms: number; over33ms: number; worstMs: number };
  memory: { jsHeap: number | null; drawCalls: number };
  renderComplete: number | null;
}

const hooks = (page: Page) => page.evaluate(() => (window as unknown as { __merfish?: Hooks }).__merfish ?? null);
const now = (page: Page) => page.evaluate(() => performance.now());
const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

async function until<T>(
  what: string,
  fn: () => Promise<T | null | undefined>,
  timeout = 600_000,
  every = 100,
): Promise<T> {
  const t0 = Date.now();

  for (;;) {
    const v = await fn();

    if (v !== null && v !== undefined) return v;
    if (Date.now() - t0 > timeout) throw new Error(`timed out after ${timeout / 1000} s waiting for ${what}`);
    await sleep(every);
  }
}

async function proxyLog(): Promise<{ path: string; bytes: number; ms: number }[]> {
  return (await fetch(`${PROXY}/__log`)).json();
}
async function proxyReset(): Promise<void> {
  await fetch(`${PROXY}/__reset`);
}
const bytesOf = (log: { path: string; bytes: number }[], needle: string) =>
  log.filter((r) => r.path.includes(needle)).reduce((a, r) => a + r.bytes, 0);

async function sample(store: { folder: string; gene: string }) {
  const browser = await chromium.launch({ args: ["--use-gl=swiftshader"] });
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  const out: Record<string, number> = {};

  try {
    await proxyReset();
    const url = `${APP}/spatialdata-viewer/from-s3?url=${encodeURIComponent(`${PROXY}/${store.folder}`)}`;

    await page.goto(url, { waitUntil: "commit" });
    const t0 = await now(page);

    const cellsAt = await until("cells", async () => (await hooks(page))?.layers.table?.firstContentAt);

    out.openToCells = cellsAt - t0;
    const log1 = await proxyLog();

    out.openBytesKB = log1.reduce((a, r) => a + r.bytes, 0) / 1024;
    out.openRequests = log1.length;

    // Molecules: the default genes, or one named gene selected by hand
    if (GENE || store.gene) {
      const gene = GENE || store.gene;
      const toggle = page.locator('[aria-label="Single molecule"]');

      await toggle.first().waitFor({ timeout: 120_000 });
      await toggle.first().click();
      await page.locator('button:has-text("Gene")').first().click();
      await page.locator('input[placeholder="Search genes..."]').fill(gene);
      await page.waitForTimeout(300);
      await proxyReset();
      const selectAt = await now(page);

      await page.locator("label", { hasText: new RegExp(`^${gene}`) }).first().click();
      const first = await until("first molecules of the gene", async () => {
        const l = (await hooks(page))?.layers.points;

        return l?.firstContentAt && l.firstContentAt > selectAt ? l.firstContentAt : null;
      });

      out.geneToFirstMolecules = first - selectAt;
      const settled = await until("the gene's settled view", async () => {
        const l = (await hooks(page))?.layers.points;

        return l?.settledAt && l.settledAt > selectAt ? l.settledAt : null;
      });

      out.geneToSettled = settled - selectAt;
      out.geneBytesKB = bytesOf(await proxyLog(), "/points/") / 1024;
    }

    // A zoom in by three notches
    const canvas = page.locator('[data-testid="sc-scene-canvas"] canvas').first();
    const box = (await canvas.boundingBox())!;

    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    await proxyReset();
    const zoomAt = await now(page);

    // 64 notches (0.95 each) in ONE event, so one replan: a 40 mm atlas
    // view becomes about 1.5 mm. (A stepped zoom queues every intermediate
    // view's reads, which today are never cancelled.)
    await page.mouse.wheel(0, -6400);
    const changed = await until("a replan after the zoom", async () => {
      const l = (await hooks(page))?.layers.points;

      return l?.changedAt && l.changedAt > zoomAt ? l.changedAt : null;
    }, 30_000);
    const settledZoom = await until("the zoomed view to settle", async () => {
      const l = (await hooks(page))?.layers.points;

      return l?.settledAt && l.settledAt >= changed ? l.settledAt : null;
    });

    out.viewSettled = settledZoom - changed;
    const zlog = await proxyLog();

    out.bytesPerViewKB = bytesOf(zlog, "/points/") / 1024;
    out.requestsPerView = zlog.filter((r) => r.path.includes("/points/")).length;

    // A one-second pan: frames
    await page.evaluate(() => {
      const h = (window as unknown as { __merfish?: Hooks }).__merfish;

      if (h) h.frames = { count: 0, drawn: 0, over16ms: 0, over33ms: 0, worstMs: 0 } as Hooks["frames"];
    });
    const x0 = box.x + box.width * 0.35;
    const y = box.y + box.height / 2;

    await page.mouse.move(x0, y);
    await page.mouse.down();
    for (let i = 1; i <= 30; i++) {
      await page.mouse.move(x0 + (box.width * 0.3 * i) / 30, y);
      await sleep(33);
    }
    await page.mouse.up();
    const f = (await hooks(page))!.frames;

    out.frameOver16Pct = f.count ? (100 * f.over16ms) / f.count : 0;
    out.worstFrameMs = f.worstMs;
    const mem = (await hooks(page))!.memory;

    if (mem.jsHeap) out.jsHeapMB = mem.jsHeap / 1048576;
    out.drawCalls = mem.drawCalls;
    const pts = (await hooks(page))!.layers.points;

    if (pts) {
      out.moleculesDrawn = pts.drawn;
      out.moleculesResident = pts.resident;
    }
  } finally {
    await browser.close();
  }

  return out;
}

(async () => {
  const store = STORES[WHICH];

  if (!store) throw new Error(`unknown store ${WHICH}; one of ${Object.keys(STORES).join(", ")}`);
  const samples: Record<string, number>[] = [];

  for (let i = 0; i < SAMPLES; i++) {
    const s = await sample(store);

    samples.push(s);
    console.log(`sample ${i + 1}:`, JSON.stringify(s));
  }
  const keys = [...new Set(samples.flatMap(Object.keys))];
  const median: Record<string, number> = {};

  for (const k of keys) {
    const v = samples.map((s) => s[k]).filter((x) => x !== undefined).sort((a, b) => a - b);

    median[k] = v[Math.floor(v.length / 2)];
  }
  const dir = path.resolve(__dirname, "..", "..", "perf", "real");

  fs.mkdirSync(dir, { recursive: true });
  const file = path.join(dir, `${WHICH}.json`);

  fs.writeFileSync(
    file,
    JSON.stringify({ store: WHICH, folder: store.folder, app: APP, measuredAt: new Date().toISOString(), median, samples }, null, 2),
  );
  console.log("median:", JSON.stringify(median, null, 1));
  console.log("wrote", path.relative(process.cwd(), file));
})().catch((e) => {
  console.error(e);
  process.exit(1);
});
