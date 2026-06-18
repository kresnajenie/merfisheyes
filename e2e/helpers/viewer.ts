import fs from "node:fs";
import { type Page, expect } from "@playwright/test";

import { type DatasetEntry, absPath } from "./databank";
import { SELECTORS, fileInputSelector } from "./selectors";

/** Shape of the window.__merfish test hook (mirror of lib/utils/test-hooks.ts). */
export interface MerfishHooks {
  enabled: boolean;
  sceneReady: number | null;
  renderComplete: number | null;
  renderCount: number;
  datasetLoadedAt: number | null;
  stats: {
    pointCount?: number;
    geneCount?: number;
    dimensions?: number;
    dataType?: string;
    pointClouds?: number;
  } | null;
  lastGene: { gene: string; count: number; at: number } | null;
}

export interface LoadTiming {
  /** Scene initialized (≈ time to initial render). */
  loadToScene: number | null;
  /** Point cloud attributes applied (≈ time to interactive). */
  loadToRender: number;
}

// Load/render budget. Large in-browser datasets (e.g. a 500k-cell h5ad ~1GB)
// can take minutes to parse, so this is overridable via E2E_TIMEOUT_MS.
const HOOK_TIMEOUT = Number(process.env.E2E_TIMEOUT_MS) || 120_000;

export async function readHooks(page: Page): Promise<MerfishHooks | null> {
  return page.evaluate(() => (window as unknown as { __merfish?: MerfishHooks }).__merfish ?? null);
}

async function hooksEnabled(page: Page): Promise<boolean> {
  const h = await readHooks(page);
  return !!h?.enabled;
}

function canvasSelector(ds: DatasetEntry): string {
  return ds.dataType === "single_molecule" ? SELECTORS.smCanvas : SELECTORS.scCanvas;
}

function viewerUrlGlob(ds: DatasetEntry): string {
  return ds.dataType === "single_molecule"
    ? "**/sm-viewer/from-local"
    : "**/viewer/from-local";
}

/**
 * Load a local dataset through the homepage upload entry, then wait for the
 * viewer to finish its initial render. Returns timing relative to the moment
 * the file was handed to the app (page monotonic clock).
 */
export async function loadLocalDataset(
  page: Page,
  ds: DatasetEntry,
): Promise<LoadTiming> {
  const home = ds.dataType === "single_molecule" ? "/?mode=sm" : "/";
  await page.goto(home, { waitUntil: "domcontentloaded" });

  const input = page.locator(fileInputSelector(ds.format));
  await expect(input).toBeAttached();

  // Mark the action start on the page's monotonic clock.
  const start = await page.evaluate(() => performance.now());
  await input.setInputFiles(absPath(ds)); // directory path for folder formats

  await page.waitForURL(viewerUrlGlob(ds), { timeout: HOOK_TIMEOUT });
  await waitForRender(page, ds);

  if (await hooksEnabled(page)) {
    const h = (await readHooks(page))!;
    return {
      loadToScene: h.sceneReady != null ? h.sceneReady - start : null,
      loadToRender: (h.renderComplete ?? (await page.evaluate(() => performance.now()))) - start,
    };
  }
  // Fallback (hooks disabled): wall-clock-ish via canvas appearance.
  const end = await page.evaluate(() => performance.now());
  return { loadToScene: null, loadToRender: end - start };
}

/** Wait until the scene has rendered (hook-based, with canvas fallback). */
export async function waitForRender(page: Page, ds: DatasetEntry): Promise<void> {
  await expect(page.locator(canvasSelector(ds)).first()).toBeVisible({ timeout: HOOK_TIMEOUT });
  if (await hooksEnabled(page)) {
    await page.waitForFunction(
      () => {
        const m = (window as unknown as { __merfish?: MerfishHooks }).__merfish;
        return !!m && m.renderComplete != null;
      },
      undefined,
      { timeout: HOOK_TIMEOUT },
    );
  } else {
    // No precise signal — give the WebGL pipeline a brief settle.
    await page.waitForTimeout(2000);
  }
}

/** Select a gene and measure the latency from click to render-complete (ms). */
export async function selectGeneAndMeasure(
  page: Page,
  ds: DatasetEntry,
  gene: string,
): Promise<number> {
  return ds.dataType === "single_molecule"
    ? selectGeneSM(page, gene)
    : selectGeneSC(page, gene);
}

async function selectGeneSC(page: Page, gene: string): Promise<number> {
  // The mode button toggles the panel, so only click it when the panel is
  // closed (otherwise a second gene selection would close it).
  const search = page.locator(SELECTORS.geneSearchInput);
  if (!(await search.isVisible())) {
    await page.locator(SELECTORS.geneModeButton).first().click();
  }
  await expect(search).toBeVisible();
  await search.fill(gene);

  const before = await page.evaluate(
    () => (window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderCount ?? 0,
  );
  const start = await page.evaluate(() => performance.now());
  await page.getByRole("radio", { name: gene, exact: true }).click();

  // Wait for the render to complete after the selection.
  await page.waitForFunction(
    (prev) => {
      const m = (window as unknown as { __merfish?: MerfishHooks }).__merfish;
      return !!m && m.renderCount > prev && m.renderComplete != null;
    },
    before,
    { timeout: HOOK_TIMEOUT },
  );
  const end = await page.evaluate(
    () => (window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderComplete ?? performance.now(),
  );
  // Confirm the badge reflects the selection.
  await expect(page.locator(SELECTORS.selectedGeneBadge)).toHaveAttribute("data-gene", gene);
  return end - start;
}

async function renderCount(page: Page): Promise<number> {
  return page.evaluate(
    () => (window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderCount ?? 0,
  );
}

async function selectGeneSM(page: Page, gene: string): Promise<number> {
  const search = page.locator(SELECTORS.smGeneSearchInput);
  if (!(await search.isVisible())) {
    await page.locator(SELECTORS.smGenesButton).first().click();
  }
  await expect(search).toBeVisible();
  await search.fill(gene);

  const cb = page.getByRole("checkbox", { name: gene, exact: true });
  // Genes may be auto-selected on load; toggle off first so the check below
  // always triggers a fresh render we can time.
  if (await cb.isChecked()) {
    const off = await renderCount(page);
    await cb.uncheck();
    await page
      .waitForFunction(
        (p) => ((window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderCount ?? 0) > p,
        off,
        { timeout: HOOK_TIMEOUT },
      )
      .catch(() => {});
  }

  const before = await renderCount(page);
  const start = await page.evaluate(() => performance.now());
  await cb.check();

  await page.waitForFunction(
    (prev) => {
      const m = (window as unknown as { __merfish?: MerfishHooks }).__merfish;
      return !!m && m.renderCount > prev && m.renderComplete != null;
    },
    before,
    { timeout: HOOK_TIMEOUT },
  );
  const end = await page.evaluate(
    () => (window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderComplete ?? performance.now(),
  );
  return end - start;
}

/**
 * Verify the canvas actually rendered something (non-blank): sample the canvas
 * screenshot and assert it isn't a single uniform color.
 */
export async function assertCanvasRendered(page: Page, ds: DatasetEntry): Promise<void> {
  const canvas = page.locator(canvasSelector(ds)).first();
  const shot = await canvas.screenshot();
  const distinct = new Set<string>();
  // PNG bytes: cheap heuristic — count distinct 4-byte groups in a sample window.
  for (let i = 0; i < shot.length - 4; i += 997) {
    distinct.add(`${shot[i]},${shot[i + 1]},${shot[i + 2]},${shot[i + 3]}`);
    if (distinct.size > 5) break;
  }
  expect(distinct.size, "canvas appears blank (uniform pixels)").toBeGreaterThan(1);
}

/**
 * Drag-and-drop a single file onto its upload card (single-file formats only:
 * h5ad, single-molecule). Returns time-to-render in ms (page clock).
 */
export async function dragDropFile(page: Page, ds: DatasetEntry): Promise<number> {
  if (ds.format !== "h5ad" && ds.format !== "single-molecule") {
    throw new Error(`drag/drop only supported for single-file formats, not ${ds.format}`);
  }
  const home = ds.dataType === "single_molecule" ? "/?mode=sm" : "/";
  await page.goto(home, { waitUntil: "domcontentloaded" });

  const selector = fileInputSelector(ds.format);
  await expect(page.locator(selector)).toBeAttached();

  const bytes = fs.readFileSync(absPath(ds));
  const b64 = bytes.toString("base64");
  const fileName = absPath(ds).split("/").pop()!;

  const start = await page.evaluate(() => performance.now());
  await page.evaluate(
    async ({ sel, data, name }) => {
      const input = document.querySelector(sel) as HTMLInputElement;
      const card = input.parentElement as HTMLElement; // card div carries onDrop
      const bin = atob(data);
      const arr = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) arr[i] = bin.charCodeAt(i);
      const file = new File([arr], name);
      const dt = new DataTransfer();
      dt.items.add(file);
      card.dispatchEvent(new DragEvent("dragover", { bubbles: true, dataTransfer: dt }));
      card.dispatchEvent(new DragEvent("drop", { bubbles: true, dataTransfer: dt }));
    },
    { sel: selector, data: b64, name: fileName },
  );

  await page.waitForURL(viewerUrlGlob(ds), { timeout: HOOK_TIMEOUT });
  await waitForRender(page, ds);
  const end = await page.evaluate(
    () => (window as unknown as { __merfish?: MerfishHooks }).__merfish?.renderComplete ?? performance.now(),
  );
  return end - start;
}
