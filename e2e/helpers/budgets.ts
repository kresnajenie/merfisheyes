import type {
  MerfishFrameStats,
  MerfishLayerStats,
  MerfishMemoryStats,
} from "@/lib/utils/test-hooks";

import { type Page, expect } from "@playwright/test";

/**
 * Reading the budgets (docs/specs/BUDGETS.md) off a page: the layer hooks
 * (time to first content, settled view), the frame window, memory, and the
 * bytes the browser fetched per element kind.
 *
 * Bytes come from the Resource Timing API, which only reports sizes for
 * responses that allow it (`Timing-Allow-Origin`, which the fixture server
 * sends). For a real bucket the proxy's log is the source instead
 * (scripts/perf/real-store.ts).
 */

type Hooks = {
  layers: Record<string, MerfishLayerStats>;
  frames: MerfishFrameStats;
  memory: MerfishMemoryStats;
  renderComplete: number | null;
};

const hooks = (page: Page) =>
  page.evaluate(
    () => (window as unknown as { __merfish?: Hooks }).__merfish ?? null,
  );

export async function layer(
  page: Page,
  kind: string,
): Promise<MerfishLayerStats | null> {
  return (await hooks(page))?.layers[kind] ?? null;
}

/** Wait until a layer has its first content, and return that time. */
export async function firstContent(
  page: Page,
  kind: string,
  timeout = 120_000,
): Promise<number> {
  await expect
    .poll(async () => (await layer(page, kind))?.firstContentAt ?? null, {
      message: `${kind} should draw something`,
      timeout,
    })
    .not.toBeNull();

  return (await layer(page, kind))!.firstContentAt!;
}

/**
 * Wait until a layer settles after `since` (a view change it has already
 * been told about), and return the settled time.
 */
export async function settledAfter(
  page: Page,
  kind: string,
  since: number,
  timeout = 120_000,
): Promise<number> {
  const t0 = Date.now();

  for (;;) {
    const l = await layer(page, kind);

    if (l?.settledAt !== null && l!.settledAt! >= since) return l!.settledAt!;
    if (Date.now() - t0 > timeout) {
      throw new Error(
        `${kind} did not settle within ${timeout / 1000} s; last state: ${JSON.stringify(l)}`,
      );
    }
    await page.waitForTimeout(100);
  }
}

/** performance.now() on the page. */
export const pageNow = (page: Page) => page.evaluate(() => performance.now());

/** Bytes fetched so far for URLs whose path contains `needle`. */
export async function bytesFetched(
  page: Page,
  needle: string,
): Promise<number> {
  return page.evaluate((n) => {
    let sum = 0;

    for (const e of performance.getEntriesByType("resource")) {
      const r = e as PerformanceResourceTiming;

      if (r.name.includes(n)) sum += r.transferSize || r.encodedBodySize || 0;
    }

    return sum;
  }, needle);
}

/** Requests so far for URLs whose path contains `needle`. */
export async function requestsMade(
  page: Page,
  needle: string,
): Promise<number> {
  return page.evaluate(
    (n) =>
      performance.getEntriesByType("resource").filter((e) => e.name.includes(n))
        .length,
    needle,
  );
}

/** Start a fresh frame window. */
export async function resetFrames(page: Page): Promise<void> {
  await page.evaluate(() => {
    const h = (window as unknown as { __merfish?: Hooks }).__merfish;

    if (!h) return;
    h.frames = {
      count: 0,
      drawn: 0,
      over16ms: 0,
      over33ms: 0,
      worstMs: 0,
      worstRenderMs: 0,
      renderMs: 0,
      since: performance.now(),
    };
  });
}

export async function frames(page: Page): Promise<MerfishFrameStats> {
  return (await hooks(page))!.frames;
}

export async function memory(page: Page): Promise<MerfishMemoryStats> {
  return (await hooks(page))!.memory;
}

/** The scene canvas of a page (cells first, then molecules). */
export async function sceneCanvas(page: Page) {
  const sc = page.locator('[data-testid="sc-scene-canvas"] canvas').first();

  if (await sc.count()) return sc;

  return page.locator('[data-testid="sm-scene-canvas"] canvas').first();
}

/**
 * Zoom the view by `ticks` wheel notches at the canvas centre (negative
 * zooms in), as ONE wheel event: OrbitControls scales by 0.95 per notch of
 * the event's delta, so 32 notches are a fifth of the width and 64 about a
 * twenty-seventh. One event means one replan, so "view settled" measures
 * one view; a hand zooming in steps also queues every intermediate view's
 * reads, which is measured separately (`zoomSteps`).
 */
export async function zoom(page: Page, ticks: number): Promise<void> {
  const canvas = await sceneCanvas(page);
  const box = (await canvas.boundingBox())!;

  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.wheel(0, ticks * 100);
}

/** The same zoom as a hand does it: `events` wheel events 50 ms apart. */
export async function zoomSteps(
  page: Page,
  ticks: number,
  events = 4,
): Promise<void> {
  const canvas = await sceneCanvas(page);
  const box = (await canvas.boundingBox())!;
  const per = Math.round(ticks / events);

  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  for (let i = 0; i < events; i++) {
    await page.mouse.wheel(0, per * 100);
    await page.waitForTimeout(50);
  }
}

/**
 * Drag the view by a fraction of the canvas width over `steps` mouse moves,
 * about one second in all — the scripted pan the frame budget is read on.
 */
export async function pan(
  page: Page,
  fraction = 0.3,
  steps = 30,
): Promise<void> {
  const canvas = await sceneCanvas(page);
  const box = (await canvas.boundingBox())!;
  const x0 = box.x + box.width * 0.35;
  const y = box.y + box.height / 2;
  const dx = (box.width * fraction) / steps;

  await page.mouse.move(x0, y);
  await page.mouse.down();
  for (let i = 1; i <= steps; i++) {
    await page.mouse.move(x0 + dx * i, y);
    await page.waitForTimeout(1000 / steps);
  }
  await page.mouse.up();
}

/** Move the pointer across the canvas for about a second (hover cost). */
export async function sweep(page: Page, steps = 40): Promise<void> {
  const canvas = await sceneCanvas(page);
  const box = (await canvas.boundingBox())!;
  const y = box.y + box.height / 2;

  for (let i = 0; i <= steps; i++) {
    await page.mouse.move(box.x + (box.width * i) / steps, y);
    await page.waitForTimeout(1000 / steps);
  }
}

/**
 * Move a HeroUI slider (`aria-label`) to the given fraction of its range
 * with the keyboard: focus the thumb and press the arrow key. One press
 * is one step, which is what a slider tick costs.
 */
export async function nudgeSlider(page: Page, label: string): Promise<void> {
  const thumb = page
    .locator(`[aria-label="${label}"]`)
    .locator('[role="slider"], input[type="range"]')
    .first();

  await thumb.focus();
  await page.keyboard.press("ArrowUp");
}
