import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Moving the 2D view (lib/webgl/scene-manager.ts): the wheel zooms about the
 * cursor, and one finger pans. Tagged @quick: the tiny store, served by a
 * local stand-in for the bucket.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

type Page = import("@playwright/test").Page;

const openStore = async (page: Page) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, { cells: true }),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );
};

/** Camera position x, y, z, then fov (degrees) and aspect. */
const camera = (page: Page) =>
  page.evaluate(() => (window as any).__merfish.camera() as number[]);

/** Poll until the camera has stopped moving (damping), then return it. */
const settledCamera = async (page: Page) => {
  let last = await camera(page);

  await expect
    .poll(async () => {
      const now = await camera(page);
      const still = now.every((v, i) => Math.abs(v - last[i]) < 1e-6);

      last = now;

      return still;
    })
    .toBe(true);

  return last;
};

/** The world point (on the data plane) under a point of the canvas, as fractions of it. */
const worldAt = (cam: number[], fx: number, fy: number) => {
  const [x, y, z, fov, aspect] = cam;
  const half = z * Math.tan((fov * Math.PI) / 360);

  return [x + (fx * 2 - 1) * half * aspect, y - (fy * 2 - 1) * half];
};

test.describe("@quick moving the 2D view", () => {
  test("the wheel zooms about the cursor", async ({ page }) => {
    await openStore(page);

    const box = (await page.locator(SELECTORS.scCanvas).first().boundingBox())!;
    const fx = 0.7;
    const fy = 0.3;
    const before = await settledCamera(page);

    await page.mouse.move(box.x + box.width * fx, box.y + box.height * fy);
    for (let i = 0; i < 4; i++) await page.mouse.wheel(0, -120);

    const after = await settledCamera(page);

    // Closer, and the point under the cursor has not moved
    expect(after[2]).toBeLessThan(before[2] * 0.9);

    const [bx, by] = worldAt(before, fx, fy);
    const [ax, ay] = worldAt(after, fx, fy);
    const tolerance = before[2] * 0.01;

    expect(Math.abs(ax - bx)).toBeLessThan(tolerance);
    expect(Math.abs(ay - by)).toBeLessThan(tolerance);
    // ... so the middle of the view has
    expect(Math.abs(after[0] - before[0])).toBeGreaterThan(tolerance);
  });

  test("one finger pans", async ({ page }) => {
    await openStore(page);

    const box = (await page.locator(SELECTORS.scCanvas).first().boundingBox())!;
    const before = await settledCamera(page);
    const cdp = await page.context().newCDPSession(page);
    const x = box.x + box.width / 2;
    const y = box.y + box.height / 2;
    const touch = (type: string, dx: number) =>
      cdp.send("Input.dispatchTouchEvent", {
        type: type as "touchStart",
        touchPoints: type === "touchEnd" ? [] : [{ x: x + dx, y }],
      });

    await touch("touchStart", 0);
    for (let dx = 20; dx <= 100; dx += 20) await touch("touchMove", dx);
    await touch("touchEnd", 0);

    const after = await settledCamera(page);

    // The view followed the finger: dragging right moves the camera left
    expect(after[0]).toBeLessThan(before[0]);
    expect(after[2]).toBeCloseTo(before[2], 3);
  });
});
