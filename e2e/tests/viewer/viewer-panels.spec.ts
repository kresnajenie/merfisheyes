import {
  expect,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * The viewer's panels and columns behave in a window of any size: the image
 * and shapes panels close on a click outside, the control rail and the
 * legends scroll instead of running off the bottom, and Share shows the link
 * with its QR code. Tagged @quick: the tiny store, served by a local stand-in
 * for the bucket.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

const openStore = async (page: import("@playwright/test").Page) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, {
        cells: true,
        points: true,
        shapes: true,
        image: true,
      }),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
};

/** A click on the scene, away from every panel and button. */
const clickScene = async (page: import("@playwright/test").Page) => {
  const size = page.viewportSize()!;

  await page.mouse.click(size.width / 2, size.height - 40);
};

/** The shortest item in the rail: it scrolls, so none is squashed. */
const shortestRailItem = (page: import("@playwright/test").Page) =>
  page
    .getByTestId("control-rail")
    .evaluate((el) =>
      Math.min(
        ...Array.from(el.children, (c) => c.getBoundingClientRect().height),
      ),
    );

/** The legends end above the molecule status panel, never under it. */
const legendsAboveStatus = async (page: import("@playwright/test").Page) => {
  const legends = (await page.getByTestId("legends").first().boundingBox())!;
  const status = (await page
    .getByTestId("points-status")
    .first()
    .boundingBox())!;

  expect(legends.y + legends.height).toBeLessThanOrEqual(status.y);
};

test("@quick the image and shapes panels close on a click outside", async ({
  page,
}) => {
  await openStore(page);

  await test.step("image", async () => {
    await page.locator(SELECTORS.imageButton).click();
    await expect(page.getByTestId("image-panel")).toBeVisible();
    await clickScene(page);
    await expect(page.getByTestId("image-panel")).toBeHidden();
  });

  await test.step("shapes, beside its button", async () => {
    const button = page.getByTestId("shapes-button");

    await button.click();
    const panel = page.getByTestId("shapes-panel");

    await expect(panel).toBeVisible();
    const b = (await button.boundingBox())!;
    const p = (await panel.boundingBox())!;

    expect(p.x).toBeGreaterThanOrEqual(b.x + b.width);
    // A click inside the panel keeps it open
    await page.getByLabel("Show cell shapes").click();
    await expect(panel).toBeVisible();
    await clickScene(page);
    await expect(panel).toBeHidden();
  });
});

test("@quick the rail and the legends scroll in a short window", async ({
  page,
}) => {
  await page.setViewportSize({ width: 1280, height: 420 });
  await openStore(page);

  const fits = async (testId: string) => {
    const el = page.getByTestId(testId).first();

    await expect(el).toBeVisible({ timeout: 60_000 });
    const box = (await el.boundingBox())!;

    expect(
      box.y + box.height,
      `${testId} inside the window`,
    ).toBeLessThanOrEqual(420);

    return el;
  };

  await test.step("the control rail", async () => {
    const rail = await fits("control-rail");

    await expect(rail).toHaveAttribute("data-fade-scroll", "overflowing");
    expect(await shortestRailItem(page)).toBeGreaterThanOrEqual(24);
    // Scrolled to the end, its last button is on screen
    await rail.evaluate((el) => {
      el.scrollTop = el.scrollHeight;
    });
    const last = rail.locator("button").last();

    await expect(last).toBeInViewport();
  });

  await test.step("the legends", async () => {
    await fits("legends");
    await legendsAboveStatus(page);
  });
});

test("@quick a molecules store's rail and legends stay in a short window", async ({
  page,
}) => {
  await page.setViewportSize({ width: 1280, height: 480 });
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      storeUrl(server, { points: true }),
    )}`,
  );
  await expect(page.locator(SELECTORS.smCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  for (const testId of ["control-rail", "legends"]) {
    const el = page.getByTestId(testId).first();

    await expect(el).toBeVisible({ timeout: 60_000 });
    const box = (await el.boundingBox())!;

    expect(
      box.y + box.height,
      `${testId} inside the window`,
    ).toBeLessThanOrEqual(480);
  }
  const rail = page.getByTestId("control-rail");

  await rail.evaluate((el) => {
    el.scrollTop = el.scrollHeight;
  });
  expect(await shortestRailItem(page)).toBeGreaterThanOrEqual(24);
  await expect(rail.locator("button").last()).toBeInViewport();
  await legendsAboveStatus(page);
});

test("@quick Share shows the link to this view and its QR code", async ({
  page,
  context,
}) => {
  await context.grantPermissions(["clipboard-read", "clipboard-write"]);
  await openStore(page);
  // The molecules have their genes, so the view has something to write
  await expect(page.getByTestId("points-status").first()).toBeVisible({
    timeout: 60_000,
  });

  await page.getByRole("button", { name: "Share", exact: true }).click();
  const panel = page.getByTestId("share-panel");

  await expect(panel).toBeVisible();
  await expect(page.getByTestId("share-qr").locator("svg")).toBeVisible();
  const link = await page
    .getByRole("textbox", { name: "Link to this view" })
    .inputValue();

  // Read after the pending writes are flushed: the view as it is now
  expect(link).toBe(page.url());
  expect(await page.evaluate(() => navigator.clipboard.readText())).toBe(link);

  await clickScene(page);
  await expect(panel).toBeHidden();
  // The button toggles it too
  await page.getByRole("button", { name: "Share", exact: true }).click();
  await expect(panel).toBeVisible();
  await page.getByRole("button", { name: "Share", exact: true }).click();
  await expect(panel).toBeHidden();
});
