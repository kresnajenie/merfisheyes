import {
  expect,
  imageUrl,
  libraryStoreUrl,
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

test("@quick the image panel scrolls through many channels in a short window", async ({
  page,
}) => {
  // The CODEX-shaped store: six channels, taller than this window
  await page.setViewportSize({ width: 1280, height: 600 });
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      libraryStoreUrl(server, "codex"),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.locator(SELECTORS.imageButton).click();
  const panel = page.getByTestId("image-panel");

  await expect(panel).toBeVisible();
  const box = (await panel.boundingBox())!;

  expect(
    box.y + box.height,
    "image panel inside the window",
  ).toBeLessThanOrEqual(600);
  const list = panel.locator("[data-fade-scroll]");

  await expect(list).toHaveAttribute("data-fade-scroll", "overflowing");
  await list.evaluate((el) => {
    el.scrollTop = el.scrollHeight;
  });
  // The last channel's last control is reachable
  await expect(panel.getByRole("slider").last()).toBeInViewport();
});

test("@quick a channel's slider range can be narrowed, typed past and kept", async ({
  page,
}) => {
  await page.goto(
    `/spatialdata-viewer/from-s3?url=${encodeURIComponent(
      libraryStoreUrl(server, "codex"),
    )}`,
  );
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.locator(SELECTORS.imageButton).click();
  const panel = page.getByTestId("image-panel");
  const box = (name: string) =>
    panel.getByRole("textbox", { name: `DAPI ${name}`, exact: true });
  const slider = panel.getByRole("slider", {
    name: "DAPI maximum",
    exact: true,
  });
  const type = async (name: string, value: string) => {
    await box(name).fill(value);
    await box(name).press("Enter");
  };

  // Narrowed: the sliders span it, and the window is pulled inside
  await type("range end", "150");
  await expect(slider).toHaveAttribute("max", "150");
  expect(Number(await box("maximum value").inputValue())).toBeLessThanOrEqual(
    150,
  );

  // A value typed past the range sets it and widens the range to hold it
  await type("maximum value", "300");
  await expect(box("maximum value")).toHaveValue("300");
  await expect(box("range end")).toHaveValue("300");
  await expect(slider).toHaveValue("300");

  // In the link: a reload keeps the range and the window
  await expect
    .poll(() => new URL(page.url()).searchParams.get("img"))
    .not.toBeNull();
  await page.reload();
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.locator(SELECTORS.imageButton).click();
  await expect(box("range end")).toHaveValue("300");
  await expect(box("maximum value")).toHaveValue("300");
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

// An image on its own has the buttons the other viewers have, in both places
// it opens
for (const shown of [
  {
    name: "/image-viewer",
    path: "image-viewer",
    url: () => imageUrl(server, 1000),
  },
  {
    name: "an image-only store",
    path: "spatialdata-viewer",
    url: () => storeUrl(server, { image: true }),
  },
]) {
  test(`@quick ${shown.name} has Hide UI, Share and Split screen`, async ({
    page,
    context,
  }) => {
    await context.grantPermissions(["clipboard-read", "clipboard-write"]);
    await page.goto(
      `/${shown.path}/from-s3?url=${encodeURIComponent(shown.url())}`,
    );
    await expect(page.locator(SELECTORS.imageCanvas).first()).toBeVisible({
      timeout: 120_000,
    });
    const share = page.getByRole("button", { name: "Share", exact: true });

    // Share: the link to this view and its QR code
    await share.click();
    await expect(page.getByTestId("share-qr").locator("svg")).toBeVisible();
    expect(
      await page
        .getByRole("textbox", { name: "Link to this view" })
        .inputValue(),
    ).toBe(page.url());
    await page.getByRole("button", { name: "Close share panel" }).click();

    // Hide UI hides the buttons; H brings them back
    await page.getByRole("button", { name: "Hide UI" }).click();
    await expect(share).toBeHidden();
    await page.keyboard.press("h");
    await expect(share).toBeVisible();

    // Split screen: the same image again on the right
    await page.getByRole("button", { name: "Split screen" }).click();
    await page.getByText("Same Dataset").click();
    await expect(page.locator(SELECTORS.imageCanvas)).toHaveCount(2, {
      timeout: 120_000,
    });
  });
}
