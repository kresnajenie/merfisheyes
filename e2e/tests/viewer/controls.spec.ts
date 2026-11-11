import {
  expect,
  getDataset,
  litFraction,
  loadLocalDataset,
  SELECTORS,
  startSpatialDataServer,
  storeUrl,
  test,
} from "../../fixtures/merfish-fixtures";
import { type FixtureServer } from "../../helpers/fixture-server";

/**
 * Driving what the scene draws from the keyboard (components/shortcut-manager.tsx)
 * and recolouring several legend chips at once
 * (lib/hooks/useLegendSelection.ts, lib/utils/color-assign.ts). Tagged
 * @quick: the tiny store, served by a local stand-in for the bucket, so no
 * database or S3 is needed.
 */

let server: FixtureServer;

test.beforeAll(async () => {
  server = await startSpatialDataServer();
});
test.afterAll(() => server.close());

/** The tiny store with everything in it, so every shortcut has a layer. */
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

test.describe("@quick keyboard shortcuts", () => {
  test("a key hides a layer and shows it again", async ({ page }) => {
    await openStore(page);

    const scene = SELECTORS.scCanvas;
    // Inside the image, which the tiny store draws as a block in the middle
    // of the scene: lit where it is drawn, nearly black where it is not.
    const insideImage = { x0: 0.42, y0: 0.35, x1: 0.6, y1: 0.7 };
    const lit = () => litFraction(page, scene, insideImage);

    await expect
      .poll(lit, { timeout: 60_000 })
      .toBeGreaterThan(0.5);
    const withImage = await lit();

    await page.keyboard.press("i");
    await expect.poll(lit, { timeout: 15_000 }).toBeLessThan(0.2);

    await page.keyboard.press("i");
    await expect
      .poll(lit, { timeout: 15_000 })
      .toBeGreaterThan(withImage * 0.9);
  });

  test("the list of shortcuts opens with ? and says what this view has", async ({
    page,
  }) => {
    await openStore(page);

    const list = page.getByRole("dialog", { name: "Keyboard shortcuts" });

    await expect(list).toBeHidden();
    await page.keyboard.press("?");
    await expect(list).toBeVisible();
    // Every layer is in this store, so none is marked missing
    await expect(list.getByText("Cells")).toBeVisible();
    await expect(list.getByText("Image")).toBeVisible();
    await expect(list.getByText("not in this view")).toHaveCount(0);

    await page.keyboard.press("Escape");
    await expect(list).toBeHidden();
  });

  test("a key typed into a search box is text, not a shortcut", async ({
    page,
  }) => {
    await openStore(page);

    await page.click(SELECTORS.geneModeButton);
    const search = page.locator(SELECTORS.geneSearchInput).first();

    await expect(search).toBeVisible();
    await search.fill("");
    await search.press("i");
    await search.press("c");

    await expect(search).toHaveValue("ic");
    await expect(
      page.getByRole("dialog", { name: "Keyboard shortcuts" }),
    ).toBeHidden();
  });
});

test.describe("@quick overlaying a preview", () => {
  test("a dropped .h5ad can have molecules checked over it", async ({
    page,
  }) => {
    // A folder of processed parts already offered this; a plain file is a
    // preview too, and its coordinates are raw, so an overlay lines up.
    await loadLocalDataset(page, getDataset("h5ad-tiny"));

    const overlay = page.locator(SELECTORS.overlayButton);

    await expect(overlay).toBeVisible({ timeout: 30_000 });
    await overlay.click();
    await expect(
      page.getByRole("dialog").getByText("Manage overlay"),
    ).toBeVisible();
  });
});

test.describe("@quick recolouring several legend chips", () => {
  test("shift-click picks a range, and Spread gives them different colours", async ({
    page,
  }) => {
    await openStore(page);

    await page.click(SELECTORS.celltypeModeButton);

    // Pick the first three celltypes in the panel, so the legend has chips
    const options = page.getByRole("checkbox");

    await expect(options.first()).toBeVisible({ timeout: 30_000 });
    const pickCount = Math.min(3, await options.count());

    for (let i = 0; i < pickCount; i++) await options.nth(i).click();

    const legend = page.getByTestId("legends").first();
    const chips = legend.getByTestId("celltype-chip");

    await expect(chips).toHaveCount(pickCount);

    // One colour for all of them, then back to different ones
    const colorOf = (i: number) =>
      chips.nth(i).evaluate((el) => getComputedStyle(el).backgroundColor);

    await chips.first().click();
    await chips.nth(pickCount - 1).click({ modifiers: ["Shift"] });

    const actions = legend.getByRole("toolbar");

    await expect(actions).toContainText(`${pickCount} picked`);

    await actions.getByRole("button", { name: "Spread" }).click();
    const spread = await Promise.all(
      Array.from({ length: pickCount }, (_, i) => colorOf(i)),
    );

    expect(new Set(spread).size).toBe(pickCount);
    // The selection is dropped once spread, so the chips read normally again
    await expect(actions).toBeHidden();
  });
});
