import path from "node:path";

import { expect, SELECTORS, test } from "../../fixtures/merfish-fixtures";

/**
 * Files the browser preview used to misread (docs/plans/PREVIEW-SAVE.md,
 * step 0; scripts/testdata/make-tiny-h5ad-edge.py), dropped on the homepage
 * as a person would. The parse itself is unit-tested
 * (tests/unit/preview-parse.test.ts); this checks the scene takes the
 * result: a cell with no position (NaN) is simply not drawn, and nothing
 * downstream — bounds, camera, hover — complains. Tagged @quick: local
 * files, no database or S3.
 */

const H5AD = path.join(__dirname, "..", "..", "..", "tests", "data", "tiny", "h5ad");

async function drop(page: import("@playwright/test").Page, name: string) {
  const errors: string[] = [];

  // Uncaught exceptions, and what three.js says about NaN positions
  // (computeBoundingSphere) — not the homepage's unrelated hydration warning
  page.on("console", (m) => {
    if (m.type() === "error" && /NaN|BoundingSphere/.test(m.text())) {
      errors.push(m.text());
    }
  });
  page.on("pageerror", (e) => errors.push(e.message));
  await page.goto("/", { waitUntil: "domcontentloaded" });
  await page.locator("#sc-file-input-h5ad").setInputFiles(path.join(H5AD, name));
  await page.waitForURL("**/viewer/from-local", { timeout: 120_000 });
  await expect(page.locator(SELECTORS.scCanvas).first()).toBeVisible({
    timeout: 120_000,
  });
  await page.waitForFunction(
    () => (window as any).__merfish?.renderComplete != null,
    undefined,
    { timeout: 120_000 },
  );

  return errors;
}

test("@quick cells with no position are kept, not drawn, and break nothing", async ({
  page,
}) => {
  const errors = await drop(page, "missing_xy.h5ad");

  // Hover across the scene: the picking grid must cope with the NaN cells
  const box = (await page.locator(SELECTORS.scCanvas).first().boundingBox())!;

  for (let i = 1; i < 10; i++) {
    await page.mouse.move(box.x + (box.width * i) / 10, box.y + box.height / 2);
  }
  await page.waitForTimeout(300);

  expect(errors).toEqual([]);
});

test("@quick an .h5ad written by anndata 0.13 lists its genes", async ({ page }) => {
  const errors = await drop(page, "nullable.h5ad");

  await page.getByRole("button", { name: "Gene", exact: true }).click();
  for (const gene of ["Gad1", "Sst", "Pvalb", "Vip"]) {
    await expect(page.getByText(gene, { exact: true }).first()).toBeVisible();
  }
  expect(errors).toEqual([]);
});
