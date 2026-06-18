import {
  test,
  expect,
  datasetsByTier,
  annotate,
  loadLocalDataset,
  waitForDeStatsReady,
  SELECTORS,
} from "../../fixtures/merfish-fixtures";

/**
 * Differentially-expressed-genes (DEG) tests — exercise the panel end-to-end
 * across every single-cell tiny fixture. Tagged @quick.
 *
 *   h5ad / xenium / merscope: worker-eager priority-column compute
 *   chunked:                  adapter-backed lazy fetch from Python output
 *
 * Both paths emit `markDeStatsReady` (see lib/utils/test-hooks.ts), so the
 * test waits on the same hook regardless of pipeline.
 */

const datasets = datasetsByTier("tiny").filter(
  (d) => d.dataType === "single_cell",
);

for (const ds of datasets) {
  test(`@quick DEG: ${ds.id} (${ds.format})`, async ({ page }, testInfo) => {
    annotate(testInfo, ds.id, ds.format);

    await test.step("load dataset locally", async () => {
      await loadLocalDataset(page, ds);
    });

    await test.step("DEG stats become ready (hook fires)", async () => {
      const column = await waitForDeStatsReady(page);
      expect(column, `${ds.id}: deStatsColumn should be set`).toBeTruthy();
      testInfo.annotations.push({
        type: "deg",
        description: `column=${column}`,
      });
    });

    await test.step("open DEG panel", async () => {
      const button = page.locator(SELECTORS.degButton);

      await expect(
        button,
        `${ds.id}: DEG button should appear when stats are available`,
      ).toBeVisible({ timeout: 10_000 });
      await button.click();
      await expect(page.locator(SELECTORS.degPanel)).toBeVisible();
    });

    await test.step("pick a target celltype + verify rows", async () => {
      const panel = page.locator(SELECTORS.degPanel);

      // Auto mode is on by default, which disables the Autocomplete. Flip to
      // Manual so the test can pick a target deterministically.
      const targetAutoBtn = panel
        .locator('button:has-text("Auto")')
        .first();
      if (await targetAutoBtn.isVisible()) {
        await targetAutoBtn.click();
      }

      // HeroUI Autocomplete: clicking the input opens the listbox.
      const target = panel.getByLabel("Target celltype", { exact: false });
      await target.click();
      const firstOption = page.getByRole("option").first();
      await firstOption.waitFor({ state: "visible", timeout: 10_000 });
      const targetCelltype = (await firstOption.textContent())?.trim() ?? "";
      await firstOption.click();

      // Wait for ranked rows to populate.
      const rowSel = `${SELECTORS.degPanel} [data-testid="deg-row"]`;

      await expect(page.locator(rowSel).first()).toBeVisible({ timeout: 10_000 });
      const rowCount = await page.locator(rowSel).count();
      expect(rowCount, `${ds.id}: ranked DEG list should be non-empty`).toBeGreaterThan(0);

      testInfo.annotations.push({
        type: "deg",
        description: `target='${targetCelltype}' rows=${rowCount}`,
      });
    });

    await test.step("clicking a DEG row updates selectedGene", async () => {
      const firstRow = page
        .locator(`${SELECTORS.degPanel} [data-testid="deg-row"]`)
        .first();
      const gene = await firstRow.getAttribute("data-gene");

      expect(gene, `${ds.id}: top row should have a data-gene attribute`).toBeTruthy();
      await firstRow.click();
      await expect(page.locator(SELECTORS.selectedGeneBadge)).toHaveAttribute(
        "data-gene",
        gene!,
        { timeout: 10_000 },
      );
    });
  });
}
