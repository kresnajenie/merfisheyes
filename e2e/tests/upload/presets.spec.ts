import { PrismaClient } from "@prisma/client";
import { type Browser, type Page, expect, test } from "@playwright/test";

import { signInDev } from "../../helpers/auth";
import { type FixtureServer } from "../../helpers/fixture-server";
import { startSpatialDataServer, visiumHdUrl } from "../../helpers/spatialdata";

/**
 * Technology presets, by someone who is not an admin
 * (docs/specs/ADR-technology-presets.md): a super admin grants the
 * permission; its holder is still refused on every admin route; they save a
 * version at /presets, preview it by its hash on a store nobody has
 * published anything for, publish it, publish another, roll back, and turn
 * the preset off.
 *
 * Tagged @upload: it needs the suite's database (no S3). E2E_AGENT=1 runs it
 * against any migrated database the dev server's DATABASE_URL points at.
 */

const ENABLED = process.env.E2E_UPLOAD === "1" || process.env.E2E_AGENT === "1";
const TECHNOLOGY = "Visium HD";

/** Every route under app/api/admin, with the method it answers to. */
const ADMIN_ROUTES: [string, string][] = [
  ["GET", "/api/admin/catalog"],
  ["POST", "/api/admin/catalog"],
  ["GET", "/api/admin/catalog/export"],
  ["POST", "/api/admin/catalog/import"],
  ["GET", "/api/admin/catalog/x"],
  ["PATCH", "/api/admin/catalog/x"],
  ["DELETE", "/api/admin/catalog/x"],
  ["GET", "/api/admin/community"],
  ["POST", "/api/admin/community/x"],
  ["GET", "/api/admin/failures"],
  ["POST", "/api/admin/failures/x/retry"],
  ["GET", "/api/admin/processing"],
  ["PATCH", "/api/admin/processing/x"],
  ["GET", "/api/admin/users"],
  ["PATCH", "/api/admin/users/x"],
];

let server: FixtureServer;
const prisma = new PrismaClient();

/** Presets are append-only through the app; the test database is put back. */
async function clearPresets() {
  await prisma.technologyPresetPublish.deleteMany({
    where: { technology: TECHNOLOGY },
  });
  await prisma.technologyPresetVersion.deleteMany({
    where: { technology: TECHNOLOGY },
  });
}

test.describe("@upload technology presets", () => {
  test.skip(
    !ENABLED,
    "Set E2E_UPLOAD=1 (docker infra) or E2E_AGENT=1 (any migrated database).",
  );

  test.beforeAll(async () => {
    server = await startSpatialDataServer();
    await clearPresets();
  });
  test.afterAll(async () => {
    await clearPresets();
    await prisma.$disconnect();
    await server.close();
  });

  test("granted, refused on admin, save, preview, publish, roll back", async ({
    page,
    browser,
  }) => {
    test.setTimeout(600_000);
    const run = Date.now();
    const editor = `e2e-presets-${run}@example.com`;
    const viewerLink = `/spatialdata-viewer/from-s3?url=${encodeURIComponent(visiumHdUrl(server))}`;

    /** The opacity a fresh visitor's open of a link shows. */
    const opensWith = async (link: string, b: Browser = browser) => {
      const context = await b.newContext();
      const p = await context.newPage();

      await p.goto(link);
      const slider = p
        .getByTestId("cell-opacity")
        .getByRole("slider", { name: "Cell opacity" });

      await expect(slider).toBeVisible({ timeout: 120_000 });
      const value = await slider.inputValue();

      await context.close();

      return { opacity: value };
    };

    await signInDev(page, editor);
    const editorId = (
      await prisma.user.findUniqueOrThrow({ where: { email: editor } })
    ).id;

    await test.step("without the permission: no page, no writes", async () => {
      await page.goto("/presets");
      await expect(page).toHaveURL(/\/$/);
      const res = await page.request.post("/api/presets/versions", {
        data: { technology: TECHNOLOGY, config: {} },
      });

      expect(res.status()).toBe(403);
    });

    await test.step("a super admin grants it on the users route", async () => {
      const context = await browser.newContext();
      const admin: Page = await context.newPage();
      const email = `e2e-presets-admin-${run}@example.com`;

      await signInDev(admin, email);
      await prisma.user.update({
        where: { email },
        data: { role: "SUPER_ADMIN" },
      });

      const bad = await admin.request.patch(`/api/admin/users/${editorId}`, {
        data: { permissions: ["everything"] },
      });

      expect(bad.status()).toBe(400);
      const res = await admin.request.patch(`/api/admin/users/${editorId}`, {
        data: { permissions: ["presets:write"] },
      });

      expect(res.ok()).toBeTruthy();
      await context.close();
    });

    await test.step("its holder is refused on every admin route", async () => {
      for (const [method, path] of ADMIN_ROUTES) {
        const res = await page.request.fetch(path, { method, data: {} });

        expect(res.status(), `${method} ${path}`).toBe(403);
      }
      await page.goto("/admin");
      await expect(page).toHaveURL(/\/$/);
    });

    await test.step("a malformed config is refused by the schema", async () => {
      const res = await page.request.post("/api/presets/versions", {
        data: { technology: TECHNOLOGY, config: { colormap: "rainbow" } },
      });

      expect(res.status()).toBe(400);
      expect((await res.json()).error).toContain("/colormap");
    });

    const saveDraft = async (opacity: string) => {
      await page.getByTestId("preset-colormap").selectOption("plasma");
      await page.getByTestId("preset-hide-zero").selectOption("false");
      await page.getByTestId("preset-opacity").fill(opacity);
      await page.getByTestId("preset-tested-on").fill(viewerLink);
      await page.getByTestId("preset-save").click();
      const message = page.getByTestId("preset-message");

      await expect(message).toContainText("Saved as draft");

      return (await message.innerText()).match(/draft ([0-9a-f]{16})/)![1];
    };
    let first = "";

    await test.step("saves a version; nothing changes for visitors yet", async () => {
      await page.goto("/presets");
      await expect(page.getByTestId("preset-in-force")).toContainText(
        "built-in",
      );
      first = await saveDraft("0.5");
      await expect(page.getByTestId(`preset-version-${first}`)).toContainText(
        "colormap plasma, zero counts shown, opacity 50%",
      );
      expect((await opensWith(viewerLink)).opacity).toBe("0.8");
    });

    await test.step("the preview link shows it, signed out", async () => {
      const link = await page
        .getByTestId(`preset-preview-${first}`)
        .getAttribute("href");

      expect(link).toContain(`preset=${first}`);
      expect((await opensWith(link!)).opacity).toBe("0.5");
    });

    await test.step("publishing it changes what a plain link opens with", async () => {
      await page.getByTestId(`preset-publish-${first}`).click();
      await expect(page.getByTestId("preset-in-force")).toContainText(first);
      expect((await opensWith(viewerLink)).opacity).toBe("0.5");
    });

    await test.step("a second version, then back to the first", async () => {
      const second = await saveDraft("0.3");

      expect(second).not.toBe(first);
      await page.getByTestId(`preset-publish-${second}`).click();
      await expect(page.getByTestId("preset-in-force")).toContainText(second);
      expect((await opensWith(viewerLink)).opacity).toBe("0.3");

      await page.getByTestId(`preset-publish-${first}`).click();
      await expect(page.getByTestId("preset-in-force")).toContainText(first);
      expect((await opensWith(viewerLink)).opacity).toBe("0.5");
    });

    await test.step("turned off: the viewer's ordinary defaults; the log has it all", async () => {
      await page.getByTestId("preset-turn-off").click();
      await expect(page.getByTestId("preset-in-force")).toContainText(
        "turned off",
      );
      expect((await opensWith(viewerLink)).opacity).toBe("1");
      await expect(
        page.getByTestId("preset-history").locator("li"),
      ).toHaveCount(4);

      const saved = await prisma.technologyPresetVersion.findUniqueOrThrow({
        where: { hash: first },
      });

      expect(saved.createdById).toBe(editorId);
      expect(saved.testedOn).toBe(viewerLink);
      expect(saved.appVersion).toMatch(/^\d+\.\d+\.\d+/);
    });
  });
});
