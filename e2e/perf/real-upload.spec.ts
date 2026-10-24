import fs from "node:fs";
import path from "node:path";

import {
  test,
  expect,
  readHooks,
  recordMetric,
} from "../fixtures/merfish-fixtures";
import { signInDev } from "../helpers/auth";

/**
 * Real datasets uploaded from the homepage's "Upload & process on server"
 * mode, timed, then opened from what was saved. Two kinds:
 *
 *   MERFISHEYES_REAL_STORES   folders of already-processed data (a prepared
 *                             SpatialData store): saved as they are, no worker
 *   MERFISHEYES_REAL_RAW      raw inputs for the worker: `sc:<folder>` (cells,
 *                             with molecules beside them) or `sm:<file>`
 *                             (molecules only; the worker reads the columns from the header)
 *
 * Comma-separated. Raw uploads need a real processing queue, so they belong
 * on staging (`E2E_BASE_URL` = a server with the staging .env), never on the
 * local test stack. Each test attaches a "summary" JSON. Tagged @perf.
 */

const list = (name: string) =>
  (process.env[name] ?? "").split(",").filter(Boolean);
const STORES = list("MERFISHEYES_REAL_STORES");
const RAW = list("MERFISHEYES_REAL_RAW");
/** Already-saved datasets to open again (ids), e.g. after a viewer fix. */
const REOPEN = list("MERFISHEYES_REOPEN");
const EMAIL =
  process.env.MERFISHEYES_UPLOAD_EMAIL ?? "e2e-real-upload@example.com";

type Page = import("@playwright/test").Page;

/** Every file under a folder, and their total size. */
function sizeOf(target: string): { files: number; bytes: number } {
  if (fs.statSync(target).isFile()) {
    return { files: 1, bytes: fs.statSync(target).size };
  }
  let files = 0;
  let bytes = 0;

  for (const rel of fs.readdirSync(target, { recursive: true }) as string[]) {
    const stat = fs.statSync(path.join(target, rel));

    if (stat.isFile()) {
      files++;
      bytes += stat.size;
    }
  }

  return { files, bytes };
}

type Mine = {
  id: string;
  title: string;
  status: string;
  numCells?: number | null;
  numGenes?: number | null;
  viewerUrl?: string | null;
};

async function mine(page: Page, title: string): Promise<Mine[]> {
  const out: Mine[] = [];

  for (const type of ["single_cell", "single_molecule", "spatialdata"]) {
    const res = await page.request.get(`/api/ingest/mine?type=${type}`);

    if (!res.ok()) continue;
    for (const d of (await res.json()).datasets as Mine[]) {
      if (d.title.startsWith(title)) out.push(d);
    }
  }

  return out;
}

async function serverMode(page: Page, molecules: boolean) {
  await page.goto("/", { waitUntil: "domcontentloaded" });
  // The homepage remembers the last mode, so it is always chosen
  await page
    .getByRole("button", {
      name: molecules ? /single molecule/i : /single cell/i,
    })
    .first()
    .click();
  await page.getByRole("radio", { name: "Upload & process on server" }).click();
}

test.describe("@perf real upload: prepared stores", () => {
  test.skip(STORES.length === 0, "Set MERFISHEYES_REAL_STORES to run.");

  for (const folder of STORES) {
    test(folder, async ({ page }, info) => {
      test.setTimeout(4 * 3_600_000);
      const title = `real-${path.basename(folder)}-${Date.now()}`;
      const size = sizeOf(folder);

      await signInDev(page, EMAIL);
      await serverMode(page, false);
      const t0 = Date.now();

      await page.locator("#sc-file-input-folder").setInputFiles(folder);
      const dialog = page.getByLabel(/Project name|Dataset name/);

      await expect(dialog).toBeVisible({ timeout: 600_000 });
      const toDialog = Date.now() - t0;

      await dialog.fill(title);
      const t1 = Date.now();

      await page.getByRole("button", { name: "Upload", exact: true }).click();
      await expect(page.getByTestId("bundle-saved")).toBeVisible({
        timeout: 4 * 3_600_000,
      });
      const toSaved = Date.now() - t1;

      await page.getByText("Open the saved dataset").click();
      await page.waitForURL(/\/ds_/, { timeout: 120_000 });
      const t2 = Date.now();

      // Drawn: cells on screen, or molecules from either molecule layer
      const drawn = () =>
        page.evaluate(() => {
          const m = (window as any).__merfish;

          return (
            (m?.stats?.pointCount ?? 0) +
            (m?.layers?.points?.drawn ?? 0) +
            (m?.layers?.molecules?.drawn ?? 0)
          );
        });

      await expect
        .poll(drawn, { timeout: 600_000, intervals: [2_000] })
        .toBeGreaterThan(0)
        .catch(() => undefined);
      const hooks = await readHooks(page);

      await recordMetric(info, {
        datasetId: title,
        format: "prepared-store",
        metric: "uploadTime",
        ms: toSaved,
      });
      await info.attach("summary", {
        contentType: "application/json",
        body: JSON.stringify({
          folder,
          url: page.url(),
          files: size.files,
          gigabytes: size.bytes / 1e9,
          toDialogSeconds: toDialog / 1000,
          uploadSeconds: toSaved / 1000,
          megabytesPerSecond: size.bytes / 1e6 / (toSaved / 1000),
          openToFirstDrawSeconds: (Date.now() - t2) / 1000,
          stats: hooks?.stats,
          layers: hooks?.layers,
          registered: await mine(page, title),
        }),
      });
      expect(await drawn(), "the saved dataset draws").toBeGreaterThan(0);
    });
  }
});

test.describe("@perf real upload: open again", () => {
  test.skip(REOPEN.length === 0, "Set MERFISHEYES_REOPEN to run.");

  for (const id of REOPEN) {
    test(id, async ({ page }, info) => {
      test.setTimeout(900_000);
      await signInDev(page, EMAIL);
      // The dataset's own viewer link, as the account page has it
      const listed = (await mine(page, "")).find((d) => d.id === id);

      expect(
        listed?.viewerUrl,
        `${id} is listed with a viewer link`,
      ).toBeTruthy();
      const t0 = Date.now();

      await page.goto(listed!.viewerUrl!);
      const drawn = () =>
        page.evaluate(() => {
          const m = (window as any).__merfish;

          return (
            (m?.stats?.pointCount ?? 0) +
            (m?.layers?.points?.drawn ?? 0) +
            (m?.layers?.molecules?.drawn ?? 0)
          );
        });

      await expect
        .poll(drawn, { timeout: 600_000, intervals: [2_000] })
        .toBeGreaterThan(0);
      const toFirstDraw = (Date.now() - t0) / 1000;

      // The molecules over the cells, when there are any, also draw
      await expect
        .poll(async () => (await readHooks(page))?.layers.points?.drawn ?? 0, {
          timeout: 120_000,
          intervals: [2_000],
        })
        .toBeGreaterThan(0)
        .catch(() => undefined);
      const hooks = await readHooks(page);

      await info.attach("summary", {
        contentType: "application/json",
        body: JSON.stringify({
          id,
          openToFirstDrawSeconds: toFirstDraw,
          stats: hooks?.stats,
          layers: hooks?.layers,
        }),
      });
    });
  }
});

test.describe("@perf real upload: raw for the worker", () => {
  test.skip(RAW.length === 0, "Set MERFISHEYES_REAL_RAW to run.");

  for (const entry of RAW) {
    const [kind, target] = [entry.slice(0, 2), entry.slice(3)];

    test(entry, async ({ page }, info) => {
      test.setTimeout(6 * 3_600_000);
      const title = `real-${path.basename(target)}-${Date.now()}`;
      const size = sizeOf(target);

      await signInDev(page, EMAIL);
      await serverMode(page, kind === "sm");
      const t0 = Date.now();

      // Server mode asks no columns: the worker detects them from the header
      await page
        .locator(
          kind === "sm" ? "#sm-file-input-file" : "#sc-file-input-folder",
        )
        .setInputFiles(target);
      await expect(page.getByText("Name your dataset")).toBeVisible({
        timeout: 600_000,
      });
      await page.getByLabel(/Project name|Dataset name/).fill(title);
      // Whether molecules beside the cells were recognised and offered
      const offeredMolecules = await page
        .getByText("Also upload the single-molecule data")
        .isVisible();
      const t1 = Date.now();

      // "Start upload", or "Upload both" / "Upload all three" with extras
      await page
        .getByRole("button", { name: /^(Start upload|Upload)/ })
        .last()
        .click();
      await page.waitForURL("**/explore", { timeout: 6 * 3_600_000 });
      // The raw bytes are up once the bar says so
      await expect(
        page.getByText(/Uploaded|Processing|Queued/).first(),
      ).toBeVisible({
        timeout: 6 * 3_600_000,
      });
      const uploadSeconds = (Date.now() - t1) / 1000;

      // Then the worker. The list leaves failed datasets out, so their ids
      // are taken while queued and each one's own status is polled
      const ids = new Set<string>();

      await expect
        .poll(
          async () => {
            for (const d of await mine(page, title)) ids.add(d.id);

            return ids.size;
          },
          { timeout: 600_000, intervals: [10_000] },
        )
        .toBeGreaterThan(0);
      let datasets: Mine[] = [];

      await expect
        .poll(
          async () => {
            datasets = await Promise.all(
              Array.from(ids, async (id) => ({
                id,
                title,
                ...(await (
                  await page.request.get(`/api/ingest/${id}/status`)
                ).json()),
              })),
            );

            return datasets.every((d) =>
              ["COMPLETE", "FAILED"].includes(d.status),
            );
          },
          { timeout: 6 * 3_600_000, intervals: [30_000] },
        )
        .toBe(true);

      await info.attach("summary", {
        contentType: "application/json",
        body: JSON.stringify({
          entry,
          files: size.files,
          gigabytes: size.bytes / 1e9,
          toDialogSeconds: (t1 - t0) / 1000,
          offeredMolecules,
          uploadSeconds,
          totalSeconds: (Date.now() - t1) / 1000,
          datasets,
        }),
      });
      expect(
        datasets.map((d) => d.status),
        "every dataset processed",
      ).toEqual(datasets.map(() => "COMPLETE"));
    });
  }
});
