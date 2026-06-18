/**
 * Capture thumbnails for BIL dataset entries using Playwright.
 *
 * Prerequisites:
 *   npm install -D playwright
 *   npx playwright install chromium
 *   aws configure  (for S3 upload)
 *
 * Usage:
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex           # single dataset
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex ace-dip-use  # multiple
 *   npx tsx scripts/capture-thumbnails.ts --all                 # all datasets
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex --base-url http://example.com
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex --dry-run  # skip S3 upload
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex --sc-only  # only single cell
 *   npx tsx scripts/capture-thumbnails.ts ace-dud-vex --sm-only  # only single molecule
 */

import { chromium } from "playwright";
import * as fs from "fs";
import * as path from "path";
import { execSync } from "child_process";

// ── Config ──────────────────────────────────────────────────────
const JSON_PATH = path.join(__dirname, "../prisma/bil-examples.json");
const SCREENSHOT_DIR = path.join(__dirname, "../.thumbnails");
const VIEWPORT = { width: 1280, height: 720 };
const WAIT_AFTER_CANVAS_MS = 5000;
const PAGE_TIMEOUT_MS = 60000;

// ── Types ───────────────────────────────────────────────────────
interface Entry {
  label: string;
  datasetType: string;
  s3BaseUrl?: string;
  datasetId?: string;
  thumbnailUrl?: string | null;
}

interface Dataset {
  bilCode: string;
  title: string;
  thumbnailUrl?: string | null;
  entries: Entry[];
  [key: string]: unknown;
}

// ── Helpers ─────────────────────────────────────────────────────

function s3BaseUrlToS3Path(s3BaseUrl: string): string {
  // https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-data2/ace-dud-vex/meyes_output/
  // → s3://merfisheyes-bil/bil-psc-data2/ace-dud-vex/meyes_output
  const match = s3BaseUrl.match(/https?:\/\/([^.]+)\.s3\.[^/]+\/(.+)/);
  if (!match) throw new Error(`Cannot parse s3BaseUrl: ${s3BaseUrl}`);
  const bucket = match[1];
  const prefix = match[2].replace(/\/+$/, "");
  return `s3://${bucket}/${prefix}`;
}

function buildViewerUrl(baseUrl: string, entry: Entry): string {
  const viewerBase =
    entry.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";
  if (entry.s3BaseUrl) {
    const url = entry.s3BaseUrl.replace(/\/+$/, "");
    return `${baseUrl}${viewerBase}/from-s3?url=${encodeURIComponent(url)}`;
  }
  if (entry.datasetId) {
    return `${baseUrl}${viewerBase}/${entry.datasetId}`;
  }
  throw new Error(`Entry "${entry.label}" has no s3BaseUrl or datasetId`);
}

function getThumbnailS3Path(entry: Entry): string {
  if (!entry.s3BaseUrl)
    throw new Error(`Entry "${entry.label}" has no s3BaseUrl`);
  const s3Dir = s3BaseUrlToS3Path(entry.s3BaseUrl);
  return `${s3Dir}/thumbnail.jpg`;
}

// ── Main ────────────────────────────────────────────────────────

async function captureEntry(
  page: any,
  baseUrl: string,
  entry: Entry,
  bilCode: string,
  dryRun: boolean,
  overwrite: boolean,
): Promise<string | null> {
  if (!entry.s3BaseUrl) {
    console.log(`    SKIP "${entry.label}" — no s3BaseUrl`);
    return null;
  }

  if (!overwrite && entry.thumbnailUrl) {
    console.log(`    SKIP "${entry.label}" — already has thumbnail`);
    return null;
  }

  const viewerUrl = buildViewerUrl(baseUrl, entry);
  const typeLabel = entry.datasetType === "single_cell" ? "SC" : "SM";
  console.log(`    [${typeLabel}] ${entry.label}`);
  console.log(`      URL: ${viewerUrl}`);

  try {
    await page.goto(viewerUrl, {
      waitUntil: "domcontentloaded",
      timeout: PAGE_TIMEOUT_MS,
    });

    // Wait for canvas to appear (Three.js renders to canvas)
    console.log(`      Waiting for canvas...`);
    await page.waitForSelector("canvas", { timeout: PAGE_TIMEOUT_MS });

    // Wait additional time for rendering to complete
    console.log(
      `      Canvas found, waiting ${WAIT_AFTER_CANVAS_MS / 1000}s for render...`,
    );
    await page.waitForTimeout(WAIT_AFTER_CANVAS_MS);

    // Zoom in by scrolling on the canvas
    const canvas = page.locator("canvas").first();
    await canvas.scrollIntoViewIfNeeded();
    const box = await canvas.boundingBox();
    if (box) {
      const cx = box.x + box.width / 2;
      const cy = box.y + box.height / 2;
      await page.mouse.move(cx, cy);
      // Scroll down = zoom in (3 ticks)
      await page.mouse.wheel(0, -1000);
      await page.waitForTimeout(1000);
    }

    // Hide all UI elements, keep only the canvas
    await page.evaluate(() => {
      document
        .querySelectorAll(
          "nav, header, footer, [class*='controls'], [class*='panel'], [class*='legend'], [class*='scale'], [class*='navbar'], [class*='toggle'], [class*='button'], [class*='modal'], [class*='toast']",
        )
        .forEach((el) => {
          (el as HTMLElement).style.display = "none";
        });
    });
    await page.waitForTimeout(200);

    // Take screenshot
    const localPath = path.join(
      SCREENSHOT_DIR,
      `${bilCode}_${entry.datasetType}_${entry.label.replace(/[^a-zA-Z0-9._-]/g, "_")}.jpg`,
    );
    await page.screenshot({ path: localPath, type: "jpeg", quality: 85 });
    console.log(`      Screenshot saved: ${localPath}`);

    // Derive S3 and HTTPS paths directly from s3BaseUrl
    const s3Path = getThumbnailS3Path(entry);
    const thumbnailUrl =
      entry.s3BaseUrl!.replace(/\/+$/, "") + "/thumbnail.jpg";

    if (dryRun) {
      console.log(`      DRY RUN — would upload to ${s3Path}`);
    } else {
      console.log(`      Uploading to ${s3Path}...`);
      execSync(
        `aws s3 cp "${localPath}" "${s3Path}" --content-type "image/jpeg"`,
        { stdio: "pipe" },
      );
      console.log(`      Uploaded: ${thumbnailUrl}`);
    }

    return thumbnailUrl;
  } catch (err: any) {
    console.error(`      ERROR: ${err.message}`);
    return null;
  }
}

async function main() {
  const args = process.argv.slice(2);
  if (args.length === 0) {
    console.log(
      "Usage: npx tsx scripts/capture-thumbnails.ts <bilCode...> [--all] [--base-url URL] [--dry-run] [--overwrite] [--sc-only] [--sm-only]",
    );
    process.exit(0);
  }

  // Parse flags
  let baseUrl = "http://localhost:3000";
  let dryRun = false;
  let overwrite = false;
  let processAll = false;
  let scOnly = false;
  let smOnly = false;
  const bilCodes: string[] = [];

  for (let i = 0; i < args.length; i++) {
    if (args[i] === "--base-url" && args[i + 1]) {
      baseUrl = args[++i];
    } else if (args[i] === "--dry-run") {
      dryRun = true;
    } else if (args[i] === "--overwrite") {
      overwrite = true;
    } else if (args[i] === "--all") {
      processAll = true;
    } else if (args[i] === "--sc-only") {
      scOnly = true;
    } else if (args[i] === "--sm-only") {
      smOnly = true;
    } else if (!args[i].startsWith("--")) {
      bilCodes.push(args[i]);
    }
  }

  // Load datasets
  const datasets: Dataset[] = JSON.parse(fs.readFileSync(JSON_PATH, "utf-8"));
  const toProcess = processAll
    ? datasets
    : datasets.filter((d) => bilCodes.includes(d.bilCode));

  if (toProcess.length === 0) {
    console.error("No matching datasets found.");
    process.exit(1);
  }

  // Ensure screenshot dir exists
  fs.mkdirSync(SCREENSHOT_DIR, { recursive: true });

  console.log(`\nCapturing thumbnails for ${toProcess.length} dataset(s)`);
  console.log(`Base URL: ${baseUrl}`);
  console.log(`Dry run: ${dryRun}\n`);

  // Launch browser
  const browser = await chromium.launch({ headless: true });
  const context = await browser.newContext({ viewport: VIEWPORT });
  const page = await context.newPage();

  let updated = 0;

  for (const ds of toProcess) {
    console.log(`\n── ${ds.bilCode}: ${ds.title.slice(0, 60)}... ──`);

    const validEntries = ds.entries.filter(
      (e): e is Entry =>
        typeof e === "object" && e !== null && !!e.label && !!e.s3BaseUrl,
    );

    const scEntries = validEntries.filter(
      (e) => e.datasetType === "single_cell",
    );
    const smEntries = validEntries.filter(
      (e) => e.datasetType === "single_molecule",
    );

    // Capture SC entries
    if (!smOnly) {
      for (const entry of scEntries) {
        const url = await captureEntry(
          page,
          baseUrl,
          entry,
          ds.bilCode,
          dryRun,
          overwrite,
        );
        if (url) {
          entry.thumbnailUrl = url;
          // Parent dataset thumbnail = first SC entry
          if (!ds.thumbnailUrl || overwrite) {
            ds.thumbnailUrl = url;
          }
          updated++;
        }
      }
    }

    // Capture SM entries
    if (!scOnly) {
      for (const entry of smEntries) {
        const url = await captureEntry(
          page,
          baseUrl,
          entry,
          ds.bilCode,
          dryRun,
          overwrite,
        );
        if (url) {
          entry.thumbnailUrl = url;
          updated++;
        }
      }
    }

    // Save after each dataset so progress isn't lost if interrupted
    fs.writeFileSync(JSON_PATH, JSON.stringify(datasets, null, 2));
  }

  await browser.close();
  console.log(`\nDone. Updated ${updated} thumbnail(s) in bil-examples.json`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
