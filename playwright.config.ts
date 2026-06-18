import { defineConfig, devices } from "@playwright/test";
import path from "path";
import dotenv from "dotenv";

// Load test infra env (Postgres + MinIO) before anything else so both the
// dev server (webServer) and globalSetup see it. Only for the upload suite.
if (process.env.E2E_UPLOAD === "1") {
  // override: .env.test is authoritative for the upload suite — it must win over
  // any placeholder DATABASE_URL/S3 env the CI workflow sets for the other jobs.
  dotenv.config({ path: path.join(__dirname, ".env.test"), override: true });
}

/**
 * Playwright config for the MERFISHEYES E2E + performance pipeline.
 *
 * Suites are selected by tag:
 *   @quick   — fast, tiny-dataset smoke run (CI on every PR)
 *   @perf    — performance/timing measurements (writes perf/results/latest.json)
 *   @upload  — upload→save→reload, needs Postgres + MinIO (see docker-compose.test.yml)
 *
 * The dev server is launched with NEXT_PUBLIC_E2E=1 so the app exposes
 * `window.__merfish` timing hooks (see lib/utils/test-hooks.ts).
 */

// Default to port 3100 so the E2E server never clashes with a dev server on 3000.
const PORT = process.env.E2E_PORT || "3100";
const BASE_URL = process.env.E2E_BASE_URL || `http://localhost:${PORT}`;
const IS_CI = !!process.env.CI;
const ENABLE_FIREFOX = process.env.E2E_FIREFOX === "1";

export default defineConfig({
  testDir: path.join(__dirname, "e2e", "tests"),
  // Perf specs live alongside; included via testDir glob below.
  fullyParallel: false, // 3D/WebGL + shared dev server → run serially for stable timings
  workers: 1,
  forbidOnly: IS_CI,
  retries: IS_CI ? 1 : 0,
  timeout: 180_000, // large datasets can take a while to parse in-browser
  expect: { timeout: 30_000 },

  globalSetup: path.join(__dirname, "e2e", "global-setup.ts"),

  reporter: [
    ["list"],
    ["html", { outputFolder: "playwright-report", open: "never" }],
    ["json", { outputFile: "perf/results/playwright-results.json" }],
    [path.join(__dirname, "e2e", "perf", "reporter.ts")],
  ],

  use: {
    baseURL: BASE_URL,
    headless: !process.env.E2E_HEADED,
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
    video: "retain-on-failure",
    actionTimeout: 30_000,
    navigationTimeout: 60_000,
  },

  projects: [
    {
      name: "chromium",
      testDir: path.join(__dirname, "e2e"),
      testMatch: /.*\.spec\.ts/,
      use: { ...devices["Desktop Chrome"] },
    },
    ...(ENABLE_FIREFOX
      ? [
          {
            name: "firefox",
            testDir: path.join(__dirname, "e2e"),
            testMatch: /.*\.spec\.ts/,
            use: { ...devices["Desktop Firefox"] },
          },
        ]
      : []),
  ],

  // Launch the app with E2E hooks enabled. Locally we reuse an already-running
  // dev server only when it was started with NEXT_PUBLIC_E2E=1 (set E2E_REUSE=1).
  webServer: process.env.E2E_NO_SERVER
    ? undefined
    : {
        command: `npm run dev -- --port ${PORT}`,
        url: BASE_URL,
        timeout: 180_000,
        reuseExistingServer: !IS_CI && process.env.E2E_REUSE === "1",
        stdout: "pipe",
        stderr: "pipe",
        env: {
          NEXT_PUBLIC_E2E: "1",
          // Separate build dir so we don't fight a dev server on .next.
          NEXT_DIST_DIR: process.env.NEXT_DIST_DIR || ".next-e2e",
        },
      },
});
