/**
 * Compare the latest perf run against the committed baseline and gate on
 * regressions.
 *
 *   npm run test:perf:compare                 # default thresholds
 *   npm run test:perf:compare -- --threshold 25 --warn 12
 *   PERF_FAIL_PCT=30 npm run test:perf:compare
 *
 * Exit code is non-zero if any metric regressed beyond the fail threshold, so
 * CI fails the build. Improvements are reported, never failed.
 */

import fs from "node:fs";
import path from "node:path";

const ROOT = path.resolve(__dirname, "..", "..");
const RESULTS = path.join(ROOT, "perf", "results", "latest.json");
const BASELINE_DIR = path.join(ROOT, "perf", "baselines");

interface PerfRecord {
  datasetId: string;
  format: string;
  metric: string;
  ms: number;
  browser?: string;
  tier?: string;
}
interface ResultsFile {
  platform: string;
  gitSha: string;
  records: PerfRecord[];
}
interface BaselineFile {
  platform: string;
  updatedAt: string;
  gitSha: string;
  metrics: Record<string, number>;
}

function arg(flag: string): string | undefined {
  const i = process.argv.indexOf(flag);
  return i >= 0 ? process.argv[i + 1] : undefined;
}

const FAIL_PCT = Number(arg("--threshold") ?? process.env.PERF_FAIL_PCT ?? 20);
const WARN_PCT = Number(arg("--warn") ?? process.env.PERF_WARN_PCT ?? 10);

export function recordKey(r: { browser?: string; datasetId: string; metric: string }): string {
  return `${r.browser ?? "unknown"}/${r.datasetId}/${r.metric}`;
}

export function baselinePath(platform: string): string {
  return path.join(BASELINE_DIR, `${platform}.baseline.json`);
}

export type PerfStatus = "OK" | "WARN" | "FAIL" | "FAST";

/** Percent change of `cur` relative to `base` (positive = slower). */
export function pctChange(base: number, cur: number): number {
  return ((cur - base) / base) * 100;
}

/** Classify a percent change against warn/fail thresholds. */
export function classify(pct: number, warnPct: number, failPct: number): PerfStatus {
  if (pct > failPct) return "FAIL";
  if (pct > warnPct) return "WARN";
  if (pct < -warnPct) return "FAST";
  return "OK";
}

function pad(s: string, n: number): string {
  return s.length >= n ? s : s + " ".repeat(n - s.length);
}
function padL(s: string, n: number): string {
  return s.length >= n ? s : " ".repeat(n - s.length) + s;
}

function main(): void {
  if (!fs.existsSync(RESULTS)) {
    console.error(`No perf results at ${path.relative(ROOT, RESULTS)}. Run \`npm run test:perf\` first.`);
    process.exit(2);
  }
  const results: ResultsFile = JSON.parse(fs.readFileSync(RESULTS, "utf8"));
  const blPath = baselinePath(results.platform);

  if (!fs.existsSync(blPath)) {
    console.log(`No baseline for platform '${results.platform}' at ${path.relative(ROOT, blPath)}.`);
    console.log("Nothing to compare. Bless the current run with `npm run test:perf:update`.");
    process.exit(0);
  }
  const baseline: BaselineFile = JSON.parse(fs.readFileSync(blPath, "utf8"));

  const rows: {
    key: string;
    base: number | null;
    cur: number;
    pct: number | null;
    status: "OK" | "WARN" | "FAIL" | "NEW" | "FAST";
  }[] = [];
  let failed = 0;
  let warned = 0;

  for (const r of results.records) {
    const key = recordKey(r);
    const base = baseline.metrics[key];
    if (base == null) {
      rows.push({ key, base: null, cur: r.ms, pct: null, status: "NEW" });
      continue;
    }
    const pct = pctChange(base, r.ms);
    const status = classify(pct, WARN_PCT, FAIL_PCT);
    if (status === "FAIL") failed++;
    else if (status === "WARN") warned++;
    rows.push({ key, base, cur: r.ms, pct, status });
  }

  // Report.
  console.log(`\nPerf comparison — platform=${results.platform}  baseline=${baseline.gitSha}  current=${results.gitSha}`);
  console.log(`Thresholds: warn >${WARN_PCT}%, fail >${FAIL_PCT}%\n`);
  console.log(
    `${pad("metric (browser/dataset/metric)", 48)} ${padL("base ms", 10)} ${padL("cur ms", 10)} ${padL("Δ%", 8)}  status`,
  );
  console.log("-".repeat(90));
  for (const row of rows.sort((a, b) => (b.pct ?? -999) - (a.pct ?? -999))) {
    console.log(
      `${pad(row.key, 48)} ${padL(row.base?.toFixed(1) ?? "—", 10)} ${padL(row.cur.toFixed(1), 10)} ` +
        `${padL(row.pct == null ? "—" : `${row.pct >= 0 ? "+" : ""}${row.pct.toFixed(1)}`, 8)}  ${row.status}`,
    );
  }
  console.log("-".repeat(90));
  console.log(`\n${rows.length} metric(s): ${failed} regressed, ${warned} warning(s).`);

  if (failed > 0) {
    console.error(`\n✗ Performance regression: ${failed} metric(s) slower than baseline by >${FAIL_PCT}%.`);
    process.exit(1);
  }
  console.log("\n✓ No performance regressions beyond threshold.");
}

// Run only when executed directly (tsx sets argv[1] to this file), not when
// imported by update-baseline.ts or the unit tests. Avoids referencing
// `require`/`module`, which aren't defined under ESM/Vitest.
if (process.argv[1] && /compare-baseline\.ts$/.test(process.argv[1])) {
  main();
}
