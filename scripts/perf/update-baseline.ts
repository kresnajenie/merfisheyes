/**
 * Promote the latest perf run to the committed baseline for this platform.
 * This is the intentional "bless" step — run it when a timing change is
 * expected/accepted, then commit the updated baseline file.
 *
 *   npm run test:perf:update
 *
 * Merges into any existing baseline (keeps metrics not present in this run).
 */

import fs from "node:fs";
import path from "node:path";

import { recordKey, baselinePath } from "./compare-baseline";

const ROOT = path.resolve(__dirname, "..", "..");
const RESULTS = path.join(ROOT, "perf", "results", "latest.json");

interface PerfRecord {
  datasetId: string;
  metric: string;
  ms: number;
  browser?: string;
}
interface ResultsFile {
  platform: string;
  gitSha: string;
  records: PerfRecord[];
}

function main(): void {
  if (!fs.existsSync(RESULTS)) {
    console.error(`No perf results at ${path.relative(ROOT, RESULTS)}. Run \`npm run test:perf\` first.`);
    process.exit(2);
  }
  const results: ResultsFile = JSON.parse(fs.readFileSync(RESULTS, "utf8"));
  const blPath = baselinePath(results.platform);

  const existing: Record<string, number> = fs.existsSync(blPath)
    ? JSON.parse(fs.readFileSync(blPath, "utf8")).metrics ?? {}
    : {};

  for (const r of results.records) existing[recordKey(r)] = Number(r.ms.toFixed(2));

  fs.mkdirSync(path.dirname(blPath), { recursive: true });
  const out = {
    platform: results.platform,
    updatedAt: new Date().toISOString(),
    gitSha: results.gitSha,
    metrics: existing,
  };
  fs.writeFileSync(blPath, JSON.stringify(out, null, 2));
  console.log(
    `Baseline updated: ${path.relative(ROOT, blPath)} (${results.records.length} metric(s) from ${results.gitSha}).`,
  );
  console.log("Commit this file to record the new baseline.");
}

main();
