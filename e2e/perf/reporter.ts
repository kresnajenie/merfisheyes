import { execSync } from "node:child_process";
import fs from "node:fs";
import path from "node:path";

import type {
  Reporter,
  TestCase,
  TestResult,
} from "@playwright/test/reporter";

/**
 * Collects "perf" JSON attachments emitted by tests (see e2e/helpers/perf.ts)
 * and writes them to perf/results/latest.json with environment metadata. The
 * baseline comparison reads this file (scripts/perf/compare-baseline.ts).
 *
 * The env key (`platform-browser`) namespaces baselines so timings from
 * different machines/browsers are never compared against each other.
 */

interface PerfRecord {
  datasetId: string;
  format: string;
  metric: string;
  ms: number;
  samples?: number[];
  tier?: string;
  browser?: string;
}

const RESULTS_DIR = path.resolve(__dirname, "..", "..", "perf", "results");
const RESULTS_FILE = path.join(RESULTS_DIR, "latest.json");

function gitSha(): string {
  try {
    return execSync("git rev-parse --short HEAD").toString().trim();
  } catch {
    return "unknown";
  }
}

export default class PerfReporter implements Reporter {
  private records: PerfRecord[] = [];

  onTestEnd(test: TestCase, result: TestResult): void {
    const browser = test.parent.project()?.name ?? "unknown";
    for (const att of result.attachments) {
      if (att.name !== "perf" || !att.body) continue;
      try {
        const rec = JSON.parse(att.body.toString("utf8")) as PerfRecord;
        this.records.push({ ...rec, browser });
      } catch {
        // ignore malformed perf attachment
      }
    }
  }

  onEnd(): void {
    if (!this.records.length) return;
    fs.mkdirSync(RESULTS_DIR, { recursive: true });
    const out = {
      platform: process.platform,
      gitSha: gitSha(),
      generatedAt: new Date().toISOString(),
      records: this.records,
    };
    fs.writeFileSync(RESULTS_FILE, JSON.stringify(out, null, 2));
    // eslint-disable-next-line no-console
    console.log(
      `\n[perf] wrote ${this.records.length} metric(s) -> ${path.relative(process.cwd(), RESULTS_FILE)}`,
    );
  }
}
