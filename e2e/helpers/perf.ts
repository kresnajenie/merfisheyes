import type { TestInfo } from "@playwright/test";

/**
 * Performance metric capture. Each metric is attached to the test as a JSON
 * attachment named "perf"; the custom reporter (e2e/perf/reporter.ts) collects
 * all of them into perf/results/latest.json for baseline comparison.
 */

export interface MetricRecord {
  datasetId: string;
  format: string;
  /** e.g. "loadToRender", "geneSelectionLatency", "uploadTime". */
  metric: string;
  /** Median (or single) value, milliseconds. */
  ms: number;
  /** All raw samples when measured repeatedly. */
  samples?: number[];
  tier?: string;
}

export function median(values: number[]): number {
  if (!values.length) return NaN;
  const s = [...values].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
}

/** Record a single metric for the current test. */
export async function recordMetric(
  testInfo: TestInfo,
  record: MetricRecord,
): Promise<void> {
  await testInfo.attach("perf", {
    body: JSON.stringify(record),
    contentType: "application/json",
  });
}

/** Record a metric from raw samples, using the median as the headline value. */
export async function recordSamples(
  testInfo: TestInfo,
  base: Omit<MetricRecord, "ms" | "samples">,
  samples: number[],
): Promise<void> {
  await recordMetric(testInfo, { ...base, ms: median(samples), samples });
}
