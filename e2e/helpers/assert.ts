import { expect, type TestInfo } from "@playwright/test";

/**
 * Rich, localized failure reporting. Every assertion carries the dataset name,
 * format, the step that failed, and expected vs actual, so a CI log line is
 * self-explanatory. Screenshots/traces are captured automatically on failure
 * via playwright.config.ts.
 */

export interface FailureContext {
  dataset: string;
  format: string;
  step: string;
  /** Optional timing context for perf-related assertions. */
  timingMs?: number;
}

function header(ctx: FailureContext, expected: unknown, actual: unknown): string {
  const lines = [
    `✗ E2E assertion failed`,
    `  dataset:  ${ctx.dataset}`,
    `  format:   ${ctx.format}`,
    `  step:     ${ctx.step}`,
    `  expected: ${JSON.stringify(expected)}`,
    `  actual:   ${JSON.stringify(actual)}`,
  ];
  if (ctx.timingMs != null) lines.push(`  timing:   ${ctx.timingMs.toFixed(1)} ms`);
  return lines.join("\n");
}

/** Tag the current test with dataset/format so reports group cleanly. */
export function annotate(testInfo: TestInfo, dataset: string, format: string): void {
  testInfo.annotations.push({ type: "dataset", description: dataset });
  testInfo.annotations.push({ type: "format", description: format });
}

export function assertEqual<T>(actual: T, expected: T, ctx: FailureContext): void {
  expect(actual, header(ctx, expected, actual)).toEqual(expected);
}

export function assertTrue(actual: boolean, ctx: FailureContext): void {
  expect(actual, header(ctx, true, actual)).toBe(true);
}

/** Assert a numeric count is within a tolerance fraction (e.g. 0 for exact). */
export function assertCount(
  actual: number,
  expected: number,
  ctx: FailureContext,
  tolerance = 0,
): void {
  const lo = expected * (1 - tolerance);
  const hi = expected * (1 + tolerance);
  expect(
    actual >= lo && actual <= hi,
    header({ ...ctx }, `${expected} (±${tolerance * 100}%)`, actual),
  ).toBe(true);
}
