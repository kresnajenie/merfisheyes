import { describe, expect, it } from "vitest";

import { calculateGenePercentile } from "@/lib/webgl/visualization-utils";

/** The sort it replaces, kept as the reference. */
function bySorting(arr: number[], p: number) {
  const s = arr.filter((v) => !isNaN(v) && v !== null).sort((a, b) => a - b);

  if (s.length === 0) return NaN;

  return s[Math.min(Math.max(0, Math.floor(s.length * p)), s.length - 1)];
}

describe("calculateGenePercentile", () => {
  it("matches sorting, on sparse data with ties and NaNs", () => {
    let seed = 7;
    const rand = () => (seed = (seed * 16807) % 2147483647) / 2147483647;

    for (let trial = 0; trial < 200; trial++) {
      const n = 1 + Math.floor(rand() * 3000);
      const arr = Array.from({ length: n }, () => {
        const r = rand();

        // mostly zeros, some repeats, a few NaNs — like a gene over cells
        return r < 0.6 ? 0 : r < 0.65 ? NaN : Math.round(rand() * 50) / 10;
      });

      for (const p of [0, 0.5, 0.95, 0.99, 1]) {
        const want = bySorting(arr, p);
        const got = calculateGenePercentile(arr, p);

        if (isNaN(want)) expect(got).toBeNaN();
        else expect(got).toBe(want);
      }
    }
  });

  it("is NaN when there are no numbers", () => {
    expect(calculateGenePercentile([NaN, NaN], 0.95)).toBeNaN();
    expect(calculateGenePercentile([], 0.95)).toBeNaN();
  });
});
