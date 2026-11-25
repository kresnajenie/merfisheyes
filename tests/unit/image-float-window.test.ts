/**
 * A float32 channel (MIBI, IMC counts) has no fixed range: its display
 * window and its sliders' range are read off its coarsest level (#218).
 */
import { describe, expect, it } from "vitest";

import { floatWindow } from "@/lib/spatialdata/ImageLayer";

describe("floatWindow", () => {
  it("windows 0 to the percentile of the non-zero pixels, on sliders that reach the largest", () => {
    // 1000 values 0.01 … 10, and as many zeros (background)
    const data = new Float32Array(2000);

    for (let i = 0; i < 1000; i++) data[i * 2] = (i + 1) / 100;

    expect(floatWindow(data, 0.9)).toEqual({
      min: 0,
      max: 9.01,
      range: [0, 10],
    });
  });

  it("keeps four digits, as the values go in links", () => {
    const { max, range } = floatWindow(new Float32Array([1 / 3, 2 / 3]), 0.5);

    expect(max).toBe(0.6667);
    expect(range).toEqual([0, 0.6667]);
  });

  it("an empty or all-zero plane gets a window that can still be drawn", () => {
    const blank = { min: 0, max: 1, range: [0, 1] };

    expect(floatWindow(new Float32Array(16), 0.999)).toEqual(blank);
    expect(floatWindow(new Float32Array([NaN, 0, Infinity]), 0.999)).toEqual(
      blank,
    );
  });
});
