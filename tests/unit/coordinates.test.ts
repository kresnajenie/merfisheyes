import { describe, it, expect } from "vitest";

import {
  normalizeCoordinates,
  // TODO: uncomment after merging debug/sync_colors_and_columns
  // normalizeCoordinatesFlat,
} from "@/lib/utils/coordinates";

describe("normalizeCoordinates (number[][])", () => {
  it("returns null for empty input", () => {
    expect(normalizeCoordinates([])).toBeNull();
  });

  it("normalizes 2D coordinates to [-1, 1]", () => {
    const coords = [
      [0, 0],
      [10, 10],
    ];
    const result = normalizeCoordinates(coords);

    expect(result).not.toBeNull();
    result!.normalized.flat().forEach((v) => {
      expect(v).toBeGreaterThanOrEqual(-1);
      expect(v).toBeLessThanOrEqual(1);
    });
  });

  it("normalizes 3D coordinates to [-1, 1]", () => {
    const coords = [
      [0, 0, 0],
      [10, 10, 10],
    ];
    const result = normalizeCoordinates(coords);

    expect(result).not.toBeNull();
    result!.normalized.flat().forEach((v) => {
      expect(v).toBeGreaterThanOrEqual(-1);
      expect(v).toBeLessThanOrEqual(1);
    });
  });

  it("returns scalingFactor > 0", () => {
    const coords = [
      [0, 0],
      [100, 200],
    ];
    const result = normalizeCoordinates(coords);

    expect(result!.scalingFactor).toBeGreaterThan(0);
  });

  it("handles all-zero coordinates without error", () => {
    const coords = [
      [0, 0],
      [0, 0],
      [0, 0],
    ];
    const result = normalizeCoordinates(coords);

    expect(result).not.toBeNull();
    // scalingFactor should be 1 (fallback for zero range)
    expect(result!.scalingFactor).toBe(1);
  });
});

// TODO: uncomment after merging debug/sync_colors_and_columns
// describe("normalizeCoordinatesFlat (Float32Array)", () => {
//   it("returns null for empty input", () => {
//     expect(normalizeCoordinatesFlat(new Float32Array(0), 2)).toBeNull();
//   });
//
//   it("normalizes 2D flat coordinates to [-1, 1]", () => {
//     const coords = new Float32Array([0, 0, 10, 10]);
//     const result = normalizeCoordinatesFlat(coords, 2);
//
//     expect(result).not.toBeNull();
//     for (let i = 0; i < result!.normalized.length; i++) {
//       expect(result!.normalized[i]).toBeGreaterThanOrEqual(-1);
//       expect(result!.normalized[i]).toBeLessThanOrEqual(1);
//     }
//   });
//
//   it("normalizes 3D flat coordinates to [-1, 1]", () => {
//     const coords = new Float32Array([0, 0, 0, 10, 10, 10]);
//     const result = normalizeCoordinatesFlat(coords, 3);
//
//     expect(result).not.toBeNull();
//     for (let i = 0; i < result!.normalized.length; i++) {
//       expect(result!.normalized[i]).toBeGreaterThanOrEqual(-1);
//       expect(result!.normalized[i]).toBeLessThanOrEqual(1);
//     }
//   });
//
//   it("returns Float32Array output", () => {
//     const coords = new Float32Array([0, 0, 10, 10]);
//     const result = normalizeCoordinatesFlat(coords, 2);
//
//     expect(result!.normalized).toBeInstanceOf(Float32Array);
//   });
//
//   it("preserves length", () => {
//     const coords = new Float32Array([1, 2, 3, 4, 5, 6]);
//     const result = normalizeCoordinatesFlat(coords, 3);
//
//     expect(result!.normalized.length).toBe(coords.length);
//   });
//
//   it("handles all-zero coordinates", () => {
//     const coords = new Float32Array([0, 0, 0, 0]);
//     const result = normalizeCoordinatesFlat(coords, 2);
//
//     expect(result).not.toBeNull();
//     expect(result!.scalingFactor).toBe(1);
//   });
// });
