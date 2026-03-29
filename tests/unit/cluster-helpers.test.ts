// TODO: uncomment after merging debug/sync_colors_and_columns
// These functions (getClusterValue, getCoord, getPointCountFromSpatial) are
// added in the memory optimization refactor and don't exist on develop yet.

import { describe, it, expect } from "vitest";

// import {
//   getClusterValue,
//   getCoord,
//   getPointCountFromSpatial,
// } from "@/lib/StandardizedDataset";

// describe("getClusterValue", () => {
//   it("uses valueIndices when available", () => {
//     const cluster = {
//       column: "leiden",
//       type: "categorical",
//       values: [],
//       valueIndices: new Uint16Array([2, 0, 1, 2]),
//       palette: null,
//       uniqueValues: ["Astrocyte", "Microglia", "Neuron"],
//     };
//
//     expect(getClusterValue(cluster, 0)).toBe("Neuron");
//     expect(getClusterValue(cluster, 1)).toBe("Astrocyte");
//     expect(getClusterValue(cluster, 2)).toBe("Microglia");
//     expect(getClusterValue(cluster, 3)).toBe("Neuron");
//   });
//
//   it("falls back to values array when no indices", () => {
//     const cluster = {
//       column: "leiden",
//       type: "categorical",
//       values: ["Neuron", "Astrocyte", "Microglia"],
//       palette: null,
//     };
//
//     expect(getClusterValue(cluster, 0)).toBe("Neuron");
//     expect(getClusterValue(cluster, 1)).toBe("Astrocyte");
//     expect(getClusterValue(cluster, 2)).toBe("Microglia");
//   });
//
//   it("converts non-string values to string", () => {
//     const cluster = {
//       column: "n_genes",
//       type: "numerical",
//       values: [100, 200, 300],
//       palette: null,
//     };
//
//     expect(getClusterValue(cluster, 0)).toBe("100");
//     expect(getClusterValue(cluster, 1)).toBe("200");
//   });
//
//   it("works with Uint32Array for large unique counts", () => {
//     const cluster = {
//       column: "test",
//       type: "categorical",
//       values: [],
//       valueIndices: new Uint32Array([0, 1]),
//       palette: null,
//       uniqueValues: ["TypeA", "TypeB"],
//     };
//
//     expect(getClusterValue(cluster, 0)).toBe("TypeA");
//     expect(getClusterValue(cluster, 1)).toBe("TypeB");
//   });
// });

// describe("getCoord", () => {
//   it("reads from Float32Array (flat layout)", () => {
//     const spatial = {
//       coordinates: new Float32Array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
//       dimensions: 3,
//       scalingFactor: 1,
//     };
//
//     expect(getCoord(spatial, 0, 0)).toBeCloseTo(1.0);
//     expect(getCoord(spatial, 0, 1)).toBeCloseTo(2.0);
//     expect(getCoord(spatial, 0, 2)).toBeCloseTo(3.0);
//     expect(getCoord(spatial, 1, 0)).toBeCloseTo(4.0);
//     expect(getCoord(spatial, 1, 1)).toBeCloseTo(5.0);
//     expect(getCoord(spatial, 1, 2)).toBeCloseTo(6.0);
//   });
//
//   it("reads from number[][] (nested layout)", () => {
//     const spatial = {
//       coordinates: [
//         [1.0, 2.0],
//         [3.0, 4.0],
//       ],
//       dimensions: 2,
//       scalingFactor: 1,
//     };
//
//     expect(getCoord(spatial, 0, 0)).toBe(1.0);
//     expect(getCoord(spatial, 0, 1)).toBe(2.0);
//     expect(getCoord(spatial, 1, 0)).toBe(3.0);
//     expect(getCoord(spatial, 1, 1)).toBe(4.0);
//   });
//
//   it("handles 2D Float32Array", () => {
//     const spatial = {
//       coordinates: new Float32Array([10, 20, 30, 40]),
//       dimensions: 2,
//       scalingFactor: 1,
//     };
//
//     expect(getCoord(spatial, 0, 0)).toBeCloseTo(10);
//     expect(getCoord(spatial, 0, 1)).toBeCloseTo(20);
//     expect(getCoord(spatial, 1, 0)).toBeCloseTo(30);
//     expect(getCoord(spatial, 1, 1)).toBeCloseTo(40);
//   });
// });

// describe("getPointCountFromSpatial", () => {
//   it("computes count from Float32Array", () => {
//     const spatial = {
//       coordinates: new Float32Array(30),
//       dimensions: 3,
//       scalingFactor: 1,
//     };
//
//     expect(getPointCountFromSpatial(spatial)).toBe(10);
//   });
//
//   it("computes count from Float32Array 2D", () => {
//     const spatial = {
//       coordinates: new Float32Array(20),
//       dimensions: 2,
//       scalingFactor: 1,
//     };
//
//     expect(getPointCountFromSpatial(spatial)).toBe(10);
//   });
//
//   it("computes count from number[][]", () => {
//     const spatial = {
//       coordinates: [
//         [0, 0],
//         [1, 1],
//         [2, 2],
//       ],
//       dimensions: 2,
//       scalingFactor: 1,
//     };
//
//     expect(getPointCountFromSpatial(spatial)).toBe(3);
//   });
//
//   it("returns 0 for empty coordinates", () => {
//     expect(
//       getPointCountFromSpatial({
//         coordinates: new Float32Array(0),
//         dimensions: 2,
//         scalingFactor: 1,
//       }),
//     ).toBe(0);
//
//     expect(
//       getPointCountFromSpatial({
//         coordinates: [],
//         dimensions: 2,
//         scalingFactor: 1,
//       }),
//     ).toBe(0);
//   });
// });

// Placeholder test so this file isn't empty
describe("cluster-helpers (placeholder)", () => {
  it("tests are commented out pending merge", () => {
    expect(true).toBe(true);
  });
});
