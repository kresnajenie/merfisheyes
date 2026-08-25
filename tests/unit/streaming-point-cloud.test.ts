import { describe, it, expect } from "vitest";
import * as THREE from "three";

import {
  appendToSmPointCloud,
  createStreamingSmPointCloud,
  streamedMoleculeCount,
} from "@/lib/webgl/point-cloud";

describe("streaming point cloud", () => {
  it("starts empty and reveals molecules as batches arrive", () => {
    const cloud = createStreamingSmPointCloud(6, "#ff0000", 1);

    expect(streamedMoleculeCount(cloud)).toBe(0);
    expect(cloud.geometry.drawRange.count).toBe(0);

    appendToSmPointCloud(cloud, Float32Array.from([0, 0, 0, 1, 2, 3]));
    expect(streamedMoleculeCount(cloud)).toBe(2);
    expect(cloud.geometry.drawRange.count).toBe(2);

    appendToSmPointCloud(cloud, Float32Array.from([4, 5, 6]));
    expect(streamedMoleculeCount(cloud)).toBe(3);

    const positions = cloud.geometry.getAttribute("position").array as Float32Array;

    expect(Array.from(positions.subarray(0, 9))).toEqual([0, 0, 0, 1, 2, 3, 4, 5, 6]);
    // The unwritten tail stays zeroed and unrendered.
    expect(Array.from(positions.subarray(9))).toEqual([0, 0, 0, 0, 0, 0, 0, 0, 0]);
  });

  it("keeps bounds over the filled molecules only", () => {
    const cloud = createStreamingSmPointCloud(1000, "#00ff00", 1);

    appendToSmPointCloud(cloud, Float32Array.from([10, 20, 30, -5, 0, 7]));

    const box = cloud.geometry.boundingBox as THREE.Box3;

    expect([box.min.x, box.min.y, box.min.z]).toEqual([-5, 0, 7]);
    expect([box.max.x, box.max.y, box.max.z]).toEqual([10, 20, 30]);
    // The 998 unwritten molecules (all zeros) must not drag the box to origin.
    expect(box.min.y).toBe(0);
    expect(box.min.x).not.toBe(0);
  });

  it("never writes past its capacity", () => {
    const cloud = createStreamingSmPointCloud(2, "#0000ff", 1);
    const total = appendToSmPointCloud(
      cloud,
      Float32Array.from([1, 1, 1, 2, 2, 2, 3, 3, 3]),
    );

    expect(total).toBe(2);
    expect(cloud.geometry.drawRange.count).toBe(2);
  });

  it("colors every vertex with the cloud color", () => {
    const cloud = createStreamingSmPointCloud(4, "#ff0000", 1);
    const colors = cloud.geometry.getAttribute("color").array as Float32Array;

    expect(Array.from(colors)).toEqual([1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0]);
  });
});
