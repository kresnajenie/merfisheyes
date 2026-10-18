import { describe, expect, it } from "vitest";

import { PickingGrid } from "@/lib/webgl/picking-grid";

/** Brute force, for the grid to agree with. */
function nearestSlow(
  positions: Float32Array,
  count: number,
  x: number,
  y: number,
  radius: number,
): number | null {
  let best: number | null = null;
  let bestD = radius;

  for (let i = 0; i < count; i++) {
    const d = Math.hypot(positions[i * 3] - x, positions[i * 3 + 1] - y);

    if (d <= bestD) {
      bestD = d;
      best = i;
    }
  }

  return best;
}

describe("PickingGrid", () => {
  it("finds the same point as a pass over every point, or none", () => {
    const count = 5000;
    const positions = new Float32Array(count * 3);
    let seed = 7;
    const rand = () => (seed = (seed * 16807) % 2147483647) / 2147483647;

    for (let i = 0; i < count; i++) {
      positions[i * 3] = rand() * 1000 - 500;
      positions[i * 3 + 1] = rand() * 200;
      positions[i * 3 + 2] = rand();
    }
    const grid = new PickingGrid(positions, count);

    for (let q = 0; q < 300; q++) {
      const x = rand() * 1100 - 550;
      const y = rand() * 220 - 10;
      const radius = rand() * 5;

      expect(grid.nearest(x, y, radius)).toBe(
        nearestSlow(positions, count, x, y, radius),
      );
    }
  });

  it("takes a distance of its own, for picking in screen space", () => {
    const positions = Float32Array.from([0, 0, 0, 10, 0, 0, 20, 0, 0]);
    const grid = new PickingGrid(positions, 3);

    // The planar nearest to (9, 0) is point 1; a screen-space distance
    // that favours point 2 wins within the cells searched
    expect(grid.nearest(9, 0, 5)).toBe(1);
    expect(grid.nearest(9, 0, 15, (i) => (i === 2 ? 1 : 3))).toBe(2);
  });

  it("copes with no points and with points on one spot", () => {
    expect(new PickingGrid(new Float32Array(0), 0).nearest(0, 0, 1)).toBeNull();
    const same = Float32Array.from([1, 1, 0, 1, 1, 0]);

    expect([0, 1]).toContain(new PickingGrid(same, 2).nearest(1, 1, 0.1));
  });
});
