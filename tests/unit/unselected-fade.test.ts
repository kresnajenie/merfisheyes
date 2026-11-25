/**
 * Cells outside the selected cell types fade as the view zooms out, so the
 * selected ones stand out against the tissue. The colouring marks them with
 * a negative alpha (the same alpha, signed), and the cell shader fades those
 * by how far the camera is from the view that shows the whole dataset.
 */
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { describe, expect, it } from "vitest";

import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  fragmentShader,
  moleculeVertexShader,
  vertexShader,
} from "@/lib/webgl/shaders";
import { updateCelltypeVisualization } from "@/lib/webgl/visualization-utils";

const dataset = {
  getPointCount: () => 4,
  clusters: [
    {
      column: "celltype",
      type: "categorical",
      values: ["a", "b", "a", "c"],
      palette: { a: "#ff0000", b: "#00ff00", c: "#0000ff" },
    },
  ],
} as unknown as StandardizedDataset;

describe("cells outside the selected cell types", () => {
  it("carry their alpha negated; selected cells do not", () => {
    const result = updateCelltypeVisualization(
      dataset,
      "celltype",
      new Set(["a"]),
      {},
      1,
      1,
    )!;
    const grey = VISUALIZATION_CONFIG.GREYED_OUT_ALPHA;

    expect(Array.from(result.alphas)).toEqual([1, -grey, 1, -grey]);
  });

  it("nothing is marked when no cell type is selected", () => {
    const result = updateCelltypeVisualization(
      dataset,
      "celltype",
      new Set(),
      {},
      1,
      1,
    )!;

    expect(Array.from(result.alphas).every((a) => a > 0)).toBe(true);
  });
});

describe("the cell shader", () => {
  it("fades a marked cell by zoom, and draws every cell with a positive alpha", () => {
    expect(vertexShader).toContain("vAlpha = abs(alpha);");
    expect(vertexShader).toContain("if (alpha < 0.0 && uHomeDistance > 0.0)");
    expect(vertexShader).toContain("uHomeDistance / -mvPosition.z");
    expect(fragmentShader).toContain("vAlpha * uOpacity");
  });

  it("leaves no placeholder in either shader, and molecules do not fade", () => {
    expect(vertexShader).not.toContain("POINT_FADE");
    expect(moleculeVertexShader).not.toContain("POINT_FADE");
    expect(moleculeVertexShader).not.toContain("uHomeDistance");
  });
});
