/**
 * The default colouring of a table with no cluster column (a Visium HD
 * table holds only array_row, array_col, in_tissue besides its keys):
 * never a positional column — coloured by array_col a capture area is
 * stripes — but one neutral colour.
 */
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { describe, expect, it } from "vitest";

import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  selectBestClusterColumn,
  selectBestClusterColumnByName,
} from "@/lib/utils/dataset-utils";
import { updateCelltypeVisualization } from "@/lib/webgl/visualization-utils";

const categorical = (names: string[]) =>
  Object.fromEntries(names.map((n) => [n, "categorical"]));

describe("the default colouring column", () => {
  it("skips positional columns, and is none when only they are left", () => {
    const visiumHd = ["array_col", "array_row", "in_tissue"];

    expect(selectBestClusterColumnByName(visiumHd, categorical(visiumHd))).toBeNull();
    const withClusters = ["array_col", "in_tissue", "clusters"];

    expect(selectBestClusterColumnByName(withClusters, categorical(withClusters))).toBe("clusters");
    // A positional column can still be chosen by hand: it is only not a default
    expect(selectBestClusterColumnByName(["leiden", "array_col"], categorical(["leiden"]))).toBe("leiden");
  });

  it("does the same for a loaded dataset", () => {
    const ds = (columns: string[]) =>
      ({ clusters: columns.map((column) => ({ column, type: "categorical" })) }) as unknown as StandardizedDataset;

    expect(selectBestClusterColumn(ds(["array_row", "array_col"]))).toBeNull();
    expect(selectBestClusterColumn(ds(["array_row", "sample"]))).toBe("sample");
  });
});

describe("cells with no column to colour by", () => {
  it("are drawn in one translucent neutral colour, not left black", () => {
    const dataset = { getPointCount: () => 3, clusters: null } as unknown as StandardizedDataset;
    const result = updateCelltypeVisualization(dataset, null, new Set(), {}, 1, 1);

    expect(result).not.toBeNull();
    const hex = VISUALIZATION_CONFIG.NEUTRAL_POINT_COLOR;
    const expected = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);

    for (let p = 0; p < 3; p++) {
      for (let c = 0; c < 3; c++) {
        expect(result!.colors[p * 3 + c]).toBeCloseTo(expected[c], 5);
      }
      expect(result!.alphas[p]).toBeCloseTo(VISUALIZATION_CONFIG.NEUTRAL_POINT_ALPHA, 6);
      expect(result!.sizes[p]).toBeGreaterThan(0);
    }
  });
});
