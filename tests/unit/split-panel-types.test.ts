import { describe, expect, it } from "vitest";

import {
  panelTypeForDataset,
  panelTypeForViewer,
} from "@/lib/stores/splitScreenStore";

/**
 * The split view's right panel opens what the dataset's own viewer would:
 * a catalog entry by its type (a SpatialData-DB card is "spatialdata"), a
 * pasted link by its viewer.
 */
describe("split panel types", () => {
  it.each([
    ["single_cell", "cell"],
    ["h5ad", "cell"],
    [null, "cell"],
    ["single_molecule", "sm"],
    ["spatialdata", "spatialdata"],
    ["image", "image"],
  ])("a %s dataset opens in a %s panel", (type, panel) => {
    expect(panelTypeForDataset(type)).toBe(panel);
  });

  it.each([
    ["viewer", "cell"],
    ["sm-viewer", "sm"],
    ["lm-viewer", "lm"],
    ["spatialdata-viewer", "spatialdata"],
    ["image-viewer", "image"],
  ])("a /%s link opens in a %s panel", (segment, panel) => {
    expect(panelTypeForViewer(segment)).toBe(panel);
  });

  it("is nothing for a page that is not a viewer", () => {
    expect(panelTypeForViewer("explore")).toBeNull();
  });
});
