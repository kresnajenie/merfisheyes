import { describe, it, expect, vi } from "vitest";

import {
  applyViewerConfig,
  extractViewerConfig,
  validateViewerConfig,
} from "@/lib/utils/viewer-config";

/**
 * The 2D/3D view mode is part of a single-cell dataset's saved defaults:
 * "Save current view as default" records it and opening the dataset restores
 * it (issue #134). Uploads send { viewMode: "3D" } through the same validator.
 */

const dataset = { id: "ds_1", clusters: [] } as any;

function cellStore(viewMode: "2D" | "3D") {
  return {
    viewMode,
    sceneRotation: 0,
    flipX: false,
    flipY: false,
    selectedColumn: null,
    transformColumn: null,
    sampleTransforms: new Map(),
    setViewMode: vi.fn(),
    setSceneRotation: vi.fn(),
    setFlipX: vi.fn(),
    setFlipY: vi.fn(),
  } as any;
}

describe("viewer config view mode", () => {
  it("saves the cell viewer's 3D view and restores it on load", () => {
    const saved = validateViewerConfig(
      // Through JSON, as it travels to the API and back from the database.
      JSON.parse(JSON.stringify(extractViewerConfig(cellStore("3D"), dataset))),
    );

    expect(saved?.viewMode).toBe("3D");

    const store = cellStore("2D");

    applyViewerConfig(saved, store, dataset);
    expect(store.setViewMode).toHaveBeenCalledWith("3D");
  });

  it("leaves the view mode alone when none was saved", () => {
    const store = cellStore("2D");

    applyViewerConfig({ version: 1, flipX: true }, store, dataset);
    expect(store.setViewMode).not.toHaveBeenCalled();
  });

  it("keeps an upload's 3D default and drops an invalid one", () => {
    expect(validateViewerConfig({ version: 1, viewMode: "3D" })).toEqual({
      version: 1,
      viewMode: "3D",
    });
    expect(validateViewerConfig({ viewMode: "4D" })?.viewMode).toBeUndefined();
  });
});
