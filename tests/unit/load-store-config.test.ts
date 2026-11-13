/**
 * A SpatialData store opened by dataset id takes its owner's saved defaults,
 * as any cells do (#206): a server upload is a store, and its "open in 3D"
 * and "Save current view as default" were read from the record and dropped.
 */
import type { SpatialDataStore } from "@/lib/spatialdata/store";
import type { ViewerConfig } from "@/lib/utils/viewer-config";

import { readdirSync, statSync } from "node:fs";
import { readFile } from "node:fs/promises";
import path from "node:path";

import { afterEach, describe, expect, it, vi } from "vitest";

const MAIN = path.resolve(__dirname, "../data/tiny/zarr/spatialdata_v3.zarr");

function listFiles(dir: string, prefix = ""): string[] {
  return readdirSync(dir).flatMap((name) => {
    const full = path.join(dir, name);
    const key = prefix ? `${prefix}/${name}` : name;

    return statSync(full).isDirectory() ? listFiles(full, key) : [key];
  });
}

/** The tiny store's table only, as the app's bucket would list it. */
function tinyStore(): SpatialDataStore {
  return {
    keys: listFiles(MAIN).filter((k) => !k.startsWith("images/")),
    store: {
      get: (key: string) =>
        readFile(path.join(MAIN, key.replace(/^\/+/, ""))).then(
          (b) => new Uint8Array(b),
          () => undefined,
        ),
    },
    sidecar: { get: async () => undefined },
  } as SpatialDataStore;
}

vi.mock("@/lib/spatialdata/store", async (importOriginal) => ({
  ...(await importOriginal<typeof import("@/lib/spatialdata/store")>()),
  openSpatialDataStore: vi.fn(async () => tinyStore()),
}));

const { loadDataset } = await import("@/lib/viewer/load-dataset");
const { createVisualizationStoreInstance } = await import(
  "@/lib/stores/createVisualizationStore"
);
const { createSingleMoleculeVisualizationStoreInstance } = await import(
  "@/lib/stores/createSingleMoleculeVisualizationStore"
);

/** Open the store as `ds_store`, whose record carries `viewerConfig`. */
async function open(viewerConfig: ViewerConfig | null) {
  vi.stubGlobal(
    "fetch",
    vi.fn(async () =>
      Response.json({
        title: "store",
        datasetType: "spatialdata",
        ownerId: "u1",
        viewerConfig,
      }),
    ),
  );
  const viz = createVisualizationStoreInstance();
  const loaded = await loadDataset(
    { kind: "spatialdata", datasetId: "ds_store" },
    {
      panel: "left",
      viz: viz.getState(),
      smViz: createSingleMoleculeVisualizationStoreInstance().getState(),
      findCells: () => undefined,
      findMolecules: () => undefined,
      onProgress: () => {},
      isCancelled: () => false,
    },
  );

  return { loaded: loaded!, viz: viz.getState() };
}

afterEach(() => {
  vi.unstubAllGlobals();
});

describe("a store opened by dataset id", () => {
  it("opens in the owner's saved view: 3D, column, rotation", async () => {
    const { loaded, viz } = await open({
      viewMode: "3D",
      priorityColumn: "area",
      sceneRotation: 90,
    });

    expect(loaded.cells).not.toBeNull();
    expect(loaded.registration.viewerConfig?.viewMode).toBe("3D");
    expect(viz.viewMode).toBe("3D");
    expect(viz.selectedColumn).toBe("area");
    expect(viz.sceneRotation).toBe(90);
  });

  it("with nothing saved, opens as before: 2D and the best column", async () => {
    const { viz } = await open(null);

    expect(viz.viewMode).toBe("2D");
    expect(viz.selectedColumn).toBe("celltype");
  });
});
