import type { LocalBundle } from "./classify-bundle";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SpatialDataImage } from "@/lib/spatialdata/image";

import { partFileMap, partFiles } from "./classify-bundle";

import {
  openOmeZarrImageFromFiles,
  openSpatialDataImage,
} from "@/lib/spatialdata/image";
import { loadSpatialDataCells } from "@/lib/spatialdata/loadSpatialDataCells";
import { openSpatialDataStoreFromFiles } from "@/lib/spatialdata/localStore";
import { elementNames } from "@/lib/spatialdata/store";

export interface OpenedBundle {
  cells?: StandardizedDataset;
  molecules?: SingleMoleculeDataset;
  image?: SpatialDataImage;
}

/**
 * Open the parts of a dropped folder for viewing. Only metadata is read here
 * (and the cells' coordinates); expression, molecules and image tiles are
 * read from the files later, as the viewer asks for them.
 */
export async function openLocalBundle(
  bundle: LocalBundle,
  onProgress: (progress: number, message: string) => void,
): Promise<OpenedBundle> {
  const opened: OpenedBundle = {};
  const id = `local_${Date.now()}`;

  if (bundle.cells?.kind === "chunked") {
    const { StandardizedDataset } = await import("@/lib/StandardizedDataset");

    opened.cells = await StandardizedDataset.fromLocalChunked(
      partFiles(bundle, bundle.cells.root),
      (p, m) => onProgress(p * 0.6, m),
    );
  } else if (bundle.cells) {
    const store = openSpatialDataStoreFromFiles(
      bundle.files,
      bundle.cells.root,
    );
    const name = (bundle.cells.root || bundle.name).replace(/\.zarr$/, "");

    // A store may hold only an image; its table is what makes it cells
    if (elementNames(store.keys, "tables").length > 0) {
      opened.cells = await loadSpatialDataCells(store, id, name, (p, m) =>
        onProgress(p * 0.6, m),
      );
    }
    opened.image = (await openSpatialDataImage(store)) ?? undefined;
  }

  if (bundle.molecules !== undefined) {
    onProgress(70, "Opening molecules...");
    const { SingleMoleculeDataset } = await import(
      "@/lib/SingleMoleculeDataset"
    );

    opened.molecules = await SingleMoleculeDataset.fromLocalChunked(
      partFiles(bundle, bundle.molecules),
    );
  }

  // A standalone image; a store's own image, if it has one, came with it
  if (bundle.image !== undefined && !opened.image) {
    onProgress(85, "Opening image...");
    opened.image = await openOmeZarrImageFromFiles(
      partFileMap(bundle, bundle.image),
      bundle.image || bundle.name,
    );
  }
  onProgress(100, "Ready");

  return opened;
}
