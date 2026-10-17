import type { LocalBundle } from "./classify-bundle";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SpatialDataImage } from "@/lib/spatialdata/image";

import { partFileMap, partFiles, storeFileMap } from "./classify-bundle";

import {
  openOmeZarrImageFromFiles,
  openSpatialDataImage,
} from "@/lib/spatialdata/image";
import { loadSpatialDataCells } from "@/lib/spatialdata/loadSpatialDataCells";
import { openSpatialDataStoreFromFiles } from "@/lib/spatialdata/localStore";
import {
  findPointsElements,
  openPointsElement,
} from "@/lib/spatialdata/points";
import { pointsMoleculeDataset } from "@/lib/spatialdata/pointsMolecules";
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
    // What sits beside the store; nothing, when the store is all there is
    const store = openSpatialDataStoreFromFiles(
      storeFileMap(bundle, bundle.cells),
      bundle.cells.root ? bundle.files : undefined,
    );
    const name = (bundle.cells.root || bundle.name).replace(/\.zarr$/, "");

    // A store may hold only an image; its table is what makes it cells
    if (elementNames(store.keys, "tables").length > 0) {
      opened.cells = await loadSpatialDataCells(store, id, name, (p, m) =>
        onProgress(p * 0.6, m),
      );
    }
    opened.image = (await openSpatialDataImage(store)) ?? undefined;
    // Molecules stored in the store itself (the Zarr points format), unless
    // a molecule folder was dropped beside it
    const [pointsName] = await findPointsElements(store);

    if (pointsName && opened.cells && bundle.molecules === undefined) {
      onProgress(65, "Opening molecules...");
      opened.molecules = await pointsMoleculeDataset(
        await openPointsElement(store, pointsName),
        `${id}_${pointsName}`,
        `${name} — ${pointsName}`,
      );
    }
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
