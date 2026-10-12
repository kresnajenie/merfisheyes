import fs from "node:fs";
import os from "node:os";
import path from "node:path";

import { REPO_ROOT } from "./databank";

const TINY = path.join(REPO_ROOT, "tests", "data", "tiny");

/** Which already-processed parts to put in a folder, and in what form. */
export interface BundleParts {
  /**
   * "chunked": the merfisheyes cell folder that keeps raw coordinates;
   * "normalized": the older one that does not; "spatialdata": the tiny store
   * (cells and an image inside it); "h5ad": raw cells that need parsing.
   */
  cells?: "chunked" | "normalized" | "spatialdata" | "h5ad";
  molecules?: boolean;
  /** Standalone OME-Zarr images, by folder name. */
  images?: string[];
}

/**
 * Build a folder of already-processed parts in a temp dir, the way someone
 * would have it on disk, out of the committed tiny fixtures. Returns the path
 * to hand to a folder input. The standalone image is stretched to the tiny
 * cells' extent (about -1000..1000) so it sits under them.
 */
export function makeBundle(parts: BundleParts, name = "sample"): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "bundle-")),
    name,
  );
  const copy = (from: string, to: string) =>
    fs.cpSync(path.join(TINY, from), path.join(folder, to), {
      recursive: true,
    });

  fs.mkdirSync(folder);
  if (parts.cells === "chunked") copy("chunked/raw", "cells");
  if (parts.cells === "normalized") copy("chunked/tiny", "cells");
  if (parts.cells === "spatialdata")
    copy("zarr/spatialdata_v3.zarr", "store.zarr");
  if (parts.cells === "h5ad") copy("h5ad/tiny.h5ad", "cells.h5ad");
  if (parts.molecules) copy("single-molecule/chunked", "molecules");
  for (const image of parts.images ?? []) writeImage(path.join(folder, image));

  return folder;
}

/** A folder that is itself the tiny OME-Zarr image. */
export function makeImageFolder(name = "dapi.zarr"): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "image-")),
    name,
  );

  writeImage(folder);

  return folder;
}

/**
 * Only the table out of the tiny SpatialData store: the table folder itself
 * ("table"), or the tables/ folder that holds it ("tables"). Neither has the
 * store's root metadata above it.
 */
export function makeTableFolder(which: "table" | "tables"): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "table-")),
    which,
  );

  fs.cpSync(
    path.join(
      TINY,
      "zarr",
      "spatialdata_v3.zarr",
      which === "table" ? "tables/table" : "tables",
    ),
    folder,
    { recursive: true },
  );

  return folder;
}

function writeImage(to: string, extent = 1000): void {
  fs.cpSync(
    path.join(TINY, "zarr", "spatialdata_v3.zarr", "images", "dapi"),
    to,
    { recursive: true },
  );
  // Written over 0-100; stretched to -extent..extent
  const file = path.join(to, "zarr.json");
  const group = JSON.parse(fs.readFileSync(file, "utf8"));

  for (const level of group.attributes.ome.multiscales[0].datasets) {
    const [scale, translation] = level.coordinateTransformations;

    scale.scale = scale.scale.map((s: number, axis: number) =>
      axis === 0 ? s : (s * extent) / 50,
    );
    translation.translation = translation.translation.map(
      (t: number, axis: number) =>
        axis === 0 ? t : (t * extent) / 50 - extent,
    );
  }
  fs.writeFileSync(file, JSON.stringify(group));
}
