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

/**
 * Which elements a SpatialData store written by `writeStore` holds. `table`
 * picks how its X is kept: "csr" with a CSC copy in layers/ (what
 * prepare_folder.py, prepare_xenium.py and prepare_visium.py write) or
 * "csc_x", X itself CSC (what prepare_merscope.py writes).
 */
export interface StoreSpec {
  table?: "csr" | "csc_x";
  image?: boolean;
  /** The tiled points element (docs/points-zarr-spec.md), points_v1.zarr's. */
  points?: boolean;
}

const STORE_SOURCES = {
  csr: path.join(TINY, "zarr", "spatialdata_v3.zarr"),
  csc_x: path.join(TINY, "zarr", "spatialdata_v3_csc_x.zarr"),
  image: path.join(TINY, "zarr", "spatialdata_v3.zarr"),
  points: path.join(TINY, "zarr", "points_v1.zarr"),
};

/**
 * Write a SpatialData store holding only the elements of `spec`, composed
 * from the committed tiny stores, with consolidated metadata naming exactly
 * what is in it (as e2e/helpers/spatialdata.ts serves them by URL).
 */
export function writeStore(to: string, spec: StoreSpec): void {
  const read = (store: string) =>
    JSON.parse(fs.readFileSync(path.join(store, "zarr.json"), "utf8"));
  const root = read(STORE_SOURCES.csr);
  const nodes: Record<string, unknown> = {};
  const take = (store: string, group: string) => {
    fs.cpSync(path.join(store, group), path.join(to, group), {
      recursive: true,
    });
    for (const [key, node] of Object.entries(
      read(store).consolidated_metadata.metadata,
    )) {
      if (key === group || key.startsWith(`${group}/`)) nodes[key] = node;
    }
  };

  fs.mkdirSync(to, { recursive: true });
  if (spec.table) take(STORE_SOURCES[spec.table], "tables");
  if (spec.image) take(STORE_SOURCES.image, "images");
  if (spec.points) take(STORE_SOURCES.points, "points");
  root.consolidated_metadata.metadata = nodes;
  fs.writeFileSync(path.join(to, "zarr.json"), JSON.stringify(root));
}

/** A dropped folder that is itself a store holding `spec`. */
export function makeStoreFolder(spec: StoreSpec, name = "sample.zarr"): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "store-")),
    name,
  );

  writeStore(folder, spec);

  return folder;
}

/**
 * A dropped folder holding each element as a folder of its own, side by side:
 * the table as a store of its own (`table.zarr`), the image as a standalone
 * OME-Zarr (`image.zarr`, over the table's 0-100 extent) and the points as a
 * store of their own (`points.zarr`).
 */
export function makeSeparateFolders(spec: StoreSpec, name = "sample"): string {
  const folder = path.join(
    fs.mkdtempSync(path.join(os.tmpdir(), "parts-")),
    name,
  );

  fs.mkdirSync(folder);
  if (spec.table) {
    writeStore(path.join(folder, "table.zarr"), { table: spec.table });
  }
  if (spec.image) {
    fs.cpSync(
      path.join(STORE_SOURCES.image, "images", "dapi"),
      path.join(folder, "image.zarr"),
      { recursive: true },
    );
  }
  if (spec.points) {
    writeStore(path.join(folder, "points.zarr"), { points: true });
  }

  return folder;
}

/**
 * Drop folders on a homepage drop zone the way a person drags them in from
 * the desktop: one drop event, one item per folder. Playwright cannot hand a
 * real directory to a DataTransfer, so the browser's entry API
 * (`webkitGetAsEntry`, `createReader().readEntries`, `file()`) is stood in for
 * with the folders' real files; everything after it — the walk in
 * file-upload.tsx and all routing — is the app's own.
 */
export async function dropFolders(
  page: import("@playwright/test").Page,
  dropZone: import("@playwright/test").Locator,
  folders: string[],
): Promise<void> {
  const trees = folders.map((folder) => ({
    name: path.basename(folder),
    files: (fs.readdirSync(folder, { recursive: true }) as string[])
      .filter((rel) => fs.statSync(path.join(folder, rel)).isFile())
      .map((rel) => ({
        path: rel.split(path.sep).join("/"),
        data: fs.readFileSync(path.join(folder, rel)).toString("base64"),
      })),
  }));

  await dropZone.evaluate((zone, trees) => {
    type Node = { name: string; files?: Map<string, Node>; file?: File };
    const dir = (name: string): Node => ({ name, files: new Map() });
    const entry = (node: Node): unknown =>
      node.file
        ? {
            name: node.name,
            isFile: true,
            isDirectory: false,
            file: (ok: (f: File) => void) => ok(node.file!),
          }
        : {
            name: node.name,
            isFile: false,
            isDirectory: true,
            createReader: () => {
              let done = false;

              return {
                readEntries: (ok: (e: unknown[]) => void) => {
                  ok(done ? [] : Array.from(node.files!.values()).map(entry));
                  done = true;
                },
              };
            },
          };
    const roots = trees.map((tree) => {
      const root = dir(tree.name);

      for (const { path: rel, data } of tree.files) {
        const parts = rel.split("/");
        let at = root;

        for (const part of parts.slice(0, -1)) {
          if (!at.files!.has(part)) at.files!.set(part, dir(part));
          at = at.files!.get(part)!;
        }
        const name = parts[parts.length - 1];
        const bytes = Uint8Array.from(atob(data), (c) => c.charCodeAt(0));

        at.files!.set(name, { name, file: new File([bytes], name) });
      }

      return root;
    });
    const transfer = new DataTransfer();

    for (const root of roots) transfer.items.add(new File([], root.name));

    // Each item answers with its folder's entry, as a real dropped folder does
    Object.defineProperty(DataTransferItem.prototype, "webkitGetAsEntry", {
      configurable: true,
      value(this: DataTransferItem) {
        const name = this.getAsFile()?.name;

        return entry(roots.find((root) => root.name === name)!);
      },
    });
    zone.dispatchEvent(
      new DragEvent("drop", {
        bubbles: true,
        cancelable: true,
        dataTransfer: transfer,
      }),
    );
  }, trees);
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
