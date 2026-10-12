import { isOmeZarrImage } from "@/lib/upload/uploadImageToS3";

/**
 * Find the already-processed parts in a folder dropped for a browser preview:
 * cells, molecules and an image, any of which may be missing. Nothing here
 * needs converting — each part is a format the viewer reads directly, a piece
 * at a time — which is what lets a folder of any size open at once.
 *
 * Parts are recognised by what is in them, not by their names:
 *
 *   cells      a SpatialData zarr store (tables/ or images/ inside it), the
 *              table of one on its own (an AnnData zarr v3: obs/, var/), or
 *              a merfisheyes cell folder (manifest.json, expr/, coords/)
 *   molecules  a merfisheyes molecule folder (manifest.json[.gz] + genes/)
 *   image      an OME-Zarr multiscale image (zarr.json naming axes c, y, x)
 *
 * Each is a folder directly inside the dropped one, or the dropped folder
 * itself. With several candidates for one part the first by name is used and
 * the others are reported as skipped.
 */

export type BundleCells =
  | { kind: "chunked"; root: string }
  | { kind: "spatialdata"; root: string }
  /** An AnnData table (Zarr v3) with no store around it. */
  | { kind: "table"; root: string };

export interface LocalBundle {
  /** The dropped folder's own name. */
  name: string;
  /** Every usable file, keyed by its path from the dropped folder. */
  files: Map<string, File>;
  cells?: BundleCells;
  /** Folder of the molecules ("" = the dropped folder itself). */
  molecules?: string;
  /** Folder of a standalone image. A SpatialData store's image is inside it. */
  image?: string;
  /** Candidates that lost to an earlier one, e.g. "image b.zarr". */
  skipped: string[];
}

/** macOS and archive-tool droppings; never part of a dataset. */
const isJunk = (key: string) =>
  key.split("/").some((s) => s.startsWith(".") || s === "__MACOSX");

const hasChunkedCells = (keys: Set<string>, p: string) =>
  keys.has(`${p}manifest.json`) &&
  keys.has(`${p}expr/index.json`) &&
  keys.has(`${p}coords/spatial.bin.gz`);

const hasMolecules = (keys: Set<string>, p: string) =>
  (keys.has(`${p}manifest.json.gz`) || keys.has(`${p}manifest.json`)) &&
  Array.from(keys).some(
    (k) => k.startsWith(`${p}genes/`) && k.endsWith(".bin.gz"),
  );

/** A SpatialData store keeps its elements in groups named by their type. */
const hasSpatialDataElements = (keys: Set<string>, p: string) =>
  Array.from(keys).some(
    (k) =>
      /^(tables|images|points|shapes|labels)\/.+\/zarr\.json$/.test(
        k.slice(p.length),
      ) && k.startsWith(p),
  );

export async function classifyBundle(files: File[]): Promise<LocalBundle> {
  const first = files[0]?.webkitRelativePath || files[0]?.name || "";
  const name = first.includes("/") ? first.split("/")[0] : "dataset";
  const all = new Map<string, File>();

  for (const file of files) {
    const path = file.webkitRelativePath || file.name;
    const key = path.includes("/") ? path.slice(path.indexOf("/") + 1) : path;

    if (!isJunk(key)) all.set(key, file);
  }

  const keys = new Set(all.keys());
  const bundle: LocalBundle = { name, files: all, skipped: [] };
  const isImage = async (prefix: string) =>
    isOmeZarrImage(
      await all
        .get(`${prefix}zarr.json`)!
        .text()
        .then(JSON.parse)
        .catch(() => null),
    );

  /** What the folder at `prefix` ("" or "name/") holds, if anything we read. */
  const partAt = async (
    prefix: string,
  ): Promise<
    | "cells-chunked"
    | "cells-spatialdata"
    | "cells-table"
    | "molecules"
    | "image"
    | null
  > => {
    // A store is known by its element folders; its own root metadata may be
    // missing when someone copied just tables/ or images/ into a folder.
    if (hasSpatialDataElements(keys, prefix)) return "cells-spatialdata";
    if (keys.has(`${prefix}zarr.json`)) {
      // The table of a store dropped on its own: AnnData written as Zarr v3.
      // (A Zarr v2 AnnData has no zarr.json and takes the older path.)
      if (
        keys.has(`${prefix}obs/zarr.json`) &&
        keys.has(`${prefix}var/zarr.json`)
      ) {
        return "cells-table";
      }

      return (await isImage(prefix)) ? "image" : null;
    }
    if (hasChunkedCells(keys, prefix)) return "cells-chunked";
    if (hasMolecules(keys, prefix)) return "molecules";

    return null;
  };

  // The dropped folder may itself be one part...
  const whole = await partAt("");

  if (whole === "image") return { ...bundle, image: "" };
  if (whole === "molecules") return { ...bundle, molecules: "" };
  if (whole === "cells-chunked") {
    return { ...bundle, cells: { kind: "chunked", root: "" } };
  }
  if (whole === "cells-spatialdata") {
    return { ...bundle, cells: { kind: "spatialdata", root: "" } };
  }
  if (whole === "cells-table") {
    return { ...bundle, cells: { kind: "table", root: "" } };
  }

  // ...or hold them as folders directly inside it.
  const folders = Array.from(
    new Set(
      Array.from(keys)
        .filter((k) => k.includes("/"))
        .map((k) => k.slice(0, k.indexOf("/"))),
    ),
  ).sort();

  for (const folder of folders) {
    const part = await partAt(`${folder}/`);

    if (part?.startsWith("cells-")) {
      if (bundle.cells) bundle.skipped.push(`cells ${folder}`);
      else {
        bundle.cells = {
          kind:
            part === "cells-chunked"
              ? "chunked"
              : part === "cells-table"
                ? "table"
                : "spatialdata",
          root: folder,
        };
      }
    } else if (part === "molecules") {
      if (bundle.molecules !== undefined) {
        bundle.skipped.push(`molecules ${folder}`);
      } else bundle.molecules = folder;
    } else if (part === "image") {
      if (bundle.image !== undefined) bundle.skipped.push(`image ${folder}`);
      else bundle.image = folder;
    }
  }

  return bundle;
}

/** Whether a classified folder holds anything the preview can open. */
export function hasParts(bundle: LocalBundle): boolean {
  return (
    bundle.cells !== undefined ||
    bundle.molecules !== undefined ||
    bundle.image !== undefined
  );
}

/**
 * The files of one part, as a folder of their own: each gets the relative
 * path it would have had if that part had been dropped by itself, which is
 * what the loaders for a single dropped folder expect.
 */
export function partFiles(bundle: LocalBundle, root: string): File[] {
  const prefix = root ? `${root}/` : "";
  const out: File[] = [];

  for (const [key, file] of bundle.files) {
    if (!key.startsWith(prefix)) continue;
    Object.defineProperty(file, "webkitRelativePath", {
      value: `${root || bundle.name}/${key.slice(prefix.length)}`,
      configurable: true,
    });
    out.push(file);
  }

  return out;
}

/** Group metadata for a folder that has none of its own. */
const groupFile = (attributes: Record<string, unknown> = {}) =>
  new File(
    [JSON.stringify({ zarr_format: 3, node_type: "group", attributes })],
    "zarr.json",
    { type: "application/json" },
  );

/**
 * The cells' zarr files laid out as a SpatialData store, which is what the
 * reader and the saved viewer expect: a store as it is, or a lone table put
 * where a store keeps its table.
 *
 * A lone table, or tables/ copied out of a store, has no group metadata
 * above it, and a store is known by its root's. Those two small files are
 * supplied; nothing else is added or changed.
 */
export function storeFileMap(
  bundle: LocalBundle,
  cells: Exclude<BundleCells, { kind: "chunked" }>,
): Map<string, File> {
  const part = partFileMap(bundle, cells.root);
  const files =
    cells.kind === "spatialdata"
      ? part
      : new Map(
          Array.from(part, ([key, file]) => [`tables/table/${key}`, file]),
        );

  if (!files.has("zarr.json")) {
    files.set(
      "zarr.json",
      groupFile({ spatialdata_attrs: { version: "0.2" } }),
    );
  }
  if (
    !files.has("tables/zarr.json") &&
    Array.from(files.keys()).some((k) => k.startsWith("tables/"))
  ) {
    files.set("tables/zarr.json", groupFile());
  }

  return files;
}

/** The files of one part keyed by their path inside it (zarr keys). */
export function partFileMap(
  bundle: LocalBundle,
  root: string,
): Map<string, File> {
  const prefix = root ? `${root}/` : "";
  const out = new Map<string, File>();

  for (const [key, file] of bundle.files) {
    if (key.startsWith(prefix)) out.set(key.slice(prefix.length), file);
  }

  return out;
}
