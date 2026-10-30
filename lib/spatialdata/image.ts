import type { SpatialDataStore } from "./store";

// Zarr v3 needs a newer zarrita than the one the h5ad-zarr path is pinned to.
import * as zarr from "zarrita-v3";

import { elementNames, findZarrRoot, publicUrlStore } from "./store";
import { fileStore } from "./localStore";
import { axisNames, omeAttrs } from "./ome";
import {
  type Map2,
  IDENTITY,
  applyMap,
  chain,
  chooseImage,
  elementTransforms,
  readCoordinateSystems,
  outputName,
  toSystem,
} from "./coordinates";

type PixelArray = zarr.Array<"uint8" | "uint16", zarr.Readable>;

export interface ImageLevel {
  array: PixelArray;
  width: number;
  height: number;
  chunkWidth: number;
  chunkHeight: number;
  /** Image units per pixel of this level (the image's own, axis-aligned frame). */
  pixelWidth: number;
  pixelHeight: number;
  /** Position, in image units, of the centre of this level's pixel (0, 0). */
  originX: number;
  originY: number;
}

export interface SpatialDataImage {
  name: string;
  channelLabels: string[];
  /** Largest value a pixel can hold (255 or 65535). */
  maxValue: number;
  /** Finest first. */
  levels: ImageLevel[];
  /**
   * Image units → the scene (the cells' coordinates). The levels stay in
   * the image's own frame, so a placement that rotates (an affine) is one
   * matrix on the layer, not a different rectangle per tile.
   */
  placement: Map2;
}

/** The image's extent in the scene: the finest level's rectangle, placed. */
export function imageBounds(image: SpatialDataImage): {
  minX: number;
  minY: number;
  maxX: number;
  maxY: number;
} {
  const l = image.levels[0];
  const x0 = l.originX - l.pixelWidth / 2;
  const y0 = l.originY - l.pixelHeight / 2;
  const x1 = x0 + l.width * l.pixelWidth;
  const y1 = y0 + l.height * l.pixelHeight;
  const corners = [
    [x0, y0],
    [x1, y0],
    [x0, y1],
    [x1, y1],
  ].map(([x, y]) => applyMap(image.placement, x, y));

  return {
    minX: Math.min(...corners.map((c) => c[0])),
    minY: Math.min(...corners.map((c) => c[1])),
    maxX: Math.max(...corners.map((c) => c[0])),
    maxY: Math.max(...corners.map((c) => c[1])),
  };
}

/**
 * Open the store's multiscale image (OME-Zarr `images/<name>`, axes c/y/x),
 * placed in the scene's units — the cells' own coordinates — through the
 * store's coordinate systems (`coordinates.ts`): the finest image that
 * reaches the target is drawn; one that cannot is not opened, and the error
 * says why.
 */
export async function openSpatialDataImage(
  sd: SpatialDataStore,
): Promise<SpatialDataImage | null> {
  const names = elementNames(sd.keys, "images");

  if (!names.length) return null;
  const cs = await readCoordinateSystems(sd);
  const chosen = chooseImage(cs, names);

  if (chosen.name === null) {
    throw new Error(`No image can be placed with the cells: ${chosen.reason}`);
  }
  const placed = cs.place(`images/${chosen.name}`);

  if (!placed.ok) throw new Error(placed.reason);
  const root = zarr.root(sd.store as zarr.Readable);
  const group = await zarr.open.v3(root.resolve(`images/${chosen.name}`), {
    kind: "group",
  });

  return readMultiscale(chosen.name, group, placed.map);
}

/**
 * A standalone image's own placement, after each level's map: a plain
 * OME-NGFF chain; or, when the entries name output systems (as
 * `spatialdata` writes them), the one into `global`, else the first.
 */
function standalonePlacement(group: zarr.Group<zarr.Readable>): Map2 {
  const el = elementTransforms("images", group.attrs as Record<string, any>);
  const systems = el.transforms
    .map(outputName)
    .filter((n): n is string => n !== null);
  const p = toSystem(
    el,
    systems.includes("global") ? "global" : (systems[0] ?? null),
  );

  if (!p.ok) throw new Error(`Image transform: ${p.reason}`);

  return p.map;
}

/**
 * Open a standalone OME-Zarr multiscale image at a public URL (the image
 * group itself, not a SpatialData store). Its transforms must already place
 * it in the scene's units — there is no cells element to undo.
 */
export async function openOmeZarrImage(url: string): Promise<SpatialDataImage> {
  const base = url.replace(/\/+$/, "");
  const root = zarr.root(publicUrlStore(base) as zarr.Readable);
  const group = await zarr.open.v3(root, { kind: "group" });
  const name = base.slice(base.lastIndexOf("/") + 1).replace(/\.zarr$/, "");

  return readMultiscale(name, group, standalonePlacement(group));
}

/**
 * Open an image dataset that lives in the app's own bucket (uploaded through
 * the app). The image group may sit at the dataset root or in a folder under
 * it, such as `dapi.zarr/` — it is found by its `zarr.json`.
 */
export async function openOmeZarrImageDataset(
  datasetId: string,
): Promise<SpatialDataImage> {
  const { prefix, store, rootBytes } = await findZarrRoot(datasetId);

  store.seed("/zarr.json", rootBytes);
  const group = await zarr.open.v3(zarr.root(store), { kind: "group" });
  const name = prefix.replace(/\/$/, "").replace(/\.zarr$/, "") || datasetId;

  return readMultiscale(name, group, standalonePlacement(group));
}

/**
 * Open an OME-Zarr image from files picked or dropped in the browser, keyed
 * by their path inside the image folder. Tiles are read from disk as the
 * view asks for them.
 */
export async function openOmeZarrImageFromFiles(
  files: Map<string, File>,
  name: string,
): Promise<SpatialDataImage> {
  const root = zarr.root(fileStore(files) as zarr.Readable);
  const group = await zarr.open.v3(root, { kind: "group" });

  return readMultiscale(
    name.replace(/\.zarr$/, ""),
    group,
    standalonePlacement(group),
  );
}

/**
 * Where an image comes from: a public URL, or (`source: "app"`) the id of an
 * image dataset in the app's own bucket. This is how a `mapping.json` names
 * the image drawn under a dataset's cells.
 */
export interface ImageRef {
  image: string;
  source?: "s3" | "app";
}

export function openImageRef(ref: ImageRef): Promise<SpatialDataImage> {
  return ref.source === "app"
    ? openOmeZarrImageDataset(ref.image)
    : openOmeZarrImage(ref.image);
}

/**
 * Levels of an OME-Zarr multiscale group (axes c/y/x) in the image's own
 * units, and the placement that takes those units to the scene.
 */
async function readMultiscale(
  name: string,
  group: zarr.Group<zarr.Readable>,
  placement: Map2 = IDENTITY,
): Promise<SpatialDataImage> {
  const ome = omeAttrs(group.attrs);
  const multiscale = ome.multiscales[0];
  const axes = axisNames(multiscale);

  if (axes !== "cyx") {
    throw new Error(`Unsupported image axes: ${[...axes].join(", ")}`);
  }

  const levels: ImageLevel[] = [];

  for (const dataset of multiscale.datasets) {
    const array = await zarr.open.v3(group.resolve(dataset.path), {
      kind: "array",
    });

    if (!array.is("uint8") && !array.is("uint16")) {
      throw new Error(`Unsupported image data type: ${array.dtype}`);
    }
    // Level pixel index → image units (a level's own transforms are a
    // scale and a translation: pixels of a coarser level are larger)
    const t = chain(dataset.coordinateTransformations ?? [], [...axes]);

    levels.push({
      array,
      height: array.shape[1],
      width: array.shape[2],
      chunkHeight: array.chunks[1],
      chunkWidth: array.chunks[2],
      pixelWidth: t.scale[0],
      pixelHeight: t.scale[1],
      originX: t.offset[0],
      originY: t.offset[1],
    });
  }

  return {
    name,
    // A label may be written as a number (spatialdata's default: 0, 1, 2)
    channelLabels: (ome.omero?.channels ?? []).map(
      (c: { label?: string | number }, i: number) =>
        c.label != null ? String(c.label) : `Channel ${i}`,
    ),
    maxValue: levels[0].array.is("uint8") ? 255 : 65535,
    levels,
    placement,
  };
}
