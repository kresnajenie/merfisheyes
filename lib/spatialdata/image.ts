import type { SpatialDataStore } from "./store";

// Zarr v3 needs a newer zarrita than the one the h5ad-zarr path is pinned to.
import * as zarr from "zarrita-v3";

import { elementNames, publicUrlStore } from "./store";

import {
  fetchDatasetKeyList,
  PresignedFetchStore,
} from "@/lib/storage/PresignedFetchStore";

type PixelArray = zarr.Array<"uint8" | "uint16", zarr.Readable>;

export interface ImageLevel {
  array: PixelArray;
  width: number;
  height: number;
  chunkWidth: number;
  chunkHeight: number;
  /** Scene units per pixel of this level. */
  pixelWidth: number;
  pixelHeight: number;
  /** Scene position of the centre of this level's pixel (0, 0). */
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
}

type Transform = { type: string; scale?: number[]; translation?: number[] };

/** Per-axis scale and offset of a chain of scale / translation / identity. */
function composeTransforms(transforms: Transform[] | undefined, axes: number) {
  const scale = new Array(axes).fill(1);
  const offset = new Array(axes).fill(0);

  for (const t of transforms ?? []) {
    if (t.type === "scale" && t.scale) {
      for (let i = 0; i < axes; i++) {
        scale[i] *= t.scale[i];
        offset[i] *= t.scale[i];
      }
    } else if (t.type === "translation" && t.translation) {
      for (let i = 0; i < axes; i++) offset[i] += t.translation[i];
    } else if (t.type !== "identity") {
      console.warn(`[spatialdata] Ignoring unsupported transform "${t.type}"`);
    }
  }

  return { scale, offset };
}

/**
 * Open the store's multiscale image (OME-Zarr `images/<name>`, axes c/y/x).
 *
 * Levels are placed in the scene's units — those of the cells, i.e. the
 * intrinsic coordinates of the store's shapes (or points) element — by
 * undoing that element's transform to the shared coordinate system.
 */
export async function openSpatialDataImage({
  store,
  keys,
}: SpatialDataStore): Promise<SpatialDataImage | null> {
  const name = elementNames(keys, "images")[0];

  if (!name) return null;

  const root = zarr.root(store as zarr.Readable);
  const group = await zarr.open.v3(root.resolve(`images/${name}`), {
    kind: "group",
  });

  // Scene units → shared coordinate system, from the cells' own element
  const [shapes] = elementNames(keys, "shapes");
  const [points] = elementNames(keys, "points");
  const cellsPath = shapes ? `shapes/${shapes}` : points && `points/${points}`;
  let toShared = { scale: [1, 1], offset: [0, 0] };

  if (cellsPath) {
    const cells = await zarr.open.v3(root.resolve(cellsPath), {
      kind: "group",
    });

    // Axes are x, y(, z)
    toShared = composeTransforms(
      (cells.attrs.coordinateTransformations as Transform[]) ?? [],
      2,
    );
  }

  return readMultiscale(name, group, toShared);
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

  return readMultiscale(name, group, { scale: [1, 1], offset: [0, 0] });
}

/**
 * Open an image dataset that lives in the app's own bucket (uploaded through
 * the app). The image group may sit at the dataset root or in a folder under
 * it, such as `dapi.zarr/` — it is found by its `zarr.json`.
 */
export async function openOmeZarrImageDataset(
  datasetId: string,
): Promise<SpatialDataImage> {
  const marker = (await fetchDatasetKeyList(datasetId))
    .filter((k) => k === "zarr.json" || k.endsWith("/zarr.json"))
    .sort((a, b) => a.length - b.length)[0];

  if (marker === undefined) {
    throw new Error("No zarr v3 image found for this dataset");
  }
  const prefix = marker.slice(0, -"zarr.json".length);
  const store = new PresignedFetchStore(datasetId, { keyPrefix: prefix });
  const group = await zarr.open.v3(zarr.root(store), { kind: "group" });
  const name = prefix.replace(/\/$/, "").replace(/\.zarr$/, "") || datasetId;

  return readMultiscale(name, group, { scale: [1, 1], offset: [0, 0] });
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

/** Levels of an OME-Zarr multiscale group (axes c/y/x), in scene units. */
async function readMultiscale(
  name: string,
  group: zarr.Group<zarr.Readable>,
  toShared: { scale: number[]; offset: number[] },
): Promise<SpatialDataImage> {
  const ome = group.attrs.ome as any;
  const multiscale = ome.multiscales[0];
  const axisNames = multiscale.axes.map((a: { name: string }) => a.name);

  if (axisNames.join("") !== "cyx") {
    throw new Error(`Unsupported image axes: ${axisNames.join(", ")}`);
  }

  const imageToShared = multiscale.coordinateTransformations as Transform[];
  const levels: ImageLevel[] = [];

  for (const dataset of multiscale.datasets) {
    const array = await zarr.open.v3(group.resolve(dataset.path), {
      kind: "array",
    });

    if (!array.is("uint8") && !array.is("uint16")) {
      throw new Error(`Unsupported image data type: ${array.dtype}`);
    }
    // Level pixel index → shared coordinates (axes c, y, x)
    const t = composeTransforms(
      [...dataset.coordinateTransformations, ...(imageToShared ?? [])],
      3,
    );

    levels.push({
      array,
      height: array.shape[1],
      width: array.shape[2],
      chunkHeight: array.chunks[1],
      chunkWidth: array.chunks[2],
      pixelWidth: t.scale[2] / toShared.scale[0],
      pixelHeight: t.scale[1] / toShared.scale[1],
      originX: (t.offset[2] - toShared.offset[0]) / toShared.scale[0],
      originY: (t.offset[1] - toShared.offset[1]) / toShared.scale[1],
    });
  }

  return {
    name,
    channelLabels: (ome.omero?.channels ?? []).map(
      (c: { label?: string }, i: number) => c.label ?? `Channel ${i}`,
    ),
    maxValue: levels[0].array.is("uint8") ? 255 : 65535,
    levels,
  };
}
