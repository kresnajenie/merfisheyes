import fs from "node:fs";
import path from "node:path";

import { REPO_ROOT } from "./databank";
import { type FixtureServer, startFixtureServer } from "./fixture-server";

const TINY = path.join(REPO_ROOT, "tests", "data", "tiny");

/** The two committed stores: the same cells and image, laid out differently. */
export type StoreLayout = "spatialdata_v3" | "spatialdata_v3_sharded";

/** Which of cells, molecules and an image a store is served with. */
export interface StoreParts {
  cells?: boolean;
  molecules?: boolean;
  image?: boolean;
}

/**
 * Serve the tiny SpatialData fixtures the way a public bucket would:
 *
 *   /store/<layout>/<parts>/store.zarr   a store holding only `<parts>`
 *   /store/<layout>/<parts>/mapping.json its link to the molecules, if any
 *   /molecules                           the processed molecules it links to
 *   /image/<extent>.zarr                 the image on its own, `<extent>` wide
 *
 * `<parts>` is any of the letters c, m, i (or `none`). The viewer only knows
 * a store by the consolidated metadata in its root `zarr.json`, so one store
 * on disk stands in for every combination by trimming that listing.
 */
export function startSpatialDataServer(): Promise<FixtureServer> {
  let server: FixtureServer;

  return startFixtureServer((urlPath) => {
    const [, kind, ...rest] = urlPath.split("/");

    if (kind === "molecules") {
      return { file: path.join(TINY, "single-molecule", "chunked", ...rest) };
    }
    if (kind === "image") {
      const [name, ...inside] = rest;
      const file = path.join(
        TINY,
        "zarr",
        "spatialdata_v3.zarr",
        "images",
        "dapi",
        ...inside,
      );

      if (inside.join("/") !== "zarr.json") return { file };

      // The same pixels, stretched from the 0-100 they were written over to
      // -extent..extent, to sit under cells that span that much.
      const extent = Number(name.replace(/\.zarr$/, ""));
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

      return { body: JSON.stringify(group) };
    }
    if (kind !== "store") return undefined;

    const [layout, parts, file, ...inside] = rest;

    if (file === "mapping.json") {
      return parts.includes("m")
        ? {
            body: JSON.stringify({
              linkColumn: "__all__",
              links: { __all__: `${server.url}/molecules` },
            }),
          }
        : undefined;
    }

    const onDisk = path.join(TINY, "zarr", `${layout}.zarr`, ...inside);

    if (inside.join("/") !== "zarr.json") return { file: onDisk };

    const root = JSON.parse(fs.readFileSync(onDisk, "utf8"));
    const keep = [
      ...(parts.includes("c") ? ["tables"] : []),
      ...(parts.includes("i") ? ["images"] : []),
    ];

    root.consolidated_metadata.metadata = Object.fromEntries(
      Object.entries(root.consolidated_metadata.metadata).filter(([key]) =>
        keep.includes(key.split("/")[0]),
      ),
    );

    return { body: JSON.stringify(root) };
  }).then((started) => (server = started));
}

/** URL of the tiny store served with only the given parts. */
export function storeUrl(
  server: FixtureServer,
  parts: StoreParts,
  layout: StoreLayout = "spatialdata_v3",
): string {
  const letters =
    (parts.cells ? "c" : "") +
    (parts.molecules ? "m" : "") +
    (parts.image ? "i" : "");

  return `${server.url}/store/${layout}/${letters || "none"}/store.zarr`;
}

/**
 * URL of the tiny image on its own, covering -extent..extent on both axes
 * (1000 puts it under the tiny h5ad's cells).
 */
export function imageUrl(server: FixtureServer, extent: number): string {
  return `${server.url}/image/${extent}.zarr`;
}
