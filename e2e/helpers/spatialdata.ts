import fs from "node:fs";
import path from "node:path";
import { gzipSync } from "node:zlib";

import { REPO_ROOT } from "./databank";
import {
  type FixtureServer,
  type FixtureServerOptions,
  type Resolved,
  startFixtureServer,
} from "./fixture-server";

const TINY = path.join(REPO_ROOT, "tests", "data", "tiny");
const ZARR = path.join(TINY, "zarr");
const MEDIUM = path.join(
  REPO_ROOT,
  "tests",
  "data",
  "generated",
  "medium",
  "points-store",
);

/** The two committed stores: the same cells and image, laid out differently. */
export type StoreLayout = "spatialdata_v3" | "spatialdata_v3_sharded";

/**
 * Which of cells, molecules, points, shapes and an image a store is served
 * with. Molecules are the tiny committed per-gene folder, or (`"large"`) a
 * generated one whose biggest genes are streamed — see `LARGE_MOLECULES`.
 * Points is the committed points element (`points_v1.zarr`) merged into the
 * store; shapes the committed polygons (`shapes_v1`).
 */
export interface StoreParts {
  cells?: boolean;
  molecules?: boolean | "large";
  points?: boolean;
  shapes?: boolean;
  image?: boolean;
}

/**
 * The generated molecule folder served at `/molecules-large/<tag>`: four
 * genes, two of them past the streaming threshold. `<tag>` is appended to
 * every gene name, so two panels on one page can show two such folders
 * without their genes' test hooks colliding. The default genes picked on
 * open are the first three: both large ones and one small one.
 */
export const LARGE_MOLECULES = {
  /** Molecules in each of the two streamed genes. */
  large: 150_000,
  /** Molecules in each of the two whole-loaded genes. */
  small: 200,
  genes: (tag: string) => [
    `Big_${tag}`,
    `Huge_${tag}`,
    `Small_${tag}`,
    `Tiny_${tag}`,
  ],
  /**
   * Gene files are trickled out so a download is observably in progress:
   * a large gene (1.8 MB) takes about 2 s.
   */
  throttle: { chunkBytes: 32 * 1024, delayMs: 40 },
};

/** Gene files and the manifest of one generated folder, gzipped. */
function largeMoleculeFolder(tag: string): Map<string, Buffer> {
  const genes = LARGE_MOLECULES.genes(tag);
  const files = new Map<string, Buffer>();
  const counts: Record<string, { assigned: number }> = {};
  // Deterministic: the same folder on every run.
  let seed = 12345;
  const random = () => {
    seed = (seed * 1664525 + 1013904223) >>> 0;

    return seed / 2 ** 32;
  };

  genes.forEach((gene, i) => {
    const n = i < 2 ? LARGE_MOLECULES.large : LARGE_MOLECULES.small;
    const coords = new Float32Array(n * 3);

    for (let m = 0; m < n; m++) {
      coords[m * 3] = random() * 2000 - 1000;
      coords[m * 3 + 1] = random() * 2000 - 1000;
      coords[m * 3 + 2] = random() * 10;
    }
    counts[gene] = { assigned: n };
    files.set(`genes/${gene}.bin.gz`, gzipSync(Buffer.from(coords.buffer)));
  });

  const total = Object.values(counts).reduce((a, c) => a + c.assigned, 0);

  files.set(
    "manifest.json.gz",
    gzipSync(
      JSON.stringify({
        version: "1.0",
        dataset_id: `large_${tag}`,
        name: `large_${tag}`,
        type: "xenium",
        statistics: {
          total_molecules: total,
          unique_genes: genes.length,
          spatial_dimensions: 3,
        },
        genes: { unique_gene_names: genes, molecule_counts: counts },
        has_unassigned: false,
        processing: {
          compression: "gzip",
          coordinate_format: "float32_flat_array",
          coordinate_range: "raw_rounded_2dp",
          scaling_factor: 1,
        },
      }),
    ),
  );

  return files;
}

/** Consolidated metadata of a store on disk, as `{ path: node }`. */
function consolidated(store: string): Record<string, unknown> {
  return JSON.parse(fs.readFileSync(path.join(store, "zarr.json"), "utf8"))
    .consolidated_metadata.metadata;
}

/** The shapes group as a consolidated node. */
function shapesNode(): unknown {
  return JSON.parse(
    fs.readFileSync(
      path.join(ZARR, "shapes_v1", "cell_boundaries", "zarr.json"),
      "utf8",
    ),
  );
}

/**
 * Serve the tiny SpatialData fixtures the way a public bucket would:
 *
 *   /store/<layout>/<parts>/store.zarr   a store holding only `<parts>`
 *   /store/<layout>/<parts>/mapping.json names the store, and its link to
 *                                        the molecules, if any
 *   /molecules                           the processed molecules it links to
 *   /molecules-large/<tag>               a generated molecule folder with
 *                                        genes big enough to be streamed
 *   /cells/<tag>                         the tiny merfisheyes cell folder,
 *                                        with a mapping.json linking it to
 *                                        /molecules-large/<tag>
 *   /image/<extent>.zarr                 the image on its own, `<extent>` wide
 *   /visium/store.zarr                   the Visium-shaped store spatialdata
 *                                        wrote (make-tiny-visium.py): three
 *                                        images in three systems, spots as
 *                                        circles, no `global`
 *   /visium-hd/store.zarr                the Visium HD-shaped one: square
 *                                        bins of side 16 tiling a grid
 *   /codex/store.zarr, /maldi/store.zarr, /cosmx/store.zarr
 *                                        stores shaped like SpatialData-DB's
 *                                        CODEX (dense X, a 6-channel image),
 *                                        MALDI-MSI (`formula_adduct` features,
 *                                        an affine image) and CosMx (no
 *                                        obsm/spatial) ones
 *                                        (make-tiny-modalities.py)
 *   /medium/store.zarr                   the generated medium store, when
 *                                        it is on disk (make-medium-store.py)
 *
 * `<parts>` is any of the letters c, m (M for the large molecules), p
 * (the points element), s (shapes), i (or `none`). The viewer only knows a
 * store by the consolidated metadata in its root `zarr.json`, so one store
 * on disk stands in for every combination by trimming that listing — and
 * the points and shapes fixtures, which live in their own folders, are
 * merged into it and routed to their files.
 *
 * `options.bytesPerSecond` caps the server's send rate (see fixture-server).
 */
export function startSpatialDataServer(
  options: FixtureServerOptions = {},
): Promise<FixtureServer> {
  let server: FixtureServer;
  const largeFolders = new Map<string, Map<string, Buffer>>();
  const largeMolecules = (tag: string, file: string): Resolved => {
    let folder = largeFolders.get(tag);

    if (!folder) {
      folder = largeMoleculeFolder(tag);
      largeFolders.set(tag, folder);
    }
    const bytes = folder.get(file);

    if (!bytes) return undefined;

    return file.startsWith("genes/")
      ? { bytes, throttle: LARGE_MOLECULES.throttle }
      : { bytes };
  };

  const resolve = (urlPath: string): Resolved => {
    const [, kind, ...rest] = urlPath.split("/");

    // A folder holding a store of cells alone, whose mapping.json lends it
    // the molecules of a store with cells too, linked by that store's folder
    // (as a store registered by its public URL is linked)
    if (kind === "linked") {
      if (rest.join("/") === "mapping.json") {
        return {
          body: JSON.stringify({
            cells: "store.zarr",
            linkColumn: "__all__",
            links: {
              __all__: storeFolderUrl(server, { cells: true, points: true }),
            },
            source: "s3",
            linkFormat: "store",
          }),
        };
      }

      return resolve(`/store/spatialdata_v3/c/${rest.join("/")}`);
    }

    if (kind === "molecules") {
      return { file: path.join(TINY, "single-molecule", "chunked", ...rest) };
    }
    if (kind === "molecules-large") {
      const [tag, ...inside] = rest;

      return largeMolecules(tag, inside.join("/"));
    }
    if (kind === "cells") {
      const [tag, ...inside] = rest;

      if (inside.join("/") === "mapping.json") {
        return {
          body: JSON.stringify({
            linkColumn: "__all__",
            links: { __all__: `${server.url}/molecules-large/${tag}` },
          }),
        };
      }

      return { file: path.join(TINY, "chunked", "raw", ...inside) };
    }
    if (kind === "medium") {
      return { file: path.join(MEDIUM, ...rest) };
    }
    if (kind === "parity") {
      return { file: path.join(TINY, "parity", ...rest) };
    }
    if (kind in LIBRARY_STORES) {
      const [, ...inside] = rest;

      return {
        file: path.join(
          ZARR,
          LIBRARY_STORES[kind as keyof typeof LIBRARY_STORES],
          ...inside,
        ),
      };
    }
    if (kind === "image") {
      const [name, ...inside] = rest;
      const file = path.join(
        ZARR,
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
    // The same image, its metadata written the way the CoQ1s stitcher wrote
    // it: NGFF 0.4 attributes at the top of a Zarr v3 group, no `ome` wrapper
    // (docs/specs/image-profile.md, docs/plans/OME-ZARR-COMPAT.md)
    if (kind === "image04") {
      const inner = resolve(`/image/${rest.join("/")}`);

      if (!inner || !("body" in inner)) return inner;
      const group = JSON.parse(inner.body);
      const { multiscales, omero } = group.attributes.ome;

      group.attributes = { multiscales, omero };
      group.attributes.multiscales[0].version = "0.4";

      return { body: JSON.stringify(group) };
    }
    if (kind !== "store") return undefined;

    const [layout, parts, file, ...inside] = rest;

    if (file === "mapping.json") {
      const molecules = parts.includes("M")
        ? `${server.url}/molecules-large/store`
        : parts.includes("m")
          ? `${server.url}/molecules`
          : null;

      return {
        body: JSON.stringify({
          cells: "store.zarr",
          ...(molecules && {
            linkColumn: "__all__",
            links: { __all__: molecules },
          }),
        }),
      };
    }

    const key = inside.join("/");

    // The points element and the shapes live in their own fixture folders
    if (key.startsWith("points/")) {
      return { file: path.join(ZARR, "points_v1.zarr", ...inside) };
    }
    if (key.startsWith("shapes/cell_boundaries/")) {
      // shapes_v1 holds the element folder itself (cell_boundaries/…)
      return { file: path.join(ZARR, "shapes_v1", ...inside.slice(1)) };
    }

    const onDisk = path.join(ZARR, `${layout}.zarr`, ...inside);

    if (key !== "zarr.json") return { file: onDisk };

    const root = JSON.parse(fs.readFileSync(onDisk, "utf8"));
    const keep = [
      ...(parts.includes("c") ? ["tables"] : []),
      ...(parts.includes("i") ? ["images"] : []),
    ];
    const nodes: Record<string, unknown> = Object.fromEntries(
      Object.entries(root.consolidated_metadata.metadata).filter(([k]) =>
        keep.includes(k.split("/")[0]),
      ),
    );

    if (parts.includes("p")) {
      for (const [k, node] of Object.entries(
        consolidated(path.join(ZARR, "points_v1.zarr")),
      )) {
        if (k === "points" || k.startsWith("points/")) nodes[k] = node;
      }
    }
    if (parts.includes("s")) {
      nodes["shapes"] = { zarr_format: 3, node_type: "group" };
      nodes["shapes/cell_boundaries"] = shapesNode();
    }
    root.consolidated_metadata.metadata = nodes;

    return { body: JSON.stringify(root) };
  };

  return startFixtureServer(resolve, options).then(
    (started) => (server = started),
  );
}

/** The letters of a parts spec, in the order the paths use. */
function letters(parts: StoreParts): string {
  return (
    (parts.cells ? "c" : "") +
    (parts.molecules === "large" ? "M" : parts.molecules ? "m" : "") +
    (parts.points ? "p" : "") +
    (parts.shapes ? "s" : "") +
    (parts.image ? "i" : "")
  );
}

/** URL of the tiny store served with only the given parts. */
export function storeUrl(
  server: FixtureServer,
  parts: StoreParts,
  layout: StoreLayout = "spatialdata_v3",
): string {
  return `${server.url}/store/${layout}/${letters(parts) || "none"}/store.zarr`;
}

/** URL of the folder that holds that store, and names it in mapping.json. */
export function storeFolderUrl(
  server: FixtureServer,
  parts: StoreParts,
  layout: StoreLayout = "spatialdata_v3",
): string {
  return storeUrl(server, parts, layout).replace(/\/store\.zarr$/, "/");
}

/** URL of the generated molecule folder tagged `tag` (see LARGE_MOLECULES). */
export function largeMoleculesUrl(server: FixtureServer, tag: string): string {
  return `${server.url}/molecules-large/${tag}`;
}

/** URL of the tiny cell folder whose mapping.json links `/molecules-large/<tag>`. */
export function cellsUrl(server: FixtureServer, tag: string): string {
  return `${server.url}/cells/${tag}`;
}

/**
 * The stores the `spatialdata` library wrote, served whole at
 * `/<kind>/store.zarr` (make-tiny-visium.py, make-tiny-modalities.py).
 */
export const LIBRARY_STORES = {
  visium: "visium_v3.zarr",
  "visium-hd": "visium_hd_v3.zarr",
  codex: "codex_v3.zarr",
  maldi: "maldi_v3.zarr",
  cosmx: "cosmx_v3.zarr",
} as const;

/** URL of the Visium-shaped store (no `global` system). */
export function visiumUrl(server: FixtureServer): string {
  return `${server.url}/visium/store.zarr`;
}

/**
 * URL of a folder holding a store of cells alone, lent the molecules of a
 * store that has cells too, linked by that store's folder.
 */
export function linkedStoreFolderUrl(server: FixtureServer): string {
  return `${server.url}/linked/`;
}

/**
 * URL of a store `scripts/prepare.py` wrote from a tiny source
 * (tests/data/tiny/parity/<name>/<name>.zarr, made by
 * make-tiny-parity-stores.py): what an upload saves. `merscope` and
 * `tiny_h5ad` carry a UMAP.
 */
export function parityStoreUrl(server: FixtureServer, name: string): string {
  return `${server.url}/parity/${name}/${name}.zarr`;
}

/** URL of the Visium HD-shaped store (square bins). */
export function visiumHdUrl(server: FixtureServer): string {
  return `${server.url}/visium-hd/store.zarr`;
}

/** URL of a store shaped like a SpatialData-DB CODEX, MALDI-MSI or CosMx one. */
export function libraryStoreUrl(
  server: FixtureServer,
  kind: "codex" | "maldi" | "cosmx",
): string {
  return `${server.url}/${kind}/store.zarr`;
}

/**
 * URL of the tiny image on its own, covering -extent..extent on both axes
 * (1000 puts it under the tiny h5ad's cells).
 */
export function imageUrl(server: FixtureServer, extent: number): string {
  return `${server.url}/image/${extent}.zarr`;
}

/**
 * The same image with NGFF 0.4 attributes at the top of a Zarr v3 group, as
 * the CoQ1s stitcher wrote (docs/plans/OME-ZARR-COMPAT.md). Same pixels.
 */
export function image04Url(server: FixtureServer, extent: number): string {
  return `${server.url}/image04/${extent}.zarr`;
}
