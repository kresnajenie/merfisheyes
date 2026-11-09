import { afterEach, describe, expect, it, vi } from "vitest";

import { imageUrlTargets } from "@/lib/overlays/imageSource";

/**
 * What a pasted link holds as an image (lib/overlays/imageSource.ts), the one
 * rule the overlay dialog and its route share. A store lends only its image,
 * whatever cells, molecules and shapes it holds beside it.
 */

const FOLDER = "https://bucket.example/spatialdata/sample";
const STORE = `${FOLDER}/sample.zarr`;
const IMAGE = "https://bucket.example/images/dapi.zarr";

/** An image group's metadata, in either layout (0.5 nests it under `ome`). */
const imageAttrs = (flat = false) => {
  const ome = { version: "0.5", multiscales: [{ axes: [], datasets: [] }] };

  return flat ? { multiscales: ome.multiscales } : { ome };
};

const storeRoot = (nodes: Record<string, unknown>) => ({
  attributes: { spatialdata_attrs: { version: "0.2" } },
  consolidated_metadata: { metadata: nodes },
});

/** A public bucket holding these files; anything else is a 404. */
function serve(files: Record<string, unknown>) {
  vi.stubGlobal(
    "fetch",
    vi.fn(async (url: string) =>
      url in files
        ? new Response(JSON.stringify(files[url]), { status: 200 })
        : new Response("", { status: 404 }),
    ),
  );
}

afterEach(() => vi.unstubAllGlobals());

describe("the image a public URL offers", () => {
  it("is the link itself, when the link is an image", async () => {
    serve({ [`${IMAGE}/zarr.json`]: { attributes: imageAttrs() } });

    expect(await imageUrlTargets(IMAGE)).toEqual([{ url: IMAGE }]);
  });

  it("is the image element of a store that holds cells and molecules too", async () => {
    serve({
      [`${STORE}/zarr.json`]: storeRoot({
        "tables/table": { attributes: {} },
        "points/transcripts": { attributes: { "points:version": "0.1" } },
        "shapes/cells": { attributes: {} },
        "images/dapi": { attributes: imageAttrs() },
      }),
    });

    expect(await imageUrlTargets(STORE)).toEqual([
      { url: `${STORE}/images/dapi`, name: "dapi" },
    ]);
  });

  it("is found through the folder that names the store", async () => {
    serve({
      [`${FOLDER}/mapping.json`]: { cells: "sample.zarr" },
      [`${STORE}/zarr.json`]: storeRoot({
        "images/dapi": { attributes: imageAttrs() },
      }),
    });

    expect(await imageUrlTargets(FOLDER)).toEqual([
      { url: `${STORE}/images/dapi`, name: "dapi" },
    ]);
  });

  it("is every image of a store, by name, for the user to choose", async () => {
    serve({
      [`${STORE}/zarr.json`]: storeRoot({
        "images/he": { attributes: imageAttrs() },
        "images/dapi": { attributes: imageAttrs() },
      }),
    });

    expect(await imageUrlTargets(STORE)).toEqual([
      { url: `${STORE}/images/dapi`, name: "dapi" },
      { url: `${STORE}/images/he`, name: "he" },
    ]);
  });

  it("reads an image whose NGFF metadata sits at the top of its attributes", async () => {
    serve({
      [`${STORE}/zarr.json`]: storeRoot({
        "images/dapi": { attributes: imageAttrs(true) },
      }),
    });

    expect(await imageUrlTargets(STORE)).toEqual([
      { url: `${STORE}/images/dapi`, name: "dapi" },
    ]);
  });

  it("is none in a store without images, or a URL that can't be read", async () => {
    serve({
      [`${STORE}/zarr.json`]: storeRoot({ "tables/table": { attributes: {} } }),
    });
    expect(await imageUrlTargets(STORE)).toEqual([]);

    serve({});
    expect(await imageUrlTargets(FOLDER)).toEqual([]);
  });

  it("does not mistake a labels element for an image", async () => {
    serve({
      [`${STORE}/zarr.json`]: storeRoot({
        "labels/segmentation": { attributes: imageAttrs() },
      }),
    });

    expect(await imageUrlTargets(STORE)).toEqual([]);
  });
});
