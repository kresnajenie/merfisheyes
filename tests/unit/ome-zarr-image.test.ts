import { afterEach, describe, expect, it, vi } from "vitest";

import { openOmeZarrImage } from "@/lib/spatialdata/image";

const BASE = "https://bucket.example/folder/dapi.zarr";

const level = (h: number, w: number) => ({
  zarr_format: 3,
  node_type: "array",
  shape: [1, h, w],
  data_type: "uint16",
  chunk_grid: { name: "regular", configuration: { chunk_shape: [1, 4, 4] } },
  chunk_key_encoding: { name: "default", configuration: { separator: "/" } },
  fill_value: 0,
  codecs: [{ name: "bytes", configuration: { endian: "little" } }],
  attributes: {},
});

const group = (axes: string[]) => ({
  zarr_format: 3,
  node_type: "group",
  attributes: {
    ome: {
      version: "0.5",
      multiscales: [
        {
          axes: axes.map((name) => ({ name })),
          datasets: [0, 1].map((l) => ({
            path: String(l),
            coordinateTransformations: [
              { type: "scale", scale: [1, 0.5 * 2 ** l, 0.5 * 2 ** l] },
              { type: "translation", translation: [0, 20, 10] },
            ],
          })),
        },
      ],
      omero: { channels: [{ label: "DAPI" }] },
    },
  },
});

/** Serve zarr.json documents by path; anything else is a missing key. */
function serve(files: Record<string, unknown>) {
  vi.stubGlobal("fetch", async (url: string) => {
    const body = files[url.slice(BASE.length)];

    // A public bucket that can't be listed answers 403 for a missing key
    return body
      ? new Response(JSON.stringify(body))
      : new Response(null, { status: 403 });
  });
}

afterEach(() => vi.unstubAllGlobals());

describe("openOmeZarrImage", () => {
  it("reads a standalone multiscale image and places it by its own transforms", async () => {
    serve({
      "/zarr.json": group(["c", "y", "x"]),
      "/0/zarr.json": level(8, 12),
      "/1/zarr.json": level(4, 6),
    });

    const image = await openOmeZarrImage(`${BASE}/`);

    expect(image.name).toBe("dapi");
    expect(image.channelLabels).toEqual(["DAPI"]);
    expect(image.maxValue).toBe(65535);
    expect(image.levels).toHaveLength(2);
    expect(image.levels[0]).toMatchObject({
      width: 12,
      height: 8,
      chunkWidth: 4,
      chunkHeight: 4,
      pixelWidth: 0.5,
      pixelHeight: 0.5,
      originX: 10,
      originY: 20,
    });
    expect(image.levels[1]).toMatchObject({ pixelWidth: 1, pixelHeight: 1 });
  });

  it("rejects axes other than c, y, x", async () => {
    serve({ "/zarr.json": group(["c", "z", "y", "x"]) });

    await expect(openOmeZarrImage(BASE)).rejects.toThrow(
      /Unsupported image axes/,
    );
  });
});
