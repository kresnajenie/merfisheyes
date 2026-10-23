import { afterEach, describe, expect, it, vi } from "vitest";

import { compose } from "@/lib/spatialdata/coordinates";
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

  it("follows a sequence of transforms, as spatialdata writes them", async () => {
    // Levels in pixels of the full image; the image placed by scale then shift
    const root = group(["c", "y", "x"]);
    const multiscale: any = root.attributes.ome.multiscales[0];

    multiscale.datasets = [0, 1].map((l) => ({
      path: `s${l}`,
      coordinateTransformations: [
        { type: "scale", scale: [1, 2 ** l, 2 ** l] },
        { type: "translation", translation: [0, l * 0.5, l * 0.5] },
      ],
    }));
    multiscale.coordinateTransformations = [
      {
        type: "sequence",
        transformations: [
          { type: "scale", scale: [1, 0.1, 0.1] },
          { type: "translation", translation: [0, -11, -12] },
        ],
      },
    ];
    serve({
      "/zarr.json": root,
      "/s0/zarr.json": level(8, 12),
      "/s1/zarr.json": level(4, 6),
    });

    const { levels, placement } = await openOmeZarrImage(BASE);

    // The levels in the image's own units; the sequence is the placement
    expect(levels[0].pixelWidth).toBeCloseTo(1);
    expect(levels[0].originX).toBeCloseTo(0);
    expect(levels[1].pixelWidth).toBeCloseTo(2);
    expect(levels[1].originX).toBeCloseTo(0.5);
    expect(placement.scale[0]).toBeCloseTo(0.1);
    expect(placement.offset).toEqual([-12, -11]);
    // Placed: pixel (0, 0) of the full image at (-12, -11), each 0.1 wide
    const placed = compose(
      { scale: [levels[1].pixelWidth, levels[1].pixelHeight], offset: [levels[1].originX, levels[1].originY] },
      placement,
    );

    expect(placed.scale[0]).toBeCloseTo(0.2);
    expect(placed.offset[0]).toBeCloseTo(-11.95);
  });

  it("rejects axes other than c, y, x", async () => {
    serve({ "/zarr.json": group(["c", "z", "y", "x"]) });

    await expect(openOmeZarrImage(BASE)).rejects.toThrow(
      /Unsupported image axes/,
    );
  });
});
