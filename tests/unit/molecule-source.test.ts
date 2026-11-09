import { afterEach, describe, expect, it, vi } from "vitest";

import { dataUrlOf, moleculeUrlFormat } from "@/lib/overlays/moleculeSource";

/**
 * What a pasted link holds as molecules (lib/overlays/moleculeSource.ts), the
 * one rule the overlay dialog, its route and a local preview's check share.
 */

const FOLDER = "https://bucket.example/spatialdata/all_sets";
const STORE = `${FOLDER}/all_sets_2d.zarr`;

const root = (nodes: Record<string, unknown>, ome = false) => ({
  attributes: ome
    ? { ome: { multiscales: [{}] } }
    : { spatialdata_attrs: { version: "0.2" } },
  consolidated_metadata: { metadata: nodes },
});
const table = { "tables/table": { attributes: {} } };
const points = { "points/transcripts": { attributes: { "points:version": "0.1" } } };

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

describe("the data URL a pasted link names", () => {
  it("is the link itself, without a trailing slash", () => {
    expect(dataUrlOf(` ${FOLDER}/ `)).toBe(FOLDER);
  });

  it("is the URL a merfisheyes viewer link opened from S3", () => {
    const viewer = `https://www.merfisheyes.com/spatialdata-viewer/from-s3?url=${encodeURIComponent(FOLDER)}`;

    expect(dataUrlOf(viewer)).toBe(FOLDER);
  });

  it("is nothing for text that is not an http(s) URL", () => {
    expect(dataUrlOf("all_sets")).toBeNull();
    expect(dataUrlOf("s3://bucket/all_sets")).toBeNull();
  });
});

describe("the molecules at a public URL", () => {
  it("are a store's points, through the folder that names it", async () => {
    serve({
      [`${FOLDER}/mapping.json`]: { cells: "all_sets_2d.zarr" },
      [`${STORE}/zarr.json`]: root({ ...table, ...points }),
    });

    expect(await moleculeUrlFormat(FOLDER)).toBe("store");
    expect(await moleculeUrlFormat(STORE)).toBe("store");
  });

  it("are a merfisheyes molecule folder, by its manifest", async () => {
    serve({ [`${FOLDER}/manifest.json.gz`]: {} });

    expect(await moleculeUrlFormat(FOLDER)).toBe("molecules");
  });

  it("are none in a store of cells alone, an image, or a URL that can't be read", async () => {
    serve({ [`${STORE}/zarr.json`]: root(table) });
    expect(await moleculeUrlFormat(STORE)).toBeNull();

    serve({ [`${STORE}/zarr.json`]: root({}, true) });
    expect(await moleculeUrlFormat(STORE)).toBeNull();

    serve({});
    expect(await moleculeUrlFormat(FOLDER)).toBeNull();
  });
});
