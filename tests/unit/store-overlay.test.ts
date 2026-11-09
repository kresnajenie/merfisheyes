import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

/**
 * A store that holds cells AND molecules can lend just its molecules as the
 * overlay of another dataset, whether it was uploaded or registered by its
 * public URL (an Explore entry). The route accepts it when storeHasPoints
 * says it has molecules, and writeOverlay links it by URL as a store.
 */

const findUnique = vi.fn();
const listObjectKeys = vi.fn();
const getObjectJson = vi.fn();
const uploadBufferToS3 = vi.fn();

vi.mock("@/lib/prisma", () => ({ prisma: { dataset: { findUnique } } }));
vi.mock("@/lib/s3", () => ({
  listObjectKeys,
  getObjectJson,
  uploadBufferToS3,
  deleteObject: vi.fn(),
}));

const { rootHasPoints, storeHasPoints, writeOverlay } = await import(
  "@/lib/ingest/overlay"
);

const FOLDER = "https://bucket.example/spatialdata/all_sets";
const STORE = `${FOLDER}/all_sets_2d.zarr`;

/** A store root as spatialdata writes it, with these elements. */
const root = (nodes: Record<string, unknown>) => ({
  attributes: { spatialdata_attrs: { version: "0.2" } },
  consolidated_metadata: { metadata: nodes },
});
const table = { "tables/table": { attributes: { "encoding-type": "anndata" } } };
const points = {
  "points/transcripts": { attributes: { "points:version": "0.1" } },
};

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

beforeEach(() => vi.clearAllMocks());
afterEach(() => vi.unstubAllGlobals());

describe("a store's molecules, from its metadata", () => {
  it("counts a points element in the format the viewer draws", () => {
    expect(rootHasPoints(root({ ...table, ...points }))).toBe(true);
  });

  it("does not count a table alone, nor spatialdata's own parquet points", () => {
    expect(rootHasPoints(root(table))).toBe(false);
    expect(
      rootHasPoints(root({ "points/transcripts": { attributes: {} } })),
    ).toBe(false);
    expect(rootHasPoints(null)).toBe(false);
  });
});

describe("storeHasPoints", () => {
  it("reads a store registered by its folder's URL where it lives", async () => {
    findUnique.mockResolvedValue({
      s3BaseUrl: FOLDER,
      ingestSource: "s3_registered",
    });
    serve({
      [`${FOLDER}/mapping.json`]: { cells: "all_sets_2d.zarr" },
      [`${STORE}/zarr.json`]: root({ ...table, ...points }),
    });

    expect(await storeHasPoints("s3_3k27Qke0ET")).toBe(true);
    // Not in the app's bucket: listing it would find nothing
    expect(listObjectKeys).not.toHaveBeenCalled();
  });

  it("says no for a registered store with cells alone", async () => {
    findUnique.mockResolvedValue({ s3BaseUrl: STORE, ingestSource: null });
    serve({ [`${STORE}/zarr.json`]: root(table) });

    expect(await storeHasPoints("s3_cellsonly")).toBe(false);
  });

  it("looks in the app's bucket for an uploaded store, as before", async () => {
    findUnique.mockResolvedValue({ s3BaseUrl: null, ingestSource: "upload" });
    listObjectKeys.mockResolvedValue([
      { key: "datasets/ds_1/ds_1.zarr/points/molecules/zarr.json" },
    ]);

    expect(await storeHasPoints("ds_1")).toBe(true);
    expect(listObjectKeys).toHaveBeenCalledWith("datasets/ds_1/");
  });
});

describe("writeOverlay", () => {
  it("links a registered store by its URL, as a store", async () => {
    findUnique.mockResolvedValue({
      id: "s3_3k27Qke0ET",
      s3BaseUrl: FOLDER,
      ingestSource: "s3_registered",
      datasetType: "spatialdata",
    });
    getObjectJson.mockResolvedValue(null);

    await writeOverlay("ds_cells", "s3_3k27Qke0ET");

    const [key, body] = uploadBufferToS3.mock.calls[0];

    expect(key).toBe("datasets/ds_cells/mapping.json");
    expect(JSON.parse(body.toString())).toEqual({
      linkColumn: "__all__",
      links: { __all__: FOLDER },
      source: "s3",
      datasetId: "s3_3k27Qke0ET",
      linkFormat: "store",
    });
  });
});
