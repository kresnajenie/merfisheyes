import { afterEach, describe, expect, it, vi } from "vitest";

import {
  datasetKind,
  isAppDatasetId,
  viewerBase,
  viewerHref,
} from "@/lib/datasets/kind";
import { isOmeZarrImage } from "@/lib/upload/uploadImageToS3";
import { fetchMappingConfig } from "@/lib/config/dataset-links";

describe("dataset kinds", () => {
  it("treats legacy and missing types as single cell", () => {
    // Older rows hold the raw format, or nothing at all
    for (const type of [
      "single_cell",
      "h5ad",
      "merscope",
      "",
      null,
      undefined,
    ]) {
      expect(datasetKind(type)).toBe("single_cell");
    }
    expect(datasetKind("single_molecule")).toBe("single_molecule");
    expect(datasetKind("image")).toBe("image");
  });

  it("opens each kind in its own viewer", () => {
    expect(viewerBase("h5ad")).toBe("/viewer");
    expect(viewerBase("single_molecule")).toBe("/sm-viewer");
    expect(viewerBase("image")).toBe("/image-viewer");
  });

  it("opens by URL when registered from one, else by id", () => {
    expect(
      viewerHref({
        id: "s3_abc",
        datasetType: "image",
        s3BaseUrl: "https://bucket.example/folder/dapi.zarr/",
      }),
    ).toBe(
      "/image-viewer/from-s3?url=https%3A%2F%2Fbucket.example%2Ffolder%2Fdapi.zarr",
    );
    expect(viewerHref({ id: "ds_abc", datasetType: "image" })).toBe(
      "/image-viewer/ds_abc",
    );
    expect(viewerHref({ datasetType: "image" })).toBeNull();
  });

  it("an app dataset is ds_ or sm_, whatever it became", () => {
    // A molecule table from the single molecule card, made a SpatialData
    // store by prepare.py, keeps its sm_ id and opens in that viewer
    expect(isAppDatasetId("sm_FtBRk77t3F")).toBe(true);
    expect(viewerHref({ id: "sm_FtBRk77t3F", datasetType: "spatialdata" })).toBe(
      "/spatialdata-viewer/sm_FtBRk77t3F",
    );
    expect(isAppDatasetId("ds_abc")).toBe(true);
    expect(isAppDatasetId("spatialdata_x")).toBe(false);
    expect(isAppDatasetId("../ds_abc")).toBe(false);
    expect(isAppDatasetId(undefined)).toBe(false);
  });
});

describe("isOmeZarrImage", () => {
  const withAxes = (axes: string[]) => ({
    zarr_format: 3,
    node_type: "group",
    attributes: {
      ome: { multiscales: [{ axes: axes.map((name) => ({ name })) }] },
    },
  });

  it("accepts a multiscale image with axes c, y, x", () => {
    expect(isOmeZarrImage(withAxes(["c", "y", "x"]))).toBe(true);
  });

  it("rejects other zarr stores", () => {
    // What the viewer can't draw must not be offered for upload as an image
    expect(isOmeZarrImage(withAxes(["z", "y", "x"]))).toBe(false);
    expect(
      isOmeZarrImage({ zarr_format: 3, node_type: "group", attributes: {} }),
    ).toBe(false);
    expect(isOmeZarrImage(null)).toBe(false);
  });
});

describe("fetchMappingConfig", () => {
  const dataset = {
    metadata: { customS3BaseUrl: "https://bucket.example/folder" },
  };
  const serve = (body: unknown) =>
    vi.stubGlobal("fetch", async () => new Response(JSON.stringify(body)));

  afterEach(() => vi.unstubAllGlobals());

  it("accepts a mapping that only names an image", async () => {
    // Beside a dataset on a public bucket, with no molecules to link
    serve({ image: "https://bucket.example/folder/dapi.zarr" });
    expect(await fetchMappingConfig(dataset)).toEqual({
      image: "https://bucket.example/folder/dapi.zarr",
    });
  });

  it("keeps the image of a mapping that also links molecules", async () => {
    const mapping = {
      linkColumn: "__all__",
      links: { __all__: "https://bucket.example/folder/molecules" },
      image: "https://bucket.example/folder/dapi.zarr",
    };

    serve(mapping);
    expect(await fetchMappingConfig(dataset)).toEqual(mapping);
  });

  it("ignores a file that is neither", async () => {
    serve({ something: "else" });
    expect(await fetchMappingConfig(dataset)).toBeNull();
  });
});
