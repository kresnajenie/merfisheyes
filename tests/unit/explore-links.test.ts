/**
 * Every catalog row has a standalone page and a card links to it; an entry on
 * that page opens the viewer of its kind. One module decides both, so the
 * grid, the carousel, the detail page and the info modal agree.
 */
import { describe, expect, it } from "vitest";

import { catalogPageHref, entryHref } from "@/lib/explore/links";

describe("catalogPageHref", () => {
  it("sends a BIL row to its BIL page", () => {
    expect(catalogPageHref({ id: "c1", bilCode: "ace-bag-tin" })).toBe(
      "/explore/bil/ace-bag-tin",
    );
  });

  it("sends every other row to /explore/[id], whatever its entries", () => {
    expect(catalogPageHref({ id: "c1", bilCode: null })).toBe("/explore/c1");
    expect(catalogPageHref({ id: "c2" })).toBe("/explore/c2");
  });
});

describe("entryHref", () => {
  const folder = "https://bucket.s3.us-west-2.amazonaws.com/spatialdata/skin";

  it("opens a SpatialData entry in the SpatialData viewer by URL", () => {
    expect(
      entryHref({
        datasetType: "spatialdata",
        s3BaseUrl: `${folder}/`,
        datasetId: "s3_abc",
      }),
    ).toBe(`/spatialdata-viewer/from-s3?url=${encodeURIComponent(folder)}`);
  });

  it("opens cells, molecules and images in their own viewers", () => {
    expect(
      entryHref({ datasetType: "single_cell", s3BaseUrl: folder, datasetId: null }),
    ).toBe(`/viewer/from-s3?url=${encodeURIComponent(folder)}`);
    expect(
      entryHref({ datasetType: "single_molecule", s3BaseUrl: null, datasetId: "sm_1" }),
    ).toBe("/sm-viewer/sm_1");
    expect(
      entryHref({ datasetType: "image", s3BaseUrl: null, datasetId: "ds_img" }),
    ).toBe("/image-viewer/ds_img");
  });

  it("prefers the public URL over the id and gives null with neither", () => {
    expect(
      entryHref({ datasetType: "single_cell", s3BaseUrl: folder, datasetId: "s3_x" }),
    ).toBe(`/viewer/from-s3?url=${encodeURIComponent(folder)}`);
    expect(
      entryHref({ datasetType: "single_cell", s3BaseUrl: null, datasetId: null }),
    ).toBeNull();
  });
});
