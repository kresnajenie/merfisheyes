/**
 * The SpatialData-DB catalog sync against the published contract: a release
 * file in their schema reads, a row becomes exactly one card, a file in a
 * different schema is refused, and a uid that left the release is unpublished.
 * Fixture: tests/data/sddb/catalog-fixture.parquet, the first rows of release
 * catalog-2026-10-20 (scripts/testdata/make-sddb-catalog-fixture.py).
 */
import { describe, expect, it } from "vitest";

import { readCatalogFile, CatalogSchemaError } from "@/lib/sddb/catalog-file";
import {
  formatBytes,
  mapCatalogRow,
  sampleRows,
  staleUids,
  storeUrlFromZarrUrl,
} from "@/lib/sddb/catalog-row";
import { READ_COLUMNS, REQUIRED_COLUMNS } from "@/lib/sddb/catalog-schema";

const FIXTURE = "tests/data/sddb/catalog-fixture.parquet";
const OTHER = "tests/data/sddb/collections-fixture.parquet";

describe("readCatalogFile", () => {
  it("reads a release file in the published schema", async () => {
    const file = await readCatalogFile(FIXTURE);

    expect(file.rows.length).toBe(3);
    expect(file.generatedAt).toBe("2026-10-20");
    expect(file.missingColumns).toEqual([]);
    for (const c of REQUIRED_COLUMNS) expect(file.rows[0][c]).toBeTruthy();
    expect(READ_COLUMNS.every((c) => c in file.rows[0])).toBe(true);
  });

  it("refuses a parquet file in another schema", async () => {
    await expect(readCatalogFile(OTHER)).rejects.toBeInstanceOf(
      CatalogSchemaError,
    );
  });
});

describe("storeUrlFromZarrUrl", () => {
  it("rewrites their bucket URL to the HTTPS form the browser reads", () => {
    expect(
      storeUrlFromZarrUrl(
        "s3://scverse-spatial-eu-central-1/.lamindb/hODmVoyNFahbWd1X.zarr",
      ),
    ).toBe(
      "https://scverse-spatial-eu-central-1.s3.eu-central-1.amazonaws.com/.lamindb/hODmVoyNFahbWd1X.zarr",
    );
  });

  it("refuses a bucket whose region it cannot tell", () => {
    expect(() => storeUrlFromZarrUrl("s3://some-bucket/x.zarr")).toThrow(
      /region/,
    );
    expect(() => storeUrlFromZarrUrl("https://x/y")).toThrow(/s3:/);
  });
});

describe("mapCatalogRow", () => {
  it("turns a release row into one card with one spatialdata entry", async () => {
    const { rows, generatedAt } = await readCatalogFile(FIXTURE);
    const m = mapCatalogRow(rows[0], {
      release: "catalog-2026-10-20",
      generatedAt,
    });

    expect(m.uid).toBe("hODmVoyNFahbWd1X0000");
    expect(m.data).toMatchObject({
      title: "GSM8641107",
      species: "Human",
      tissue: "Pancreas",
      platform: "Visium",
      source: "sddb",
      sddbUid: "hODmVoyNFahbWd1X0000",
      license: "unknown",
      licenseNoncommercial: false,
      sourceRelease: "catalog-2026-10-20",
      sourceGeneratedAt: "2026-10-20",
      numCells: 3571,
      numGenes: 19144,
      publicationLink: "https://doi.org/10.1016/j.celrep.2025.116827",
      isPublished: true,
      isBil: false,
      isCommunity: false,
    });
    expect(m.data.tags).toContain("spatialdata-db");
    expect(m.data.metadata.cardUrl).toMatch(/^https:\/\/lamin\.ai\//);
    expect(m.data.metadata.storeSize).toBe("79.2 MB");
    expect(m.entry).toEqual({
      label: "Visium store",
      datasetType: "spatialdata",
      s3BaseUrl:
        "https://scverse-spatial-eu-central-1.s3.eu-central-1.amazonaws.com/.lamindb/hODmVoyNFahbWd1X.zarr",
    });
  });

  it("maps their technology and organism spellings onto Explore's", () => {
    const base = {
      uid: "u",
      zarr_url: "s3://scverse-spatial-eu-central-1/.lamindb/u.zarr",
      organism: "mouse",
      assay: "Vizgen MERSCOPE",
      validation_status: "pass",
    };
    const m = mapCatalogRow(
      { ...base, technology: "Merscope", tissue: "brain", license: "CC-BY-4.0", license_noncommercial: false } as never,
      { release: "r", generatedAt: null },
    );

    expect(m.data.platform).toBe("MERSCOPE");
    expect(m.data.species).toBe("Mouse");
    expect(m.data.tissue).toBe("Brain");
    expect(m.data.title).toBe("u");
    expect(m.entry.label).toBe("MERSCOPE store");
    expect(
      mapCatalogRow({ ...base, technology: "VisiumHD" } as never, { release: "r", generatedAt: null }).data.platform,
    ).toBe("Visium HD");
  });

  it("refuses a row without a uid or a store", () => {
    expect(() =>
      mapCatalogRow({ technology: "Visium" } as never, { release: "r", generatedAt: null }),
    ).toThrow(/uid/);
    expect(() =>
      mapCatalogRow({ uid: "u", technology: "Visium" } as never, { release: "r", generatedAt: null }),
    ).toThrow(/zarr_url/);
  });
});

describe("staleUids", () => {
  it("names the synced uids the release no longer carries", () => {
    expect(staleUids(["a", "b", "c"], ["b", "c", "d"])).toEqual(["a"]);
    expect(staleUids([], ["x"])).toEqual([]);
  });
});

describe("sampleRows", () => {
  const rows = [
    { technology: "Visium", id: 1 },
    { technology: "Visium", id: 2 },
    { technology: "Visium", id: 3 },
    { technology: "Xenium", id: 4 },
    { technology: "Merscope", id: 5 },
    { technology: "Xenium", id: 6 },
  ];

  it("spreads a small sample over the technologies, in file order", () => {
    expect(sampleRows(rows, 4).map((r) => r.id)).toEqual([1, 4, 5, 2]);
    expect(sampleRows(rows, 2).map((r) => r.id)).toEqual([1, 4]);
  });

  it("returns every row when asked for more than there are", () => {
    expect(sampleRows(rows, 10)).toHaveLength(6);
    expect(sampleRows([], 3)).toEqual([]);
  });
});

describe("formatBytes", () => {
  it("rounds to a readable unit", () => {
    expect(formatBytes(79186969)).toBe("79.2 MB");
    expect(formatBytes(3.1e12)).toBe("3.1 TB");
    expect(formatBytes(12_000)).toBe("12 KB");
  });
});
