/**
 * mapping.json names its parts relative to its own folder
 * (scripts/prepare_folder.py writes "molecules", "dapi.zarr"); the viewer
 * turns them into URLs against the folder it fetched the mapping from.
 */
import { afterEach, describe, expect, it, vi } from "vitest";

import {
  fetchMappingConfig,
  resolveMappingPaths,
} from "@/lib/config/dataset-links";

const FOLDER = "https://bucket.s3.eu-west-1.amazonaws.com/data/sample";

describe("resolveMappingPaths", () => {
  it("resolves names against the folder and leaves URLs alone", () => {
    expect(
      resolveMappingPaths(
        {
          linkColumn: "__all__",
          links: { __all__: "molecules", other: "https://x/y" },
          image: "./dapi.zarr",
        },
        `${FOLDER}/`,
      ),
    ).toEqual({
      linkColumn: "__all__",
      links: { __all__: `${FOLDER}/molecules`, other: "https://x/y" },
      image: `${FOLDER}/dapi.zarr`,
    });
  });

  it("leaves app ids alone", () => {
    const mapping = {
      linkColumn: "__all__",
      links: { __all__: "sm_abc123" },
      source: "app" as const,
      image: "ds_img456",
      imageSource: "app" as const,
    };

    expect(resolveMappingPaths(mapping, FOLDER)).toEqual(mapping);
  });
});

describe("fetchMappingConfig", () => {
  afterEach(() => vi.unstubAllGlobals());

  it("returns the mapping with its names resolved", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => ({
        ok: true,
        json: async () => ({
          cells: "sample.zarr",
          linkColumn: "__all__",
          links: { __all__: "molecules" },
          image: "dapi.zarr",
        }),
      })),
    );
    const mapping = await fetchMappingConfig({
      metadata: { customS3BaseUrl: FOLDER },
    });

    expect(mapping).toMatchObject({
      links: { __all__: `${FOLDER}/molecules` },
      image: `${FOLDER}/dapi.zarr`,
    });
  });

  it("an image alone is resolved too", async () => {
    vi.stubGlobal(
      "fetch",
      vi.fn(async () => ({
        ok: true,
        json: async () => ({ cells: "sample.zarr", image: "dapi.zarr" }),
      })),
    );
    expect(
      await fetchMappingConfig({ metadata: { customS3BaseUrl: FOLDER } }),
    ).toEqual({ image: `${FOLDER}/dapi.zarr`, imageSource: undefined });
  });
});
