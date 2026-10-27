import { afterEach, describe, expect, it, vi } from "vitest";

import {
  applyImageState,
  encodeImageState,
} from "@/lib/spatialdata/imageUrlState";
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";

const STORE = { attributes: { spatialdata_attrs: { version: "0.2" } } };
const IMAGE = { attributes: { ome: { multiscales: [{}] } } };

/** A bucket that can't be listed: only these URLs answer. */
function bucket(objects: Record<string, unknown>) {
  vi.stubGlobal(
    "fetch",
    vi.fn(async (url: string) =>
      url in objects
        ? new Response(JSON.stringify(objects[url]))
        : new Response("", { status: 403 }),
    ),
  );
}

afterEach(() => vi.unstubAllGlobals());

describe("resolveZarrUrl", () => {
  const base = "https://bucket.example/data/skin";

  it("takes a store or an image as it is", async () => {
    bucket({
      [`${base}/skin.zarr/zarr.json`]: STORE,
      [`${base}/dapi.zarr/zarr.json`]: IMAGE,
    });

    // The root read to tell what it is comes back, so it is not read twice
    expect(await resolveZarrUrl(`${base}/skin.zarr/`)).toEqual({
      kind: "spatialdata",
      url: `${base}/skin.zarr`,
      root: STORE,
    });
    expect(await resolveZarrUrl(`${base}/dapi.zarr`)).toEqual({
      kind: "image",
      url: `${base}/dapi.zarr`,
      root: IMAGE,
    });
  });

  it("finds a folder's store by the name in its mapping.json", async () => {
    bucket({
      [`${base}/mapping.json`]: { cells: "cells_v2.zarr" },
      [`${base}/cells_v2.zarr/zarr.json`]: STORE,
    });

    expect(await resolveZarrUrl(`${base}/`)).toEqual({
      kind: "spatialdata",
      url: `${base}/cells_v2.zarr`,
      root: STORE,
    });
  });

  it("accepts a full URL there", async () => {
    const elsewhere = "https://other.example/stores/a.zarr";

    bucket({
      [`${base}/mapping.json`]: { cells: `${elsewhere}/` },
      [`${elsewhere}/zarr.json`]: STORE,
    });

    expect((await resolveZarrUrl(base))?.url).toBe(elsewhere);
  });

  it("otherwise tries a store named after the folder", async () => {
    bucket({ [`${base}/skin.zarr/zarr.json`]: STORE });

    expect((await resolveZarrUrl(base))?.url).toBe(`${base}/skin.zarr`);
  });

  it("gives up on a folder with no store it can name", async () => {
    bucket({ [`${base}/other.zarr/zarr.json`]: STORE });

    expect(await resolveZarrUrl(base)).toBeNull();
  });
});

describe("image settings in a link", () => {
  /** The part of an ImageLayer the link reads and writes. */
  function layer(windowed: boolean[]) {
    const channels = windowed.map((_, i) => ({
      label: `ch${i}`,
      color: "#3d7dff",
      visible: i === 0,
      min: 0,
      max: 255,
    }));

    return {
      channels,
      maxValue: 255,
      visible: true,
      hasWindow: (i: number) => windowed[i],
      setChannel(i: number, patch: object) {
        if ("min" in patch) windowed[i] = true;
        Object.assign(channels[i], patch);
      },
    } as any;
  }

  it("round-trips visibility, colour and window", () => {
    const from = layer([true, false]);

    from.visible = false;
    from.setChannel(0, { color: "#ff8800", min: 10, max: 120 });
    from.setChannel(1, { visible: true });

    const to = layer([false, false]);

    expect(applyImageState(to, encodeImageState(from))).toBe(true);
    expect(to.visible).toBe(false);
    expect(to.channels[0]).toMatchObject({
      color: "#ff8800",
      min: 10,
      max: 120,
    });
    expect(to.channels[1].visible).toBe(true);
    // A channel that had no window yet still gets its own from the data
    expect(to.hasWindow(1)).toBe(false);
  });

  it("round-trips a narrowed slider range", () => {
    const from = layer([true]);

    from.setChannel(0, { min: 10, max: 120, range: [5, 150] });
    const to = layer([false]);

    expect(applyImageState(to, encodeImageState(from))).toBe(true);
    expect(to.channels[0]).toMatchObject({ min: 10, max: 120, range: [5, 150] });
  });

  it("ignores a range that does not hold its window", () => {
    const from = layer([true]);

    // Written by hand: a window of 10–120 inside a range of 20–100
    from.setChannel(0, { min: 10, max: 120, range: [20, 100] });
    const to = layer([false]);

    expect(applyImageState(to, encodeImageState(from))).toBe(false);
    expect(to.channels[0].range).toBeUndefined();
  });

  it("ignores settings written for another image", () => {
    const to = layer([false]);

    expect(applyImageState(to, encodeImageState(layer([true, true])))).toBe(
      false,
    );
    expect(applyImageState(to, "not-a-state")).toBe(false);
    expect(to.channels[0]).toMatchObject({ visible: true, max: 255 });
  });
});
