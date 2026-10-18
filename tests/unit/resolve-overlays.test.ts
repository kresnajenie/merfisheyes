/**
 * What is drawn with a cell dataset besides its cells, resolved with fakes
 * for the network: a dropped folder's parts, the hardcoded registry, a
 * mapping.json with an image, with an `__all__` link, with a per-cell link
 * column, and the app's own mapping route.
 */
import type { OverlayEvents, OverlayIO } from "@/lib/viewer/resolve-overlays";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { describe, expect, it, vi } from "vitest";

import { resolveOverlays } from "@/lib/viewer/resolve-overlays";

const image = { name: "dapi", levels: [] } as any;
const molecules = (id: string) =>
  ({ id, uniqueGenes: ["A"], metadata: {} }) as any;

function dataset(extra: Partial<StandardizedDataset> = {}): StandardizedDataset {
  return {
    id: "ds",
    normalized: false,
    metadata: {},
    clusters: [],
    allClusterColumnNames: ["batch"],
    ...extra,
  } as unknown as StandardizedDataset;
}

function fakes(mapping: any = null, appMapping: any = null) {
  const io: OverlayIO = {
    fetchMapping: vi.fn(async () => mapping),
    fetchAppMapping: vi.fn(async () => appMapping),
    openImage: vi.fn(async () => image),
    loadMolecules: vi.fn(async (ref: string) => molecules(ref)),
    loadColumn: vi.fn(async () => true),
  };
  const got: Record<string, unknown[]> = {};
  const on: OverlayEvents = {
    element: (p) => (got.element ??= []).push(p),
    link: (c) => (got.link ??= []).push(c),
    loading: (b) => (got.loading ??= []).push(b),
    source: (u) => (got.source ??= []).push(u),
    column: () => (got.column ??= []).push(true),
    warn: (m) => (got.warn ??= []).push(m),
    stale: () => false,
  };

  return { io, on, got };
}

const settle = () => new Promise((r) => setTimeout(r, 0));

describe("resolveOverlays", () => {
  it("hands over the parts dropped beside the cells, and skips an image under normalized cells", async () => {
    const { io, on, got } = fakes();
    const m = molecules("local");

    await resolveOverlays(
      dataset(),
      { image, shapes: null, molecules: m },
      null,
      io,
      on,
    );
    expect(got.element).toEqual([
      { image, shapes: undefined, molecules: m },
    ]);
    expect(got.link).toEqual([null]);
    expect(io.fetchMapping).not.toHaveBeenCalled();

    const n = fakes();

    await resolveOverlays(
      dataset({ normalized: true } as any),
      { image, shapes: null, molecules: null },
      null,
      n.io,
      n.on,
    );
    expect(n.got.element).toEqual([{ image: undefined, shapes: undefined }]);
    expect(n.got.warn?.[0]).toMatch(/not drawn/);
  });

  it("takes the hardcoded registry as the link and looks no further", async () => {
    const { io, on, got } = fakes();
    const registry = { linkColumn: "batch", links: { a: "u" } };

    await resolveOverlays(dataset(), null, registry, io, on);
    expect(got.link).toEqual([registry]);
    expect(got.element).toBeUndefined();
    expect(io.fetchMapping).not.toHaveBeenCalled();
  });

  it("draws the image a public folder's mapping names, with no link", async () => {
    const { io, on, got } = fakes({ image: "https://x/dapi.zarr" });

    await resolveOverlays(
      dataset({ metadata: { customS3BaseUrl: "https://x" } } as any),
      null,
      null,
      io,
      on,
    );
    await settle();
    expect(got.element).toEqual([{ image }]);
    expect(got.link).toEqual([null]);
  });

  it("loads the molecules of an __all__ link and says where they came from", async () => {
    const mapping = { linkColumn: "__all__", links: { __all__: "https://x/sm" } };
    const { io, on, got } = fakes(mapping);

    await resolveOverlays(
      dataset({ metadata: { customS3BaseUrl: "https://x" } } as any),
      null,
      null,
      io,
      on,
    );
    expect(got.link).toEqual([null, mapping]);
    expect(got.loading).toEqual([true, false]);
    expect(got.element).toEqual([{ molecules: molecules("https://x/sm") }]);
    expect(got.source).toEqual(["https://x/sm"]);
    expect(io.loadMolecules).toHaveBeenCalledWith("https://x/sm", "s3");
  });

  it("fetches a per-cell link column so the split opens at once", async () => {
    const mapping = { linkColumn: "batch", links: { a: "u" } };
    const { io, on, got } = fakes(mapping);

    await resolveOverlays(
      dataset({ metadata: { customS3BaseUrl: "https://x" } } as any),
      null,
      null,
      io,
      on,
    );
    expect(got.link).toEqual([null, mapping]);
    expect(io.loadColumn).toHaveBeenCalledWith(expect.anything(), "batch");
    expect(got.column).toEqual([true]);
    expect(io.loadMolecules).not.toHaveBeenCalled();

    // A column the dataset does not have is left alone
    const missing = fakes({ linkColumn: "nope", links: { a: "u" } });

    await resolveOverlays(
      dataset({ metadata: { customS3BaseUrl: "https://x" } } as any),
      null,
      null,
      missing.io,
      missing.on,
    );
    expect(missing.io.loadColumn).not.toHaveBeenCalled();
  });

  it("asks the app's mapping route for an uploaded dataset", async () => {
    const mapping = {
      linkColumn: "__all__",
      links: { __all__: "sm_1" },
      source: "app",
      image: "img_1",
      imageSource: "app",
    };
    const { io, on, got } = fakes(null, mapping);

    await resolveOverlays(dataset(), null, null, io, on);
    await settle();
    expect(io.fetchAppMapping).toHaveBeenCalledWith("ds");
    expect(io.loadMolecules).toHaveBeenCalledWith("sm_1", "app");
    expect(io.openImage).toHaveBeenCalledWith({ image: "img_1", source: "app" });
    expect(got.element).toEqual(
      expect.arrayContaining([{ molecules: molecules("sm_1") }, { image }]),
    );
    expect(got.source).toBeUndefined();
  });

  it("hands nothing over once the caller has moved on", async () => {
    const mapping = { linkColumn: "__all__", links: { __all__: "https://x/sm" } };
    const { io, on, got } = fakes(mapping);

    on.stale = () => true;
    await resolveOverlays(
      dataset({ metadata: { customS3BaseUrl: "https://x" } } as any),
      null,
      null,
      io,
      on,
    );
    expect(got.element).toBeUndefined();
  });
});
