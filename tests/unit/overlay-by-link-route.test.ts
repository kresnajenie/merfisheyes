import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * POST /api/datasets/[id]/overlay with { url }: an owner overlays public
 * molecules on their cells by pasting the link, without making a dataset of
 * them. The link is checked where it lives, as the viewer will open it.
 */

const auth = vi.fn();
const ownedDataset = vi.fn();
const getOverlaySmId = vi.fn();
const getOverlayUrl = vi.fn();
const writeOverlayUrl = vi.fn();
const writeOverlay = vi.fn();
const moleculeUrlFormat = vi.fn();

vi.mock("@/lib/auth", () => ({ auth }));
vi.mock("@/lib/ingest/overlay", () => ({
  ownedDataset,
  getOverlaySmId,
  getOverlayUrl,
  writeOverlayUrl,
  writeOverlay,
  hasMolecules: vi.fn(),
  getOverlayImage: vi.fn(),
  writeOverlayImage: vi.fn(),
  imageOverlayFor: vi.fn(),
  removeOverlay: vi.fn(),
  groupInProject: vi.fn(),
}));
vi.mock("@/lib/overlays/moleculeSource", async (actual) => ({
  ...(await actual<typeof import("@/lib/overlays/moleculeSource")>()),
  moleculeUrlFormat,
}));

const { POST } = await import("@/app/api/datasets/[datasetId]/overlay/route");

const FOLDER = "https://bucket.example/spatialdata/all_sets";

const post = (body: unknown) =>
  POST(
    new Request("http://localhost/api/datasets/ds_cells/overlay", {
      method: "POST",
      body: JSON.stringify(body),
    }) as never,
    { params: Promise.resolve({ datasetId: "ds_cells" }) },
  );

describe("overlaying molecules by their link", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    auth.mockResolvedValue({ user: { id: "u1", role: "USER" } });
    ownedDataset.mockResolvedValue({
      id: "ds_cells",
      ownerId: "u1",
      datasetType: "spatialdata",
      title: "cells",
    });
    getOverlaySmId.mockResolvedValue(null);
    getOverlayUrl.mockResolvedValue(null);
  });

  it("writes the link of a store with molecules", async () => {
    moleculeUrlFormat.mockResolvedValue("store");

    const res = await post({ url: `${FOLDER}/` });

    expect(res.status).toBe(200);
    expect(writeOverlayUrl).toHaveBeenCalledWith("ds_cells", FOLDER, "store");
    expect(writeOverlay).not.toHaveBeenCalled();
  });

  it("takes the data URL out of a merfisheyes viewer link", async () => {
    moleculeUrlFormat.mockResolvedValue("molecules");

    await post({
      url: `https://www.merfisheyes.com/sm-viewer/from-s3?url=${encodeURIComponent(FOLDER)}`,
    });

    expect(moleculeUrlFormat).toHaveBeenCalledWith(FOLDER);
    expect(writeOverlayUrl).toHaveBeenCalledWith("ds_cells", FOLDER, "molecules");
  });

  it("refuses a link with no molecules it can draw, and writes nothing", async () => {
    moleculeUrlFormat.mockResolvedValue(null);

    const res = await post({ url: FOLDER });

    expect(res.status).toBe(400);
    expect((await res.json()).error).toMatch(/No molecules at that link/);
    expect(writeOverlayUrl).not.toHaveBeenCalled();
  });

  it("asks before replacing an overlay, by dataset or by link", async () => {
    moleculeUrlFormat.mockResolvedValue("store");
    getOverlaySmId.mockResolvedValue("sm_old");

    expect((await post({ url: FOLDER })).status).toBe(409);
    expect(writeOverlayUrl).not.toHaveBeenCalled();

    expect((await post({ url: FOLDER, force: true })).status).toBe(200);
    expect(writeOverlayUrl).toHaveBeenCalledTimes(1);
  });

  it("only on the user's own cells", async () => {
    ownedDataset.mockResolvedValue(null);
    moleculeUrlFormat.mockResolvedValue("store");

    expect((await post({ url: FOLDER })).status).toBe(404);
    expect(writeOverlayUrl).not.toHaveBeenCalled();
  });

  it("needs a dataset or a link", async () => {
    expect((await post({})).status).toBe(400);
    expect((await post({ url: "not a link" })).status).toBe(400);
  });
});
