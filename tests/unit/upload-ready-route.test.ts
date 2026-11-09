import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * POST /api/uploads/ready: a save uploaded as it is (an image, a SpatialData
 * store, a processed folder) never reaches the worker, so the browser says
 * when it is done and the signed-in user gets one "ready" email for the whole
 * save (docs/plans/PREVIEW-SAVE.md, D6).
 */

const requireUser = vi.fn();
const findMany = vi.fn();
const sendDatasetReadyEmail = vi.fn();

vi.mock("@/lib/admin-auth", () => ({ requireUser }));
vi.mock("@/lib/prisma", () => ({ prisma: { dataset: { findMany } } }));
vi.mock("@/lib/ses", () => ({ sendDatasetReadyEmail }));

const { POST } = await import("@/app/api/uploads/ready/route");

const post = (body: unknown) =>
  POST(
    new Request("http://localhost/api/uploads/ready", {
      method: "POST",
      body: JSON.stringify(body),
    }) as never,
  );

const row = (id: string, datasetType: string, extra = {}) => ({
  id,
  title: `set1 (${datasetType})`,
  datasetType,
  numCells: 0,
  numGenes: 0,
  platform: null,
  ...extra,
});

describe("a save uploaded as it is", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    requireUser.mockResolvedValue({
      error: null,
      session: { user: { id: "u1", email: "owner@example.com" } },
    });
  });

  it("sends one email for a folder's three parts, the cells first", async () => {
    // The database answers in its own order
    findMany.mockResolvedValue([
      row("ds_img", "image"),
      row("sm_mol", "single_molecule", { numCells: 1200000, numGenes: 300 }),
      row("ds_cells", "spatialdata", { numCells: 5000, numGenes: 300 }),
    ]);

    const res = await post({ datasetIds: ["ds_cells", "sm_mol", "ds_img"] });

    expect(await res.json()).toMatchObject({ emailed: true, counted: 3 });
    expect(sendDatasetReadyEmail).toHaveBeenCalledTimes(1);
    expect(sendDatasetReadyEmail).toHaveBeenCalledWith(
      expect.objectContaining({
        email: "owner@example.com",
        datasetId: "ds_cells",
        datasetType: "spatialdata",
        molecules: false,
        metadata: expect.objectContaining({ numCells: 5000, numGenes: 300 }),
        alsoSaved: [
          {
            datasetId: "sm_mol",
            datasetType: "single_molecule",
            datasetName: "set1 (single_molecule)",
          },
          {
            datasetId: "ds_img",
            datasetType: "image",
            datasetName: "set1 (image)",
          },
        ],
      }),
    );
  });

  it("counts nothing for an image on its own", async () => {
    findMany.mockResolvedValue([row("ds_img", "image")]);

    await post({ datasetIds: ["ds_img"] });

    const [args] = sendDatasetReadyEmail.mock.calls[0];

    expect(args.metadata.numCells).toBeUndefined();
    expect(args.metadata.numGenes).toBeUndefined();
    expect(args.alsoSaved).toEqual([]);
  });

  it("only counts the user's own finished datasets, and mails only them", async () => {
    findMany.mockResolvedValue([]);

    const res = await post({
      datasetIds: ["ds_someone_elses"],
      // An address in the request is never used
      email: "attacker@example.com",
    });

    expect(findMany).toHaveBeenCalledWith(
      expect.objectContaining({
        where: {
          id: { in: ["ds_someone_elses"] },
          ownerId: "u1",
          status: "COMPLETE",
        },
      }),
    );
    expect(await res.json()).toMatchObject({ emailed: false, counted: 0 });
    expect(sendDatasetReadyEmail).not.toHaveBeenCalled();
  });

  it("still succeeds when the email does not go out", async () => {
    findMany.mockResolvedValue([row("ds_img", "image")]);
    sendDatasetReadyEmail.mockRejectedValueOnce(new Error("SES down"));

    const res = await post({ datasetIds: ["ds_img"] });

    expect(res.status).toBe(200);
    expect(await res.json()).toMatchObject({
      success: true,
      emailed: false,
      counted: 1,
    });
  });

  it("needs dataset ids, and not too many", async () => {
    expect((await post({})).status).toBe(400);
    expect((await post({ datasetIds: [] })).status).toBe(400);
    expect(
      (await post({ datasetIds: Array.from({ length: 11 }, (_, i) => `d${i}`) }))
        .status,
    ).toBe(400);
    expect(findMany).not.toHaveBeenCalled();
  });

  it("needs a signed-in user", async () => {
    requireUser.mockResolvedValue({
      error: new Response("Unauthorized", { status: 401 }),
      session: null,
    });

    expect((await post({ datasetIds: ["ds_img"] })).status).toBe(401);
    expect(sendDatasetReadyEmail).not.toHaveBeenCalled();
  });
});
