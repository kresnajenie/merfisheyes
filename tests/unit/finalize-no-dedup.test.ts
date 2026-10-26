import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * Duplicate checking is off (docs/plans/PREVIEW-SAVE.md, D4): an upload whose
 * raw files hash the same as another dataset's still completes, nothing is
 * looked up across datasets, and no other owner's dataset is named or
 * touched. The hash is still stored, as a record of what was uploaded.
 */

const FINGERPRINT = "a".repeat(64);

const tx = {
  uploadSession: { deleteMany: vi.fn(), create: vi.fn() },
  uploadFile: { createMany: vi.fn() },
  dataset: { update: vi.fn() },
};
const prisma = {
  dataset: {
    findUnique: vi.fn(),
    // Another user's dataset with the very same fingerprint: it must not be
    // asked for
    findFirst: vi.fn(async () => ({
      id: "ds_someone_else",
      title: "Someone else's data",
      datasetType: "spatialdata",
    })),
    update: vi.fn(),
  },
  $transaction: vi.fn(async (fn: (t: typeof tx) => Promise<void>) => fn(tx)),
};
const sendDatasetReadyEmail = vi.fn();
const writeOverlay = vi.fn();

vi.mock("@/lib/prisma", () => ({ prisma }));
vi.mock("@/lib/s3", () => ({
  listObjectKeys: vi.fn(async () => [
    { key: "datasets/ds_mine/mapping.json", size: 36 },
    { key: "datasets/ds_mine/ds_mine.zarr/zarr.json", size: 900 },
  ]),
  getObjectJson: vi.fn(),
  uploadBufferToS3: vi.fn(),
}));
vi.mock("@/lib/ses", () => ({ sendDatasetReadyEmail }));
vi.mock("@/lib/ingest/overlay", () => ({
  imageOverlayFor: vi.fn(async () => ({})),
  writeOverlay,
}));
vi.mock("@/lib/datasets/read-genes", () => ({ updateDatasetGenes: vi.fn() }));

const { finalizeCompletedDataset } = await import("@/lib/ingest/finalize");

describe("finalize with duplicate checking off", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    prisma.dataset.findUnique.mockResolvedValue({
      id: "ds_mine",
      title: "My upload",
      fingerprint: "raw_ds_mine",
      datasetType: "single_cell",
      formatVersion: null,
      manifestJson: null,
      numCells: null,
      numGenes: null,
      ownerId: "user_me",
      owner: { email: "me@example.com", role: "USER" },
      // Uploaded with molecules beside it: the old duplicate path would have
      // moved this overlay onto the other dataset
      processingParams: { pipeline: "prepare", linkedSmDatasetId: "sm_mine" },
    });
  });

  it("completes an upload whose files match another dataset's", async () => {
    const result = await finalizeCompletedDataset("ds_mine", {
      fingerprint: FINGERPRINT,
      stats: { numCells: 300, numGenes: 40 },
    });

    expect(result.registeredFiles).toBe(2);
    expect(tx.dataset.update).toHaveBeenCalledWith(
      expect.objectContaining({
        where: { id: "ds_mine" },
        data: expect.objectContaining({ status: "COMPLETE" }),
      }),
    );
    expect(
      prisma.dataset.update.mock.calls.some(
        ([arg]) =>
          (arg as { data: { status?: string } }).data.status === "FAILED",
      ),
    ).toBe(false);
  });

  it("never looks for other datasets, nor writes to one", async () => {
    await finalizeCompletedDataset("ds_mine", { fingerprint: FINGERPRINT });

    expect(prisma.dataset.findFirst).not.toHaveBeenCalled();
    expect(writeOverlay).not.toHaveBeenCalled();
    for (const [arg] of [
      ...prisma.dataset.update.mock.calls,
      ...tx.dataset.update.mock.calls,
    ]) {
      expect((arg as { where: { id: string } }).where.id).toBe("ds_mine");
    }
  });

  it("stores the fingerprint and sends the ready email, not a duplicate one", async () => {
    await finalizeCompletedDataset("ds_mine", { fingerprint: FINGERPRINT });

    expect(prisma.dataset.update).toHaveBeenCalledWith({
      where: { id: "ds_mine" },
      data: { fingerprint: FINGERPRINT },
    });
    expect(sendDatasetReadyEmail).toHaveBeenCalledTimes(1);
    expect(sendDatasetReadyEmail.mock.calls[0][0]).toMatchObject({
      email: "me@example.com",
      datasetId: "ds_mine",
    });
  });
});
