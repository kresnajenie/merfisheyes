import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * The routes an upload calls as it goes, and the "ready" emails, answer only
 * the dataset's signed-in owner (requireUploadOwner in lib/admin-auth.ts).
 * They used to trust whoever knew a dataset id and an upload id — and the
 * email routes sent our mail to any address in the request body.
 */

const session = vi.fn<() => Promise<unknown>>();
const db = {
  dataset: { findUnique: vi.fn(), findFirst: vi.fn(), update: vi.fn() },
  uploadSession: { findUnique: vi.fn(), update: vi.fn() },
  uploadFile: { updateMany: vi.fn(), count: vi.fn() },
  $transaction: vi.fn(),
};
const sendDatasetReadyEmail = vi.fn();
const sesSend = vi.fn();
const createPresignedPost = vi.fn();

vi.mock("@/lib/auth", () => ({ auth: () => session() }));
vi.mock("@/lib/prisma", () => ({ prisma: db }));
vi.mock("@/lib/ses", () => ({ sendDatasetReadyEmail }));
vi.mock("@aws-sdk/client-ses", () => ({
  SESClient: class {
    send = sesSend;
  },
  SendEmailCommand: class {
    constructor(readonly input: unknown) {}
  },
}));
vi.mock("@aws-sdk/s3-presigned-post", () => ({ createPresignedPost }));

const OWNER = { user: { id: "owner", email: "owner@example.com", role: "USER" } };
const STRANGER = { user: { id: "stranger", email: "x@example.com", role: "USER" } };
const ADMIN = { user: { id: "admin", email: "admin@example.com", role: "ADMIN" } };

const post = (body: unknown) =>
  new Request("http://localhost/x", { method: "POST", body: JSON.stringify(body) }) as never;
const params = <T>(p: T) => ({ params: Promise.resolve(p) });

// Each upload route, called for dataset ds_1 with upload up_1.
const ROUTES: [string, () => Promise<Response>][] = [
  ["datasets/[id]/complete", async () =>
    (await import("@/app/api/datasets/[datasetId]/complete/route")).POST(post({ uploadId: "up_1" }), params({ datasetId: "ds_1" }))],
  ["datasets/[id]/files/complete", async () =>
    (await import("@/app/api/datasets/[datasetId]/files/complete/route")).POST(post({ uploadId: "up_1", fileKeys: ["a"] }), params({ datasetId: "ds_1" }))],
  ["datasets/[id]/files/[key]/complete", async () =>
    (await import("@/app/api/datasets/[datasetId]/files/[fileKey]/complete/route")).POST(post({ uploadId: "up_1" }), params({ datasetId: "ds_1", fileKey: "a" }))],
  ["datasets/[id]/upload-progress", async () =>
    (await import("@/app/api/datasets/[datasetId]/upload-progress/route")).POST(post({ uploadId: "up_1", delta: 1 }), params({ datasetId: "ds_1" }))],
  ["single-molecule/[id]/complete", async () =>
    (await import("@/app/api/single-molecule/[id]/complete/route")).POST(post({ uploadId: "up_1" }), params({ id: "ds_1" }))],
  ["single-molecule/[id]/files/[key]/complete", async () =>
    (await import("@/app/api/single-molecule/[id]/files/[key]/complete/route")).POST(post({ uploadId: "up_1" }), params({ id: "ds_1", key: "a" }))],
];

const writes = () =>
  db.uploadFile.updateMany.mock.calls.length +
  db.uploadSession.update.mock.calls.length +
  db.dataset.update.mock.calls.length +
  db.$transaction.mock.calls.length;

beforeEach(() => {
  vi.clearAllMocks();
  session.mockResolvedValue(null);
  db.dataset.findUnique.mockResolvedValue({ ownerId: "owner", adminOwned: false, datasetType: "spatialdata" });
  db.uploadSession.findUnique.mockResolvedValue({ datasetId: "ds_1" });
  db.uploadFile.updateMany.mockResolvedValue({ count: 1 });
});

describe("upload routes", () => {
  it.each(ROUTES)("%s: signed out → 401, nothing written", async (_, call) => {
    expect((await call()).status).toBe(401);
    expect(writes()).toBe(0);
  });

  it.each(ROUTES)("%s: someone else's dataset → 403, nothing written", async (_, call) => {
    session.mockResolvedValue(STRANGER);
    expect((await call()).status).toBe(403);
    expect(writes()).toBe(0);
  });

  it.each(ROUTES)("%s: an upload of another dataset → 404, nothing written", async (_, call) => {
    session.mockResolvedValue(OWNER);
    db.uploadSession.findUnique.mockResolvedValue({ datasetId: "ds_other" });
    expect((await call()).status).toBe(404);
    expect(writes()).toBe(0);
  });

  it("lets the owner through, and an admin on an admin-owned dataset", async () => {
    const [, call] = ROUTES[2];

    session.mockResolvedValue(OWNER);
    expect((await call()).status).toBe(200);
    expect(db.uploadFile.updateMany).toHaveBeenCalledOnce();

    session.mockResolvedValue(ADMIN);
    db.dataset.findUnique.mockResolvedValue({ ownerId: null, adminOwned: true });
    expect((await call()).status).toBe(200);

    db.dataset.findUnique.mockResolvedValue({ ownerId: "owner", adminOwned: false });
    expect((await call()).status).toBe(403); // an admin is not every owner
  });
});

describe("the ready emails", () => {
  const send = async () =>
    (await import("@/app/api/send-email/route")).POST(
      post({ email: "victim@example.com", datasetId: "ds_1", datasetName: "x" }),
    );
  const sendMolecules = async () =>
    (await import("@/app/api/send-email-single-molecule/route")).POST(
      post({ email: "victim@example.com", datasetId: "ds_1" }),
    );

  it("refuse a stranger and the signed-out", async () => {
    expect((await send()).status).toBe(401);
    expect((await sendMolecules()).status).toBe(401);
    session.mockResolvedValue(STRANGER);
    expect((await send()).status).toBe(403);
    expect((await sendMolecules()).status).toBe(403);
    expect(sendDatasetReadyEmail).not.toHaveBeenCalled();
    expect(sesSend).not.toHaveBeenCalled();
  });

  it("go to the owner's own address, whatever the body says", async () => {
    session.mockResolvedValue(OWNER);
    expect((await send()).status).toBe(200);
    expect(sendDatasetReadyEmail).toHaveBeenCalledWith(
      expect.objectContaining({ email: "owner@example.com", datasetId: "ds_1" }),
    );
    db.dataset.findUnique.mockResolvedValue({ ownerId: "owner", adminOwned: false, title: "t", numCells: 1, numGenes: 1, datasetType: "single_molecule" });
    expect((await sendMolecules()).status).toBe(200);
    expect(sesSend.mock.calls[0][0].input.Destination.ToAddresses).toEqual(["owner@example.com"]);
  });
});

describe("check-duplicate", () => {
  const check = async (route: "datasets" | "single-molecule") =>
    (route === "datasets"
      ? await import("@/app/api/datasets/check-duplicate/[fingerprint]/route")
      : await import("@/app/api/single-molecule/check-duplicate/[fingerprint]/route")
    ).GET(new Request("http://localhost/x") as never, params({ fingerprint: "f".repeat(16) }));

  it("needs sign-in and looks among the caller's own datasets only", async () => {
    db.dataset.findFirst.mockResolvedValue(null);
    for (const r of ["datasets", "single-molecule"] as const) {
      session.mockResolvedValue(null);
      expect((await check(r)).status).toBe(401);
      session.mockResolvedValue(OWNER);
      expect((await check(r)).status).toBe(200);
      expect(db.dataset.findFirst.mock.lastCall![0].where).toMatchObject({ ownerId: "owner" });
    }
  });
});

describe("initiate-zarr", () => {
  it("signs no upload policy for the signed-out", async () => {
    const { POST } = await import("@/app/api/datasets/initiate-zarr/route");
    const res = await POST(post({ fingerprint: "f", metadata: { numCells: 1, numGenes: 1 }, totalFiles: 1 }));

    expect(res.status).toBe(401);
    expect(createPresignedPost).not.toHaveBeenCalled();
  });
});
