import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * POST /api/uploads/failed: a browser upload that failed emails the signed-in
 * user the reason, and marks their unfinished dataset FAILED
 * (kresnajenie/merfisheyes#166).
 */

const requireUser = vi.fn();
const findUnique = vi.fn();
const sendDatasetFailedEmail = vi.fn();
const markDatasetFailed = vi.fn();

vi.mock("@/lib/admin-auth", () => ({ requireUser }));
vi.mock("@/lib/prisma", () => ({ prisma: { dataset: { findUnique } } }));
vi.mock("@/lib/ses", () => ({ sendDatasetFailedEmail }));
vi.mock("@/lib/ingest/mark-failed", () => ({ markDatasetFailed }));

const { POST } = await import("@/app/api/uploads/failed/route");

const post = (body: unknown) =>
  POST(
    new Request("http://localhost/api/uploads/failed", {
      method: "POST",
      body: JSON.stringify(body),
    }) as never,
  );

describe("a failed browser upload", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    requireUser.mockResolvedValue({
      error: null,
      session: { user: { id: "u1", email: "owner@example.com" } },
    });
  });

  it("marks the user's unfinished dataset FAILED, which emails them", async () => {
    findUnique.mockResolvedValue({ ownerId: "u1", status: "UPLOADING" });

    const res = await post({ datasetId: "sm_1", error: "Failed to upload genes/A.bin.gz to S3" });

    expect(await res.json()).toMatchObject({ marked: true });
    expect(markDatasetFailed).toHaveBeenCalledWith("sm_1", "Failed to upload genes/A.bin.gz to S3");
    expect(sendDatasetFailedEmail).not.toHaveBeenCalled();
  });

  it("emails the reason when no dataset was made", async () => {
    const res = await post({ title: "set1", error: "Storage is unavailable" });

    expect(await res.json()).toMatchObject({ marked: false });
    expect(sendDatasetFailedEmail).toHaveBeenCalledWith(
      expect.objectContaining({
        email: "owner@example.com",
        datasetName: "set1",
        errorMessage: "Storage is unavailable",
      }),
    );
  });

  it("never touches someone else's dataset, nor a finished one", async () => {
    findUnique.mockResolvedValueOnce({ ownerId: "u2", status: "UPLOADING" });
    await post({ datasetId: "ds_other", error: "x" });
    findUnique.mockResolvedValueOnce({ ownerId: "u1", status: "COMPLETE" });
    await post({ datasetId: "ds_done", error: "x" });

    expect(markDatasetFailed).not.toHaveBeenCalled();
    // Only ever to the signed-in user
    for (const [args] of sendDatasetFailedEmail.mock.calls) {
      expect(args.email).toBe("owner@example.com");
    }
  });

  it("needs a signed-in user", async () => {
    requireUser.mockResolvedValue({
      error: new Response("Unauthorized", { status: 401 }),
      session: null,
    });

    expect((await post({ error: "x" })).status).toBe(401);
    expect(sendDatasetFailedEmail).not.toHaveBeenCalled();
  });
});
