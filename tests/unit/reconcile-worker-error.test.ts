import { beforeEach, describe, expect, it, vi } from "vitest";

/**
 * A failed job whose callback never arrived: the dataset should carry the
 * worker's own reason (left in S3 by worker/entrypoint.py's fail()), not
 * Batch's "container exited with exit code 1".
 */

const send = vi.fn();
const getObjectJson = vi.fn();
const markDatasetFailed = vi.fn();

vi.mock("@/lib/batch", () => ({ batchClient: { send } }));
vi.mock("@/lib/s3", () => ({ getObjectJson }));
vi.mock("@/lib/ingest/finalize", () => ({ finalizeCompletedDataset: vi.fn() }));
vi.mock("@/lib/ingest/mark-failed", () => ({ markDatasetFailed }));

const { reconcileWithBatch } = await import("@/lib/ingest/reconcile");

const failedJob = {
  jobs: [
    {
      status: "FAILED",
      statusReason: "Essential container in task exited",
      attempts: [{ container: { exitCode: 1 } }],
    },
  ],
};
const dataset = {
  id: "sm_x",
  status: "PROCESSING",
  batchJobId: "job-1",
  createdAt: new Date(0),
  processingProgress: null,
} as Parameters<typeof reconcileWithBatch>[0];

describe("reconcileWithBatch on a failed job", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    send.mockResolvedValue(failedJob);
  });

  it("records the reason the worker left", async () => {
    getObjectJson.mockResolvedValue({
      error: "Could not find the molecule columns: …",
    });

    await reconcileWithBatch(dataset);
    expect(getObjectJson).toHaveBeenCalledWith(
      "datasets/sm_x/worker-error.json",
    );
    expect(markDatasetFailed.mock.calls[0][1]).toBe(
      "Processing failed: Could not find the molecule columns: …",
    );
  });

  it("falls back to what Batch knows when the worker left nothing", async () => {
    getObjectJson.mockResolvedValue(null);

    await reconcileWithBatch(dataset);
    expect(markDatasetFailed.mock.calls[0][1]).toBe(
      "Processing failed: container exited with exit code 1 (Essential container in task exited)",
    );
  });
});
