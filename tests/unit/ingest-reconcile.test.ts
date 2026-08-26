import { describe, it, expect } from "vitest";

import { isStale, failureDetail } from "@/lib/ingest/reconcile";
import { classifyFailure } from "@/lib/ingest/error-classification";

/**
 * isStale is the gate that decides whether /status calls AWS Batch at all.
 * Too eager and every poll costs a DescribeJobs; too lax and a stranded dataset
 * is never recovered.
 */

const MINUTE = 60_000;

function ds(over: Partial<Parameters<typeof isStale>[0]> = {}) {
  return {
    id: "ds_x",
    status: "PROCESSING",
    batchJobId: "job-1",
    createdAt: new Date(Date.now() - 10 * MINUTE),
    processingProgress: { updatedAt: new Date().toISOString() },
    ...over,
  } as Parameters<typeof isStale>[0];
}

describe("isStale", () => {
  it("is false while progress is arriving", () => {
    expect(isStale(ds())).toBe(false);
  });

  it("is true once progress has gone quiet past the threshold", () => {
    expect(
      isStale(ds({ processingProgress: { updatedAt: new Date(Date.now() - 5 * MINUTE).toISOString() } })),
    ).toBe(true);
  });

  it("falls back to createdAt when no progress was ever reported", () => {
    // A job that dies before its first callback leaves processingProgress null;
    // without the createdAt fallback such a dataset could never be recovered.
    expect(isStale(ds({ processingProgress: null }))).toBe(true);
  });

  it("gives a freshly queued dataset time to start", () => {
    expect(
      isStale(ds({
        status: "QUEUED",
        processingProgress: null,
        createdAt: new Date(Date.now() - 5_000),
      })),
    ).toBe(false);
  });

  it("ignores datasets that already reached a terminal state", () => {
    for (const status of ["COMPLETE", "FAILED"]) {
      expect(isStale(ds({ status, processingProgress: null }))).toBe(false);
    }
  });

  it("ignores a dataset that was never submitted to Batch", () => {
    // No job id means nothing to ask about — e.g. SubmitJob itself failed.
    expect(isStale(ds({ batchJobId: null, processingProgress: null }))).toBe(false);
  });

  it("tolerates an unparseable progress timestamp", () => {
    expect(
      isStale(ds({ processingProgress: { updatedAt: "not-a-date" } })),
    ).toBe(true);
  });
});

describe("failureDetail", () => {
  // Real payload from an OOM-killed job (a6879df6, dev): the task-level
  // statusReason is the useless generic string; the container reason carries
  // the actual cause. Storing the generic one hid OOMs from the classifier,
  // the admin failures table, and the owner's failure email.
  it("prefers the container reason over the generic task statusReason", () => {
    expect(
      failureDetail({
        statusReason: "Essential container in task exited",
        attempts: [
          {
            container: {
              exitCode: 1,
              reason: "OutOfMemoryError: container killed due to memory usage",
            },
          },
        ],
      }),
    ).toBe("OutOfMemoryError: container killed due to memory usage");
  });

  it("keeps the exit code visible when no container reason exists", () => {
    expect(
      failureDetail({
        statusReason: "Essential container in task exited",
        attempts: [{ container: { exitCode: 137 } }],
      }),
    ).toBe(
      "container exited with exit code 137 (Essential container in task exited)",
    );
  });

  it("falls back through statusReason to a placeholder", () => {
    expect(failureDetail({ statusReason: "Task failed to start" })).toBe(
      "Task failed to start",
    );
    expect(failureDetail({})).toBe("no reason reported");
  });

  it("classifies the container OOM reason as an out-of-memory failure", () => {
    const detail = failureDetail({
      statusReason: "Essential container in task exited",
      attempts: [
        {
          container: {
            exitCode: 1,
            reason: "OutOfMemoryError: container killed due to memory usage",
          },
        },
      ],
    });
    const classified = classifyFailure(`Processing failed: ${detail}`);

    expect(classified.label).toBe("Out of memory");
    expect(classified.category).toBe("PLATFORM");
  });
});
