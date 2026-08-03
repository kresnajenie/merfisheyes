import { describe, it, expect } from "vitest";

import { isStale } from "@/lib/ingest/reconcile";

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
