import { DescribeJobsCommand } from "@aws-sdk/client-batch";

import { batchClient } from "@/lib/batch";
import { getObjectJson } from "@/lib/s3";
import { finalizeCompletedDataset } from "@/lib/ingest/finalize";
import { markDatasetFailed } from "@/lib/ingest/mark-failed";

/**
 * Recover a dataset whose worker callback never arrived.
 *
 * Status normally flows worker → /callback → DB. When that path breaks — a
 * misconfigured callback URL, deployment protection blocking the POST, the
 * container dying before it could report — the dataset sits at QUEUED or
 * PROCESSING forever and the UI polls indefinitely. AWS Batch still knows what
 * happened to the job, so we ask it.
 *
 * This only ever runs against a *stale* dataset (see STALE_AFTER_MS). A healthy
 * job posts progress every few seconds, so the common case costs no API calls
 * at all.
 */

// Long enough to cover Fargate cold start (~40s) plus queue wait, so a job that
// is merely slow to start is never mistaken for a lost one.
const STALE_AFTER_MS = 120_000;

/** Batch states that mean "still legitimately in flight". */
const PENDING_STATES = new Set([
  "SUBMITTED",
  "PENDING",
  "RUNNABLE",
  "STARTING",
  "RUNNING",
]);

interface ReconcilableDataset {
  id: string;
  status: string;
  batchJobId: string | null;
  createdAt: Date;
  processingProgress: unknown;
}

function lastSignalAt(dataset: ReconcilableDataset): number {
  const progress = dataset.processingProgress as { updatedAt?: string } | null;
  const updated = progress?.updatedAt ? Date.parse(progress.updatedAt) : NaN;

  return Number.isNaN(updated) ? dataset.createdAt.getTime() : updated;
}

export function isStale(dataset: ReconcilableDataset): boolean {
  if (dataset.status !== "QUEUED" && dataset.status !== "PROCESSING") {
    return false;
  }
  if (!dataset.batchJobId) return false;

  return Date.now() - lastSignalAt(dataset) > STALE_AFTER_MS;
}

export type ReconcileOutcome =
  | { changed: false; jobStatus?: string }
  | { changed: true; status: "COMPLETE" | "FAILED"; reason?: string };

interface FailedJobShape {
  statusReason?: string;
  attempts?: {
    statusReason?: string;
    container?: { reason?: string; exitCode?: number };
  }[];
}

/**
 * The most specific failure reason a Batch job record offers.
 *
 * Order matters: the container-level reason (e.g. "OutOfMemoryError:
 * container killed due to memory usage") is what actually happened, while the
 * task-level statusReason is usually the generic "Essential container in task
 * exited" — storing that hid every OOM from the fault classifier, the admin
 * failures table, and the owner's failure email. A bare exit code is kept as
 * a fallback signal (the classifier recognises "exit code 137" as OOM).
 */
export function failureDetail(job: FailedJobShape): string {
  const attempt = job.attempts?.[(job.attempts?.length ?? 0) - 1];
  const container = attempt?.container;

  if (container?.reason) return container.reason;
  if (container?.exitCode != null) {
    const generic = job.statusReason || attempt?.statusReason;

    return `container exited with exit code ${container.exitCode}${
      generic ? ` (${generic})` : ""
    }`;
  }

  return job.statusReason || attempt?.statusReason || "no reason reported";
}

/**
 * Ask Batch what actually happened, and write the answer through if the job
 * has finished. Returns whether the dataset row changed.
 *
 * Never throws: reconciliation is a safety net on a read path, so a Batch or S3
 * problem must not take down /status itself.
 */
export async function reconcileWithBatch(
  dataset: ReconcilableDataset,
): Promise<ReconcileOutcome> {
  if (!dataset.batchJobId) return { changed: false };

  try {
    const described = await batchClient.send(
      new DescribeJobsCommand({ jobs: [dataset.batchJobId] }),
    );
    const job = described.jobs?.[0];

    // Batch drops job records after ~24h. A dataset still QUEUED/PROCESSING
    // with no job record is never going to progress on its own.
    if (!job) {
      await markFailed(
        dataset.id,
        `Processing job ${dataset.batchJobId} is no longer known to AWS Batch (records expire after ~24h). The job's outcome could not be recovered; re-upload to retry.`,
      );

      return { changed: true, status: "FAILED", reason: "job record expired" };
    }

    if (PENDING_STATES.has(job.status ?? "")) {
      return { changed: false, jobStatus: job.status };
    }

    if (job.status === "FAILED") {
      const detail = failureDetail(job);

      await markFailed(dataset.id, `Processing failed: ${detail}`);

      return { changed: true, status: "FAILED", reason: detail };
    }

    if (job.status === "SUCCEEDED") {
      // The work is done and the output is in S3 — only the callback was lost.
      // Recover the statistics from the dataset's own manifest so the row isn't
      // completed with zeroed counts. Single-molecule writes a gzipped manifest
      // with molecule/gene counts; single-cell a plain one with cell/gene counts.
      const isSM = dataset.id.startsWith("sm_");
      const manifest = await getObjectJson<any>(
        `datasets/${dataset.id}/manifest.json${isSM ? ".gz" : ""}`,
      );

      if (!manifest) {
        await markFailed(
          dataset.id,
          `Processing job ${dataset.batchJobId} succeeded but produced no manifest. Re-upload to retry.`,
        );

        return { changed: true, status: "FAILED", reason: "no manifest" };
      }

      const s = manifest?.statistics ?? {};

      await finalizeCompletedDataset(dataset.id, {
        stats: isSM
          ? { numCells: s.total_molecules, numGenes: s.unique_genes }
          : { numCells: s.total_cells, numGenes: s.total_genes },
      });
      console.warn(
        `Ingest reconcile: ${dataset.id} completed from Batch state; its worker callback never arrived.`,
      );

      return { changed: true, status: "COMPLETE" };
    }

    return { changed: false, jobStatus: job.status };
  } catch (e: any) {
    console.error(`Ingest reconcile failed for ${dataset.id}:`, e?.message);

    return { changed: false };
  }
}

async function markFailed(datasetId: string, message: string): Promise<void> {
  // Routes through the shared helper so a reconcile-detected failure classifies
  // the fault and emails the owner, exactly like a worker-reported one.
  await markDatasetFailed(datasetId, message);
}
