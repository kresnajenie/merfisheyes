import { prisma } from "@/lib/prisma";
import { listObjectKeys } from "@/lib/s3";
import {
  isBatchConfigured,
  normalizeComputeTier,
  submitIngestJob,
  type ComputeTier,
} from "@/lib/batch";

/**
 * Submit a fresh processing job for a dataset whose raw upload is still in S3.
 *
 * Used by the admin "Retry" action and by automatic OOM escalation. The raw
 * objects under raw/{id}/ are the input the worker re-downloads, and they
 * expire 7 days after upload (S3 lifecycle) — after that a retry is impossible
 * and the owner has to re-upload.
 *
 * Resets the row to QUEUED (clearing the previous failure) and records the
 * tier + attempt count so a later escalation knows where it stands. Throws
 * with a human-readable message when the retry can't happen.
 */
export interface ResubmitResult {
  batchJobId: string;
  computeTier: ComputeTier;
  attempt: number;
}

export async function resubmitProcessing(
  datasetId: string,
  computeTier: ComputeTier,
): Promise<ResubmitResult> {
  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: {
      id: true,
      datasetType: true,
      processingParams: true,
      processingAttempts: true,
      ingestSource: true,
    },
  });

  if (!dataset) throw new Error("Dataset not found");
  if (dataset.ingestSource !== "server") {
    throw new Error("Only server-processed uploads can be retried");
  }
  if (!isBatchConfigured()) {
    throw new Error("Processing is not configured on this server");
  }

  const rawKeys = await listObjectKeys(`raw/${datasetId}/`);

  if (rawKeys.length === 0) {
    throw new Error(
      "The raw upload is no longer available (raw files expire 7 days after upload) — the owner needs to re-upload",
    );
  }

  const tier = normalizeComputeTier(computeTier);
  const attempt = (dataset.processingAttempts ?? 0) + 1;

  // QUEUED first so the state is right even if SubmitJob throws; the stale
  // reconciliation won't touch a row without a batchJobId.
  await prisma.dataset.update({
    where: { id: datasetId },
    data: {
      status: "QUEUED",
      errorMessage: null,
      faultCategory: null,
      completedAt: null,
      batchJobId: null,
      computeTier: tier,
      processingAttempts: attempt,
      processingProgress: {
        stage: "Queued for retry",
        attempt,
        updatedAt: new Date().toISOString(),
      },
    },
  });

  const batchJobId = await submitIngestJob({
    datasetId,
    kind:
      dataset.datasetType === "single_molecule"
        ? "single_molecule"
        : "single_cell",
    processingParams: dataset.processingParams,
    computeTier: tier,
  });

  await prisma.dataset.update({
    where: { id: datasetId },
    data: { batchJobId },
  });

  return { batchJobId, computeTier: tier, attempt };
}
