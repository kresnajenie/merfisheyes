import { prisma } from "@/lib/prisma";
import { sendDatasetFailedEmail, sendDatasetRetryingEmail } from "@/lib/ses";
import { classifyFailure } from "@/lib/ingest/error-classification";
import { planOomEscalation } from "@/lib/ingest/escalation";
import { resubmitProcessing } from "@/lib/ingest/resubmit";
import { computeTierMemoryGb, normalizeComputeTier } from "@/lib/batch";

/**
 * Mark a dataset FAILED and notify its owner — or, for an out-of-memory
 * failure with a bigger compute tier still available, retry it there
 * automatically and tell the owner that instead.
 *
 * The single chokepoint every processing failure flows through (the worker
 * callback and the /status reconciliation both call it), so the DB write, the
 * auto fault-classification, the OOM escalation, and the emails never drift
 * apart — mirrors how COMPLETE goes through finalizeCompletedDataset.
 *
 * Never throws on the email: a mail problem must not stop the row from being
 * marked FAILED. An escalation that can't be submitted (raw upload expired,
 * Batch error) falls through to the normal FAILED path.
 */
export async function markDatasetFailed(
  datasetId: string,
  message: string,
  opts: { progress?: Record<string, unknown> } = {},
): Promise<void> {
  const errorMessage = (message || "Processing failed").slice(0, 2000);
  const classification = classifyFailure(errorMessage);

  const current = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: {
      title: true,
      computeTier: true,
      processingAttempts: true,
      owner: { select: { email: true, computeTier: true } },
    },
  });

  // ── Automatic OOM escalation ───────────────────────────────────────
  // Legacy rows (submitted before computeTier was recorded) ran on the
  // owner's tier, so fall back to that.
  const ranOn = normalizeComputeTier(
    current?.computeTier ?? current?.owner?.computeTier,
  );
  const plan = current
    ? planOomEscalation({
        label: classification.label,
        currentTier: ranOn,
        attempts: Math.max(1, current.processingAttempts ?? 0),
      })
    : null;

  if (plan) {
    try {
      const result = await resubmitProcessing(datasetId, plan.nextTier);

      console.log(
        `Ingest markFailed: ${datasetId} ran out of memory on ${ranOn}; ` +
          `retrying on ${result.computeTier} (attempt ${result.attempt}, job ${result.batchJobId})`,
      );

      if (current?.owner?.email) {
        try {
          await sendDatasetRetryingEmail({
            email: current.owner.email,
            datasetId,
            datasetName: current.title ?? undefined,
            reason: "out_of_memory",
            previousMemoryGb: computeTierMemoryGb(ranOn),
            nextMemoryGb: computeTierMemoryGb(result.computeTier),
          });
        } catch (e: any) {
          console.error(
            `Ingest markFailed: retry email for ${datasetId} failed:`,
            e?.message,
          );
        }
      }

      return;
    } catch (e: any) {
      // Can't escalate (raw expired, SubmitJob failed…) — record the failure
      // normally, with the reason appended so it isn't lost.
      console.error(
        `Ingest markFailed: OOM escalation for ${datasetId} failed:`,
        e?.message,
      );
    }
  }

  // ── Plain failure ──────────────────────────────────────────────────
  const dataset = await prisma.dataset.update({
    where: { id: datasetId },
    data: {
      status: "FAILED",
      errorMessage,
      completedAt: new Date(),
      faultCategory: classification.category,
      processingProgress: {
        ...(opts.progress ?? {}),
        failedAt: new Date().toISOString(),
      },
    },
    select: {
      id: true,
      title: true,
      owner: { select: { email: true } },
    },
  });

  if (dataset.owner?.email) {
    try {
      await sendDatasetFailedEmail({
        email: dataset.owner.email,
        datasetId: dataset.id,
        datasetName: dataset.title ?? undefined,
        errorMessage,
        userHint: classification.userHint,
      });
    } catch (e: any) {
      console.error(
        `Ingest markFailed: failure email for ${datasetId} failed:`,
        e?.message,
      );
    }
  }
}
