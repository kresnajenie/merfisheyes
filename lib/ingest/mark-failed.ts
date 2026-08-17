import { prisma } from "@/lib/prisma";
import { sendDatasetFailedEmail } from "@/lib/ses";
import { classifyFailure } from "@/lib/ingest/error-classification";

/**
 * Mark a dataset FAILED and notify its owner.
 *
 * The single chokepoint every processing failure flows through (the worker
 * callback and the /status reconciliation both call it), so the DB write, the
 * auto fault-classification, and the failure email never drift apart — mirrors
 * how COMPLETE goes through finalizeCompletedDataset.
 *
 * Never throws on the email: a mail problem must not stop the row from being
 * marked FAILED.
 */
export async function markDatasetFailed(
  datasetId: string,
  message: string,
  opts: { progress?: Record<string, unknown> } = {},
): Promise<void> {
  const errorMessage = (message || "Processing failed").slice(0, 2000);
  const classification = classifyFailure(errorMessage);

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
