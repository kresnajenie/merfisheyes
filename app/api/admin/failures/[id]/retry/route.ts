import { NextRequest, NextResponse } from "next/server";

import { requireAdmin } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { resubmitProcessing } from "@/lib/ingest/resubmit";
import { planOomEscalation } from "@/lib/ingest/escalation";
import { classifyFailure } from "@/lib/ingest/error-classification";
import {
  computeTierMemoryGb,
  isComputeTier,
  normalizeComputeTier,
} from "@/lib/batch";
import { sendDatasetRetryingEmail } from "@/lib/ses";

/**
 * POST /api/admin/failures/[id]/retry — re-run processing for a FAILED
 * dataset from its raw upload (admin only).
 *
 * Body: { computeTier?: "standard" | "large" | "xlarge" }. When omitted the
 * tier is the one the failed run used, stepped up one level if that failure
 * was an out-of-memory. The owner is emailed that processing restarted.
 */
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error } = await requireAdmin();

  if (error) return error;

  const { id } = await params;
  const body = await request.json().catch(() => ({}));

  const dataset = await prisma.dataset.findUnique({
    where: { id },
    select: {
      id: true,
      title: true,
      status: true,
      batchJobId: true,
      errorMessage: true,
      computeTier: true,
      processingAttempts: true,
      owner: { select: { email: true, computeTier: true } },
    },
  });

  if (!dataset) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }
  // FAILED rows, or a QUEUED row whose SubmitJob never happened.
  if (
    dataset.status !== "FAILED" &&
    !(dataset.status === "QUEUED" && !dataset.batchJobId)
  ) {
    return NextResponse.json(
      { error: `Dataset is ${dataset.status}, not retryable` },
      { status: 409 },
    );
  }

  const ranOn = normalizeComputeTier(
    dataset.computeTier ?? dataset.owner?.computeTier,
  );
  let tier = ranOn;

  if (isComputeTier(body?.computeTier)) {
    tier = body.computeTier;
  } else {
    const plan = planOomEscalation({
      label: classifyFailure(dataset.errorMessage ?? "").label,
      currentTier: ranOn,
      attempts: 0, // explicit admin retries aren't capped by the auto limit
    });

    if (plan) tier = plan.nextTier;
  }

  try {
    const result = await resubmitProcessing(id, tier);

    if (dataset.owner?.email) {
      try {
        await sendDatasetRetryingEmail({
          email: dataset.owner.email,
          datasetId: id,
          datasetName: dataset.title ?? undefined,
          reason: "admin",
          previousMemoryGb: computeTierMemoryGb(ranOn),
          nextMemoryGb: computeTierMemoryGb(result.computeTier),
        });
      } catch (e: any) {
        console.error(`Admin retry: email for ${id} failed:`, e?.message);
      }
    }

    return NextResponse.json({ ok: true, ...result });
  } catch (e: any) {
    return NextResponse.json(
      { error: e?.message || "Retry failed" },
      { status: 400 },
    );
  }
}
