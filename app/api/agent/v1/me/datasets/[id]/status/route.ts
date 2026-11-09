import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { notFound } from "@/lib/agent/respond";
import { datasetCardOut } from "@/lib/agent/shape";
import { ownDatasetsWhere } from "@/lib/agent/visibility";
import { isStale, reconcileWithBatch } from "@/lib/ingest/reconcile";

// GET /api/agent/v1/me/datasets/{id}/status — meyes_get_processing_status.
// Reconciles a stale row with AWS Batch as /api/ingest/[id]/status does.
export const GET = agentRoute<{ id: string }>(
  ["projects:read"],
  async ({ caller, params }) => {
    if (caller.kind !== "user")
      throw new Error("unreachable: scope check admits users only");
    const where = { AND: [{ id: params.id }, ownDatasetsWhere(caller)] };
    let dataset = await prisma.dataset.findFirst({ where });

    if (!dataset) return notFound("dataset of yours");
    if (isStale(dataset) && (await reconcileWithBatch(dataset)).changed) {
      dataset = (await prisma.dataset.findFirst({ where })) ?? dataset;
    }

    return NextResponse.json({
      ...datasetCardOut(dataset),
      progress: dataset.processingProgress,
      errorMessage: dataset.errorMessage,
      faultCategory: dataset.faultCategory,
    });
  },
);
