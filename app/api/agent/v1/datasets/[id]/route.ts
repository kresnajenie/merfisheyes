import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { notFound } from "@/lib/agent/respond";
import { datasetDetailOut } from "@/lib/agent/shape";
import { findDataset, ownDatasetsWhere } from "@/lib/agent/visibility";

// GET /api/agent/v1/datasets/{id} — meyes_get_dataset for an app
// dataset, when visible (docs/specs/MCP-PROTOCOL.md §2).
export const GET = agentRoute<{ id: string }>(
  ["catalog:read"],
  async ({ caller, params }) => {
    const dataset = await findDataset(caller, params.id);

    if (!dataset) return notFound("dataset");
    const owned =
      caller.kind === "user" &&
      (await prisma.dataset.count({
        where: { AND: [{ id: dataset.id }, ownDatasetsWhere(caller)] },
      })) > 0;

    return NextResponse.json(await datasetDetailOut(dataset, owned));
  },
  { expensive: true },
);
