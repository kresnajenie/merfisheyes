import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { notFound } from "@/lib/agent/respond";
import { datasetCardOut } from "@/lib/agent/shape";
import { projectsWhere } from "@/lib/agent/visibility";

// GET /api/agent/v1/projects/{id} — meyes_get_project
export const GET = agentRoute<{ id: string }>(
  ["projects:read"],
  async ({ caller, params }) => {
    if (caller.kind !== "user")
      throw new Error("unreachable: scope check admits users only");
    const project = await prisma.project.findFirst({
      where: { AND: [{ id: params.id }, projectsWhere(caller)] },
      include: {
        datasets: { include: { dataset: true }, orderBy: { sortOrder: "asc" } },
      },
    });

    if (!project) return notFound("project");

    return NextResponse.json({
      id: project.id,
      title: project.title,
      description: project.description,
      species: project.species,
      tissue: project.tissue,
      platform: project.platform,
      tags: project.tags,
      metadata: project.metadata,
      createdAt: project.createdAt,
      updatedAt: project.updatedAt,
      datasets: project.datasets.map((m) => datasetCardOut(m.dataset)),
    });
  },
);
