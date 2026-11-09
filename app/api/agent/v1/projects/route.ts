import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { page, paging } from "@/lib/agent/respond";
import { projectsWhere } from "@/lib/agent/visibility";

// GET /api/agent/v1/projects?q= — meyes_list_projects: the caller's
// projects in reach, each with its dataset ids.
export const GET = agentRoute(["projects:read"], async ({ url, caller }) => {
  if (caller.kind !== "user")
    throw new Error("unreachable: scope check admits users only");
  const p = paging(url);
  const q = url.searchParams.get("q")?.trim();
  const where = {
    AND: [
      projectsWhere(caller),
      ...(q ? [{ title: { contains: q, mode: "insensitive" as const } }] : []),
    ],
  };
  const [rows, total] = await Promise.all([
    prisma.project.findMany({
      where,
      orderBy: { updatedAt: "desc" },
      skip: p.offset,
      take: p.limit,
      include: {
        datasets: {
          select: { datasetId: true },
          orderBy: { sortOrder: "asc" },
        },
      },
    }),
    prisma.project.count({ where }),
  ]);

  return NextResponse.json(
    page(
      rows.map((r) => ({
        id: r.id,
        title: r.title,
        species: r.species,
        tissue: r.tissue,
        platform: r.platform,
        datasetIds: r.datasets.map((d) => d.datasetId),
        updatedAt: r.updatedAt,
      })),
      total,
      p,
    ),
  );
});
