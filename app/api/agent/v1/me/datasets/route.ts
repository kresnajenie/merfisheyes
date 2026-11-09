import type { Prisma } from "@prisma/client";

import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { page, paging } from "@/lib/agent/respond";
import { datasetCardOut } from "@/lib/agent/shape";
import { ownDatasetsWhere } from "@/lib/agent/visibility";

// GET /api/agent/v1/me/datasets?q=&status= — meyes_list_my_datasets.
// The caller's datasets in reach, newest first; an upload still sending
// bytes is not listed (it has nothing to show yet).
export const GET = agentRoute(["projects:read"], async ({ url, caller }) => {
  if (caller.kind !== "user")
    throw new Error("unreachable: scope check admits users only");
  const p = paging(url);
  const q = url.searchParams.get("q")?.trim();
  const status = url.searchParams.get("status")?.toUpperCase();
  const where: Prisma.DatasetWhereInput = {
    AND: [
      ownDatasetsWhere(caller),
      { status: { not: "UPLOADING" } },
      ...(q ? [{ title: { contains: q, mode: "insensitive" as const } }] : []),
      ...(status &&
      ["QUEUED", "PROCESSING", "COMPLETE", "FAILED"].includes(status)
        ? [
            {
              status: status as "QUEUED" | "PROCESSING" | "COMPLETE" | "FAILED",
            },
          ]
        : []),
    ],
  };
  const [rows, total] = await Promise.all([
    prisma.dataset.findMany({
      where,
      orderBy: { createdAt: "desc" },
      skip: p.offset,
      take: p.limit,
    }),
    prisma.dataset.count({ where }),
  ]);

  return NextResponse.json(page(rows.map(datasetCardOut), total, p));
});
