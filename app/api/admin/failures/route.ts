import type { Prisma } from "@prisma/client";

import { NextRequest, NextResponse } from "next/server";

import { requireAdmin } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import {
  classifyFailure,
  OUT_OF_MEMORY_LABEL,
} from "@/lib/ingest/error-classification";
import { nextComputeTier, normalizeComputeTier } from "@/lib/batch";

/**
 * GET /api/admin/failures — paginated, sortable list of FAILED datasets for
 * the admin failures page. Returns the FULL error message (the dashboard's
 * inline table truncated it, which made errors unreadable).
 *
 * Query params: page (1-based), limit, sort (one of SORTS), dir (asc|desc).
 * Default sort: most recent upload first.
 */
const SORT_FIELDS = ["created", "failedAt", "title", "type", "owner", "fault"];

function orderBy(
  sort: string,
  dir: "asc" | "desc",
): Prisma.DatasetOrderByWithRelationInput[] {
  switch (sort) {
    case "title":
      return [{ title: dir }, { createdAt: "desc" }];
    case "type":
      return [{ datasetType: dir }, { createdAt: "desc" }];
    case "owner":
      return [{ owner: { email: dir } }, { createdAt: "desc" }];
    case "fault":
      // Nulls (untagged) group together at one end; recency breaks ties.
      return [{ faultCategory: dir }, { createdAt: "desc" }];
    case "failedAt":
      return [{ completedAt: dir }, { createdAt: dir }];
    case "created":
    default:
      return [{ createdAt: dir }];
  }
}

export async function GET(request: NextRequest) {
  const { error } = await requireAdmin();

  if (error) return error;

  const params = request.nextUrl.searchParams;
  const page = Math.max(1, Number(params.get("page") ?? "1") || 1);
  const limit = Math.min(
    100,
    Math.max(1, Number(params.get("limit") ?? "20") || 20),
  );
  const sortRaw = params.get("sort") ?? "created";
  const sort = SORT_FIELDS.includes(sortRaw) ? sortRaw : "created";
  const dir = params.get("dir") === "asc" ? "asc" : "desc";

  const where = { status: "FAILED" as const };

  const [rows, total] = await Promise.all([
    prisma.dataset.findMany({
      where,
      orderBy: orderBy(sort, dir),
      skip: (page - 1) * limit,
      take: limit,
      select: {
        id: true,
        title: true,
        datasetType: true,
        errorMessage: true,
        faultCategory: true,
        createdAt: true,
        completedAt: true,
        computeTier: true,
        processingAttempts: true,
        owner: { select: { name: true, email: true, computeTier: true } },
      },
    }),
    prisma.dataset.count({ where }),
  ]);

  return NextResponse.json({
    total,
    page,
    limit,
    failures: rows.map((f) => {
      const label = classifyFailure(f.errorMessage).label;
      // Tier the failed run used (legacy rows: the owner's tier), and the
      // tier a retry would default to (one step up after an OOM).
      const ranOn = normalizeComputeTier(f.computeTier ?? f.owner?.computeTier);
      const suggested =
        label === OUT_OF_MEMORY_LABEL
          ? (nextComputeTier(ranOn) ?? ranOn)
          : ranOn;

      return {
        ...f,
        autoLabel: label,
        ranOnTier: ranOn,
        suggestedTier: suggested,
      };
    }),
  });
}
