import { NextResponse } from "next/server";

import { requireAdmin } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { classifyFailure } from "@/lib/ingest/error-classification";

/**
 * GET /api/admin/processing — upload/processing stats for the admin dashboard.
 *
 * Aggregates the real upload `Dataset` rows (NOT the Explore catalog): status
 * totals + success rate, the recent failures to act on (with owner + fault),
 * per-user activity, and recent-activity windows.
 */
export async function GET() {
  const { error } = await requireAdmin();

  if (error) return error;

  const now = Date.now();
  const days = (n: number) => new Date(now - n * 24 * 3600 * 1000);
  const d7 = days(7);
  const d30 = days(30);

  const [statusGroups, ownerStatusGroups, recentFailedRaw, created7, created30] =
    await Promise.all([
      // Status totals across all upload datasets.
      prisma.dataset.groupBy({ by: ["status"], _count: { _all: true } }),
      // Per-(owner, status) counts → folded into per-user activity below.
      prisma.dataset.groupBy({
        by: ["ownerId", "status"],
        _count: { _all: true },
      }),
      // The failures to act on (most recent first).
      prisma.dataset.findMany({
        where: { status: "FAILED" },
        orderBy: [{ completedAt: "desc" }, { createdAt: "desc" }],
        take: 100,
        select: {
          id: true,
          title: true,
          datasetType: true,
          errorMessage: true,
          faultCategory: true,
          createdAt: true,
          completedAt: true,
          owner: { select: { name: true, email: true } },
        },
      }),
      prisma.dataset.count({ where: { createdAt: { gte: d7 } } }),
      prisma.dataset.count({ where: { createdAt: { gte: d30 } } }),
    ]);

  // Status totals as a plain map.
  const statusCounts: Record<string, number> = {
    UPLOADING: 0,
    QUEUED: 0,
    PROCESSING: 0,
    COMPLETE: 0,
    FAILED: 0,
  };

  for (const g of statusGroups) statusCounts[g.status] = g._count._all;

  const total = Object.values(statusCounts).reduce((a, b) => a + b, 0);
  const complete = statusCounts.COMPLETE;
  const failed = statusCounts.FAILED;
  const finished = complete + failed;
  const successRate = finished > 0 ? complete / finished : null;

  // Fold per-(owner, status) into per-user rows.
  const byOwner = new Map<
    string,
    { total: number; complete: number; failed: number; inProgress: number }
  >();

  for (const g of ownerStatusGroups) {
    if (!g.ownerId) continue; // skip unowned/admin-collective rows
    const row = byOwner.get(g.ownerId) ?? {
      total: 0,
      complete: 0,
      failed: 0,
      inProgress: 0,
    };
    const c = g._count._all;

    row.total += c;
    if (g.status === "COMPLETE") row.complete += c;
    else if (g.status === "FAILED") row.failed += c;
    else if (g.status === "PROCESSING" || g.status === "QUEUED")
      row.inProgress += c;
    byOwner.set(g.ownerId, row);
  }

  // Attach user name/email to each owner row.
  const ownerIds = [...byOwner.keys()];
  const users = ownerIds.length
    ? await prisma.user.findMany({
        where: { id: { in: ownerIds } },
        select: { id: true, name: true, email: true },
      })
    : [];
  const userById = new Map(users.map((u) => [u.id, u]));

  const perUser = ownerIds
    .map((id) => {
      const u = userById.get(id);
      const c = byOwner.get(id)!;

      return {
        userId: id,
        name: u?.name ?? null,
        email: u?.email ?? null,
        ...c,
      };
    })
    .sort((a, b) => b.total - a.total)
    .slice(0, 50);

  // Distinct active uploaders in each window.
  const [active7, active30] = await Promise.all([
    prisma.dataset.findMany({
      where: { createdAt: { gte: d7 }, ownerId: { not: null } },
      distinct: ["ownerId"],
      select: { ownerId: true },
    }),
    prisma.dataset.findMany({
      where: { createdAt: { gte: d30 }, ownerId: { not: null } },
      distinct: ["ownerId"],
      select: { ownerId: true },
    }),
  ]);

  // Attach the auto-classified label to each failure (the stored faultCategory
  // is the admin-facing tag; the classifier label is the "why" hint).
  const failures = recentFailedRaw.map((f) => ({
    ...f,
    autoLabel: classifyFailure(f.errorMessage).label,
  }));

  return NextResponse.json({
    statusCounts,
    total,
    successRate,
    failedCount: failed,
    uploaders: perUser.length,
    recent: {
      created7,
      created30,
      active7: active7.length,
      active30: active30.length,
    },
    failures,
    perUser,
  });
}
