import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireAdmin } from "@/lib/admin-auth";

const includeEntries = { entries: { orderBy: { sortOrder: "asc" as const } } };

const REVIEW_STATUSES = ["PENDING", "APPROVED", "REJECTED"] as const;

type ReviewStatus = (typeof REVIEW_STATUSES)[number];

// GET /api/admin/community?status=PENDING — community submissions to review.
export async function GET(req: NextRequest) {
  const { error } = await requireAdmin();

  if (error) return error;

  const statusParam = req.nextUrl.searchParams.get("status") ?? "PENDING";
  const status: ReviewStatus = REVIEW_STATUSES.includes(
    statusParam as ReviewStatus,
  )
    ? (statusParam as ReviewStatus)
    : "PENDING";

  const items = await prisma.catalogDataset.findMany({
    where: { isCommunity: true, reviewStatus: status },
    include: includeEntries,
    orderBy: { updatedAt: "desc" },
  });

  // Attach submitter info (name/email) for each row.
  const submitterIds = Array.from(
    new Set(
      items
        .map((i) => i.submittedById)
        .filter((v): v is string => Boolean(v)),
    ),
  );
  const submitters = await prisma.user.findMany({
    where: { id: { in: submitterIds } },
    select: { id: true, name: true, email: true },
  });
  const submitterById = new Map(submitters.map((u) => [u.id, u]));

  return NextResponse.json({
    items: items.map((item) => ({
      ...item,
      submitter: item.submittedById
        ? (submitterById.get(item.submittedById) ?? null)
        : null,
    })),
  });
}
