import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";

// GET /api/explore/[id] — one catalog dataset, FULL row (including the genes
// array the list payload omits). Used by detail views that started from a
// slim list item, e.g. the in-viewer explore modal.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { id } = await params;

  const dataset = await prisma.catalogDataset.findUnique({
    where: { id },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });

  if (!dataset || !dataset.isPublished) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  // Same visibility rules as the listing: internal needs admin; community
  // needs approval.
  if (dataset.isInternal) {
    const session = await auth();
    const isAdmin =
      session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";

    if (!isAdmin) {
      return NextResponse.json({ error: "Not found" }, { status: 404 });
    }
  }
  if (dataset.isCommunity && dataset.reviewStatus !== "APPROVED") {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return NextResponse.json(dataset);
}
