import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { catalogRowWhere, isAdminRole } from "@/lib/explore/visibility";

// GET /api/explore/[id] — one catalog dataset, FULL row (including the genes
// array the list payload omits). Used by detail views that started from a
// slim list item, e.g. the in-viewer explore modal. Visibility is the
// listing's (lib/explore/visibility.ts): internal needs admin; community
// needs approval.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { id } = await params;
  const session = await auth();

  const dataset = await prisma.catalogDataset.findFirst({
    where: { AND: [{ id }, catalogRowWhere(isAdminRole(session?.user?.role))] },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });

  if (!dataset) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return NextResponse.json(dataset);
}
