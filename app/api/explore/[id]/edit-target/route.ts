import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";

// GET /api/explore/[id]/edit-target — where (if anywhere) the signed-in
// viewer can edit the source behind this catalog row. Explore detail pages
// render server-side without auth() to stay cacheable, so the owner's Edit
// button resolves its target through this per-session endpoint instead.
// Returns { editHref: string | null }.
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { id } = await params;
  const session = await auth();

  if (!session?.user) {
    return NextResponse.json({ editHref: null });
  }

  const catalog = await prisma.catalogDataset.findUnique({
    where: { id },
    select: { sourceProjectId: true, sourceDatasetId: true },
  });

  if (!catalog) {
    return NextResponse.json({ editHref: null });
  }

  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";

  // Same permission rules as the account pages themselves.
  if (catalog.sourceProjectId) {
    const project = await prisma.project.findUnique({
      where: { id: catalog.sourceProjectId },
      select: { ownerId: true },
    });

    if (project && (project.ownerId === session.user.id || isAdmin)) {
      return NextResponse.json({
        editHref: `/account/projects/${catalog.sourceProjectId}`,
      });
    }
  }

  if (catalog.sourceDatasetId) {
    const dataset = await prisma.dataset.findUnique({
      where: { id: catalog.sourceDatasetId },
      select: { ownerId: true, adminOwned: true },
    });

    if (
      dataset &&
      (dataset.ownerId === session.user.id || (isAdmin && dataset.adminOwned))
    ) {
      return NextResponse.json({
        editHref: `/account/datasets/${catalog.sourceDatasetId}`,
      });
    }
  }

  return NextResponse.json({ editHref: null });
}
