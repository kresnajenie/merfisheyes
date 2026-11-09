import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireUser } from "@/lib/admin-auth";

// GET /api/account/agents — the agents the signed-in person has allowed
// (OAuth grants not revoked), for the account page.
export async function GET() {
  const { error, session } = await requireUser();

  if (error) return error;
  const grants = await prisma.oAuthGrant.findMany({
    where: { userId: session.user.id, revokedAt: null },
    include: { client: { select: { name: true, kind: true } } },
    orderBy: { createdAt: "desc" },
  });
  const titles = new Map(
    (
      await prisma.project.findMany({
        where: {
          id: { in: grants.flatMap((g) => g.projectIds) },
          ownerId: session.user.id,
        },
        select: { id: true, title: true },
      })
    ).map((p) => [p.id, p.title]),
  );

  return NextResponse.json(
    grants.map((g) => ({
      id: g.id,
      name: g.client.name,
      verified: g.client.kind === "metadata_document",
      publisher:
        g.client.kind === "metadata_document" ? new URL(g.clientId).host : null,
      scopes: g.scopes,
      projects: g.allProjects
        ? "all"
        : g.projectIds.map((id) => titles.get(id) ?? "(deleted project)"),
      connectedAt: g.createdAt,
      lastUsedAt: g.lastUsedAt,
    })),
  );
}
