import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { requireUser } from "@/lib/admin-auth";
import { revokeGrant } from "@/lib/oauth/flow";

// DELETE /api/account/agents/{grant id} — disconnect an agent: the grant and
// every token issued under it stop working at once.
export async function DELETE(
  _req: Request,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;
  const { id } = await params;
  const grant = await prisma.oAuthGrant.findFirst({
    where: { id, userId: session.user.id },
  });

  if (!grant) return NextResponse.json({ error: "Not found" }, { status: 404 });
  await revokeGrant(grant.id);

  return NextResponse.json({ ok: true });
}
