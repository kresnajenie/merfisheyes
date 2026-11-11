import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { agentRoute } from "@/lib/agent/route";
import { ALL_SCOPES, PUBLIC_SCOPE } from "@/lib/oauth/config";

// GET /api/agent/v1/me — meyes_whoami: the person, the scopes and the
// projects this token reaches.
export const GET = agentRoute([], async ({ caller }) => {
  if (caller.kind === "anonymous") {
    return NextResponse.json({
      signedIn: false,
      scopes: [PUBLIC_SCOPE],
      projects: [],
    });
  }
  const user = await prisma.user.findUnique({
    where: { id: caller.userId },
    select: { name: true, email: true },
  });

  return NextResponse.json({
    signedIn: true,
    name: user?.name ?? null,
    email: user?.email ?? null,
    via: caller.via,
    client: caller.clientName,
    scopes: caller.scopes ? [...caller.scopes] : ALL_SCOPES,
    projects: caller.projects,
  });
});
