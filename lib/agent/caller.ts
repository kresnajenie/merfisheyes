import { auth } from "@/lib/auth";
import { isAdminRole } from "@/lib/explore/visibility";
import { PUBLIC_SCOPE, type Scope } from "@/lib/oauth/config";
import { unauthorized } from "@/lib/agent/respond";

/**
 * Who is asking the agent API (docs/specs/MCP-PROTOCOL.md §1): an OAuth
 * bearer token, the site's own session cookie, or nobody.
 */

export type Caller =
  | { kind: "anonymous" }
  | {
      kind: "user";
      userId: string;
      role: string;
      isAdmin: boolean;
      via: "oauth" | "session";
      /** null = every scope (a session). */
      scopes: Set<string> | null;
      projects: "all" | string[];
      grantId: string | null;
      clientName: string | null;
    };

export function hasScope(caller: Caller, scope: Scope): boolean {
  if (scope === PUBLIC_SCOPE) return true;
  if (caller.kind === "anonymous") return false;

  return caller.scopes === null || caller.scopes.has(scope);
}

/** The caller, or the 401 to send back: a bearer token that is present but
 * not valid for this API is refused, never treated as anonymous. */
export async function resolveCaller(req: Request): Promise<Caller | Response> {
  const header = req.headers.get("authorization");

  if (header) {
    // No access tokens are issued yet (the authorization server is step 3 of
    // docs/plans/MCP.md); a token is refused, never ignored.
    return unauthorized(
      "invalid_token",
      "the access token is invalid, expired, revoked, or for another resource",
    );
  }
  const session = await auth();

  if (session?.user?.id) {
    return {
      kind: "user",
      userId: session.user.id,
      role: session.user.role ?? "USER",
      isAdmin: isAdminRole(session.user.role),
      via: "session",
      scopes: null,
      projects: "all",
      grantId: null,
      clientName: null,
    };
  }

  return { kind: "anonymous" };
}
