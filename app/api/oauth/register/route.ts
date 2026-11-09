import { ClientError, registerClient } from "@/lib/oauth/clients";
import {
  oauthError,
  oauthJson,
  preflight,
  rateLimited,
} from "@/lib/oauth/http";
import { LIMITS, checkRateLimits, clientIp } from "@/lib/rate-limit";

// POST /api/oauth/register — dynamic client registration (RFC 7591).
// Deprecated in MCP in favour of Client ID Metadata Documents; kept for
// clients that do not support those yet (docs/specs/MCP-PROTOCOL.md §6.2).
export async function POST(req: Request) {
  const limited = await checkRateLimits(LIMITS.register(clientIp(req)));

  if (!limited.ok) return rateLimited(limited.retryAfter);
  const body = await req.json().catch(() => null);

  try {
    return oauthJson(await registerClient(body), 201);
  } catch (e) {
    if (e instanceof ClientError) return oauthError(e.code, e.message);
    console.error("oauth register:", e);

    return oauthError("server_error", "registration failed; try again", 500);
  }
}

export const OPTIONS = preflight;
