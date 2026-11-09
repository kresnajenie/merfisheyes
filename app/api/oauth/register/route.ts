import { ClientError, registerClient } from "@/lib/oauth/clients";
import { oauthError, oauthJson, preflight } from "@/lib/oauth/http";

// POST /api/oauth/register — dynamic client registration (RFC 7591).
// Deprecated in MCP in favour of Client ID Metadata Documents; kept for
// clients that do not support those yet (docs/specs/MCP-PROTOCOL.md §6.2).
export async function POST(req: Request) {
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
