
import { resolveCaller } from "@/lib/agent/caller";
import { agentGet } from "@/lib/agent/in-process";
import { unauthorized } from "@/lib/agent/respond";
import { mcpHandler } from "@/lib/mcp/server";
import { issuer, mcpResource } from "@/lib/oauth/config";

/**
 * `/mcp` — MERFISHEYES as a remote MCP server (MCP.md pillar E): an agent
 * adds this address as a connector and installs nothing. Streamable HTTP,
 * stateless (a fresh server per request, so it runs on serverless
 * functions). Signing in is required (decided 2026-11-16): a request without
 * a token is answered 401 pointing at this resource's metadata, which is how
 * a connector discovers our OAuth server and signs the person in. A token
 * must be for this resource; the site's session cookie is never accepted
 * here. The tools are the npm package's own (one implementation, D23),
 * reaching the agent API in process as the verified caller.
 */
async function handle(req: Request): Promise<Response> {
  if (!req.headers.get("authorization"))
    return unauthorized(
      null,
      "sign in to MERFISHEYES to use this connector",
      undefined,
      mcpResource(),
    );
  const caller = await resolveCaller(req, mcpResource());

  if (caller instanceof Response) return caller;
  const site = issuer();

  return mcpHandler({
    site,
    get: (path, params, signal) => agentGet(caller, path, params, signal),
  }).fetch(req);
}

export const GET = handle;
export const POST = handle;
export const DELETE = handle;
export const dynamic = "force-dynamic";
