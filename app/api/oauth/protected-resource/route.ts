import { DEFAULT_SCOPES, SCOPES, issuer, resources } from "@/lib/oauth/config";
import { oauthError, oauthJson, preflight } from "@/lib/oauth/http";

// GET /.well-known/oauth-protected-resource[/<path>] (rewritten here) — RFC
// 9728 metadata for the resource at <path>, or for the agent API at the root
// (docs/specs/MCP-PROTOCOL.md §6.1).
export function GET(req: Request) {
  const path = new URL(req.url).searchParams.get("path");
  const resource = path
    ? resources().find(
        (r) => new URL(r).pathname === `/${path.replace(/^\/+|\/+$/g, "")}`,
      )
    : resources()[0];

  if (!resource)
    return oauthError("not_found", "no protected resource at that path", 404);

  return oauthJson({
    resource,
    authorization_servers: [issuer()],
    scopes_supported: DEFAULT_SCOPES,
    bearer_methods_supported: ["header"],
    resource_name: "MERFISHEYES agent API",
    resource_documentation: `${issuer()}/docs`,
    scopes_documentation: SCOPES,
  });
}

export const OPTIONS = preflight;
export const dynamic = "force-dynamic";
