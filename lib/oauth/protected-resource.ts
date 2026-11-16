import { DEFAULT_SCOPES, SCOPES, issuer, resources } from "@/lib/oauth/config";
import { oauthError, oauthJson } from "@/lib/oauth/http";

/**
 * RFC 9728 metadata for the protected resource at `path` (the part after
 * `/.well-known/oauth-protected-resource`), or for the agent API when there
 * is none (docs/specs/MCP-PROTOCOL.md §6.1). A client checks that
 * `resource` is exactly the URL it is connecting to, so an unknown path is
 * 404, never another resource's metadata.
 */
export function protectedResourceMetadata(path: string | null): Response {
  const clean = (path ?? "").replace(/^\/+|\/+$/g, "");
  const resource = clean
    ? resources().find((r) => new URL(r).pathname === `/${clean}`)
    : resources()[0];

  if (!resource)
    return oauthError("not_found", "no protected resource at that path", 404);

  return oauthJson({
    resource,
    authorization_servers: [issuer()],
    scopes_supported: DEFAULT_SCOPES,
    bearer_methods_supported: ["header"],
    resource_name: resource.endsWith("/mcp")
      ? "MERFISHEYES (MCP connector)"
      : "MERFISHEYES agent API",
    resource_documentation: `${issuer()}/docs`,
    scopes_documentation: SCOPES,
  });
}
