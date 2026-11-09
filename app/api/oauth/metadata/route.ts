import { ALL_SCOPES, issuer } from "@/lib/oauth/config";
import { oauthJson, preflight } from "@/lib/oauth/http";

// GET /.well-known/oauth-authorization-server (rewritten here) — RFC 8414
// metadata (docs/specs/MCP-PROTOCOL.md §6.1).
export function GET() {
  const base = issuer();

  return oauthJson({
    issuer: base,
    authorization_endpoint: `${base}/oauth/authorize`,
    token_endpoint: `${base}/api/oauth/token`,
    registration_endpoint: `${base}/api/oauth/register`,
    revocation_endpoint: `${base}/api/oauth/revoke`,
    response_types_supported: ["code"],
    response_modes_supported: ["query"],
    grant_types_supported: ["authorization_code", "refresh_token"],
    code_challenge_methods_supported: ["S256"],
    token_endpoint_auth_methods_supported: ["none"],
    revocation_endpoint_auth_methods_supported: ["none"],
    scopes_supported: ALL_SCOPES,
    client_id_metadata_document_supported: true,
    authorization_response_iss_parameter_supported: true,
    service_documentation: `${base}/docs`,
  });
}

export const OPTIONS = preflight;
export const dynamic = "force-dynamic";
