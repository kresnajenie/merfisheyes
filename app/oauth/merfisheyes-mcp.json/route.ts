import { issuer } from "@/lib/oauth/config";

/**
 * GET /oauth/merfisheyes-mcp.json — the Client ID Metadata Document of our
 * own agent package, merfisheyes-mcp (packages/merfisheyes-mcp): its
 * client_id is this URL (docs/specs/MCP-PROTOCOL.md §6.2), so the consent
 * page names the site as its publisher instead of "registered itself". It
 * runs on the person's computer, hence a loopback redirect on any port.
 */
export function GET() {
  const base = issuer();

  return Response.json(
    {
      client_id: `${base}/oauth/merfisheyes-mcp.json`,
      client_name: "MERFISHEYES MCP",
      client_uri: `${base}/docs`,
      redirect_uris: ["http://127.0.0.1/callback"],
      grant_types: ["authorization_code", "refresh_token"],
      response_types: ["code"],
      token_endpoint_auth_method: "none",
    },
    { headers: { "Cache-Control": "public, max-age=3600" } },
  );
}

export const dynamic = "force-dynamic";
