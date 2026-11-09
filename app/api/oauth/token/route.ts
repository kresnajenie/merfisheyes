import { OAuthError, exchangeCode, refresh } from "@/lib/oauth/flow";
import { formParams, oauthError, oauthJson, preflight } from "@/lib/oauth/http";

// POST /api/oauth/token — the code exchange and the refresh
// (docs/specs/MCP-PROTOCOL.md §6.4). Public clients: client_id, no secret.
export async function POST(req: Request) {
  const p = await formParams(req);
  const clientId = p.get("client_id") ?? "";

  if (!clientId)
    return oauthError("invalid_client", "client_id is required", 401);
  try {
    switch (p.get("grant_type")) {
      case "authorization_code": {
        const code = p.get("code");
        const redirectUri = p.get("redirect_uri");
        const codeVerifier = p.get("code_verifier");

        if (!code || !redirectUri || !codeVerifier) {
          return oauthError(
            "invalid_request",
            "code, redirect_uri and code_verifier are required",
          );
        }

        return oauthJson(
          await exchangeCode({
            code,
            clientId,
            redirectUri,
            codeVerifier,
            resource: p.get("resource"),
          }),
        );
      }
      case "refresh_token": {
        const refreshToken = p.get("refresh_token");

        if (!refreshToken)
          return oauthError("invalid_request", "refresh_token is required");

        return oauthJson(
          await refresh({
            refreshToken,
            clientId,
            scope: p.get("scope"),
            resource: p.get("resource"),
          }),
        );
      }
      default:
        return oauthError(
          "unsupported_grant_type",
          "grant_type is authorization_code or refresh_token",
        );
    }
  } catch (e) {
    if (e instanceof OAuthError) return oauthError(e.code, e.message, e.status);
    console.error("oauth token:", e);

    return oauthError(
      "server_error",
      "the token could not be issued; try again",
      500,
    );
  }
}

export const OPTIONS = preflight;
