/**
 * The authorization server's fixed facts (docs/specs/MCP-PROTOCOL.md §6).
 * The issuer is the site origin, exactly as clients compare it: no trailing
 * slash, no path. Each deployment needs its own (OAUTH_ISSUER, else
 * NEXT_PUBLIC_APP_URL): a token from staging must not be accepted by
 * production, and a client checks that the issuer it discovered is the one
 * that redirected back (RFC 9207).
 */

export function issuer(): string {
  const raw =
    process.env.OAUTH_ISSUER ||
    process.env.NEXT_PUBLIC_APP_URL ||
    "http://localhost:3000";

  return raw.replace(/\/+$/, "");
}

/** The protected resources this server issues tokens for (RFC 8707). The
 * agent API now; the remote MCP endpoint later (plan step 9). */
export function resources(): string[] {
  return [`${issuer()}/api/agent/v1`];
}

export function resourceMetadataUrl(resource = resources()[0]): string {
  const path = new URL(resource).pathname;

  return `${issuer()}/.well-known/oauth-protected-resource${path === "/" ? "" : path}`;
}

export const SCOPES = {
  "catalog:read": "Search and read the public catalog",
  "projects:read": "See your projects and datasets",
  "projects:write": "Create, edit and delete your projects",
  "datasets:write": "Edit, submit and delete your datasets",
  "viewer:control": "Control a MERFISHEYES tab you connect",
  export: "Save tables from MERFISHEYES to your computer",
} as const;

export type Scope = keyof typeof SCOPES;

export const ALL_SCOPES = Object.keys(SCOPES) as Scope[];

/** What a client gets when it asks for nothing (and what the resource
 * metadata advertises): the minimum for basic use; more is asked for when a
 * tool needs it (step-up, MCP authorization §scope challenge). */
export const DEFAULT_SCOPES: Scope[] = ["catalog:read", "projects:read"];

/** Everyone holds it, token or not; consent never lists it. */
export const PUBLIC_SCOPE: Scope = "catalog:read";

export const LIFETIME = {
  codeSeconds: 5 * 60,
  accessSeconds: 60 * 60,
  refreshSeconds: 30 * 24 * 60 * 60,
  /** How long a fetched Client ID Metadata Document is trusted at most. */
  clientMetadataSeconds: 24 * 60 * 60,
  /** How long the signed consent form stays valid. */
  consentSeconds: 10 * 60,
};

export function isScope(s: string): s is Scope {
  return s in SCOPES;
}

/** Parse a space-separated scope string; unknown scopes are an error the
 * caller reports (invalid_scope), not silently dropped. */
export function parseScopes(raw: string | null | undefined): {
  scopes: Scope[];
  unknown: string[];
} {
  const words = (raw ?? "").split(/\s+/).filter(Boolean);
  const scopes = [...new Set(words.filter(isScope))];

  return { scopes, unknown: words.filter((w) => !isScope(w)) };
}
