import { randomBytes } from "crypto";

import { prisma } from "@/lib/prisma";
import { LIFETIME } from "@/lib/oauth/config";
import {
  LimitedFetchError,
  limitedFetch,
  maxAgeSeconds,
} from "@/lib/net/limited-fetch";

/**
 * Agent applications (docs/specs/MCP-PROTOCOL.md §6.2). A client is either
 * identified by the HTTPS URL of its Client ID Metadata Document — fetched,
 * checked and cached here — or registered dynamically (RFC 7591, deprecated
 * in MCP, kept for older clients). Every client is public: no secret.
 */

export interface Client {
  id: string;
  kind: "registered" | "metadata_document";
  name: string;
  uri: string | null;
  redirectUris: string[];
}

export class ClientError extends Error {
  constructor(
    readonly code:
      | "invalid_client"
      | "invalid_client_metadata"
      | "invalid_redirect_uri",
    message: string,
  ) {
    super(message);
  }
}

const LOOPBACK_HOSTS = new Set(["127.0.0.1", "[::1]", "localhost"]);

/** A redirect URI may be HTTPS, or HTTP on a loopback host (MCP
 * authorization, communication security). No fragment. */
export function isAllowedRedirectUri(raw: string): boolean {
  let u: URL;

  try {
    u = new URL(raw);
  } catch {
    return false;
  }
  if (u.hash || u.username || u.password) return false;
  if (u.protocol === "https:") return true;

  return u.protocol === "http:" && LOOPBACK_HOSTS.has(u.hostname);
}

export function isLoopback(raw: string): boolean {
  try {
    return LOOPBACK_HOSTS.has(new URL(raw).hostname);
  } catch {
    return false;
  }
}

/** Exact match, except that a loopback URI matches on any port: a native
 * client picks a free port per run (RFC 8252 §7.3). */
export function redirectMatches(
  registered: string[],
  presented: string,
): boolean {
  if (registered.includes(presented)) return true;
  if (!isLoopback(presented)) return false;
  const p = new URL(presented);

  return registered.some((r) => {
    if (!isLoopback(r)) return false;
    const q = new URL(r);

    return (
      q.protocol === p.protocol &&
      q.hostname === p.hostname &&
      q.pathname === p.pathname &&
      q.search === p.search
    );
  });
}

/** A client_id that names a metadata document: an HTTPS URL with a path. */
export function isMetadataDocumentId(clientId: string): boolean {
  try {
    const u = new URL(clientId);

    return (
      u.protocol === "https:" &&
      u.pathname !== "/" &&
      !u.hash &&
      !u.username &&
      !u.password
    );
  } catch {
    return false;
  }
}

function checkMetadata(doc: unknown, clientId: string) {
  if (!doc || typeof doc !== "object" || Array.isArray(doc)) {
    throw new ClientError(
      "invalid_client",
      "the client metadata document is not a JSON object",
    );
  }
  const d = doc as Record<string, unknown>;

  if (d.client_id !== clientId) {
    throw new ClientError(
      "invalid_client",
      "the document's client_id is not its URL",
    );
  }
  if (typeof d.client_name !== "string" || !d.client_name.trim()) {
    throw new ClientError("invalid_client", "the document has no client_name");
  }
  const method = d.token_endpoint_auth_method ?? "none";

  if (method !== "none") {
    throw new ClientError(
      "invalid_client",
      `token_endpoint_auth_method ${String(method)} is not supported; public clients only`,
    );
  }
  const uris = d.redirect_uris;

  if (
    !Array.isArray(uris) ||
    uris.length === 0 ||
    !uris.every((u) => typeof u === "string" && isAllowedRedirectUri(u))
  ) {
    throw new ClientError(
      "invalid_client",
      "redirect_uris must be HTTPS or loopback URLs",
    );
  }

  return {
    name: d.client_name.trim().slice(0, 200),
    uri: typeof d.client_uri === "string" ? d.client_uri.slice(0, 1000) : null,
    logoUri: typeof d.logo_uri === "string" ? d.logo_uri.slice(0, 1000) : null,
    redirectUris: uris as string[],
  };
}

async function fetchMetadataDocument(clientId: string): Promise<Client> {
  let res;

  try {
    res = await limitedFetch(clientId, {
      maxBytes: 64 * 1024,
      timeoutMs: 5000,
    });
  } catch (e) {
    throw new ClientError(
      "invalid_client",
      e instanceof LimitedFetchError
        ? `cannot fetch the client metadata document: ${e.message}`
        : "cannot fetch the client metadata document",
    );
  }
  if (res.status !== 200) {
    throw new ClientError(
      "invalid_client",
      `the client metadata document answered ${res.status}`,
    );
  }
  let doc: unknown;

  try {
    doc = JSON.parse(res.body);
  } catch {
    throw new ClientError(
      "invalid_client",
      "the client metadata document is not JSON",
    );
  }
  const meta = checkMetadata(doc, clientId);
  const ttl = maxAgeSeconds(
    res.headers,
    3600,
    300,
    LIFETIME.clientMetadataSeconds,
  );
  const row = {
    kind: "metadata_document",
    ...meta,
    expiresAt: new Date(Date.now() + ttl * 1000),
  };

  await prisma.oAuthClient.upsert({
    where: { id: clientId },
    create: { id: clientId, ...row },
    update: row,
  });

  return {
    id: clientId,
    kind: "metadata_document",
    name: meta.name,
    uri: meta.uri,
    redirectUris: meta.redirectUris,
  };
}

/** The client a client_id names, or a ClientError. */
export async function getClient(clientId: string): Promise<Client> {
  if (!clientId || clientId.length > 1000)
    throw new ClientError("invalid_client", "missing or overlong client_id");
  const row = await prisma.oAuthClient.findUnique({ where: { id: clientId } });

  if (isMetadataDocumentId(clientId)) {
    if (row && row.expiresAt && row.expiresAt > new Date()) {
      return {
        id: row.id,
        kind: "metadata_document",
        name: row.name,
        uri: row.uri,
        redirectUris: row.redirectUris,
      };
    }

    return fetchMetadataDocument(clientId);
  }
  if (!row || row.kind !== "registered")
    throw new ClientError("invalid_client", "unknown client_id");

  return {
    id: row.id,
    kind: "registered",
    name: row.name,
    uri: row.uri,
    redirectUris: row.redirectUris,
  };
}

/** Dynamic client registration (RFC 7591 §2–3). Returns the response body. */
export async function registerClient(body: unknown) {
  if (!body || typeof body !== "object")
    throw new ClientError("invalid_client_metadata", "expected a JSON object");
  const b = body as Record<string, unknown>;
  const uris = b.redirect_uris;

  if (!Array.isArray(uris) || uris.length === 0 || uris.length > 10) {
    throw new ClientError(
      "invalid_redirect_uri",
      "redirect_uris: one to ten URLs",
    );
  }
  if (!uris.every((u) => typeof u === "string" && isAllowedRedirectUri(u))) {
    throw new ClientError(
      "invalid_redirect_uri",
      "every redirect URI must be HTTPS or http on a loopback host",
    );
  }
  const method = b.token_endpoint_auth_method ?? "none";

  if (method !== "none") {
    throw new ClientError(
      "invalid_client_metadata",
      "only public clients (token_endpoint_auth_method none) are supported",
    );
  }
  const grantTypes = (b.grant_types as string[] | undefined) ?? [
    "authorization_code",
    "refresh_token",
  ];

  if (
    !Array.isArray(grantTypes) ||
    grantTypes.some((g) => g !== "authorization_code" && g !== "refresh_token")
  ) {
    throw new ClientError(
      "invalid_client_metadata",
      "grant_types: authorization_code and refresh_token only",
    );
  }
  const name =
    typeof b.client_name === "string" && b.client_name.trim()
      ? b.client_name.trim().slice(0, 200)
      : "Unnamed agent";
  const id = `mcl_${randomBytes(16).toString("base64url")}`;
  const row = await prisma.oAuthClient.create({
    data: {
      id,
      kind: "registered",
      name,
      uri:
        typeof b.client_uri === "string" ? b.client_uri.slice(0, 1000) : null,
      redirectUris: uris as string[],
    },
  });

  return {
    client_id: row.id,
    client_id_issued_at: Math.floor(row.createdAt.getTime() / 1000),
    client_name: row.name,
    redirect_uris: row.redirectUris,
    grant_types: grantTypes,
    response_types: ["code"],
    token_endpoint_auth_method: "none",
  };
}
