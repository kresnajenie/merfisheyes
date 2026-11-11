import { prisma } from "@/lib/prisma";
import {
  DEFAULT_SCOPES,
  LIFETIME,
  PUBLIC_SCOPE,
  type Scope,
  issuer,
  parseScopes,
  resources,
} from "@/lib/oauth/config";
import {
  ClientError,
  getClient,
  isLoopback,
  redirectMatches,
} from "@/lib/oauth/clients";
import {
  isValidChallenge,
  newSecret,
  pkceMatches,
  sha256,
} from "@/lib/oauth/crypto";

/**
 * The authorization code flow with PKCE and resource indicators, and the
 * refresh-token rotation (docs/specs/MCP-PROTOCOL.md §6.4). Route handlers
 * parse HTTP and call these; the rules live here.
 */

export class OAuthError extends Error {
  constructor(
    readonly code: string,
    message: string,
    readonly status = 400,
  ) {
    super(message);
  }
}

export interface AuthorizeRequest {
  clientId: string;
  clientName: string;
  clientKind: "registered" | "metadata_document";
  /** Host of a metadata-document client_id; null for a registered client. */
  clientHost: string | null;
  redirectUri: string;
  redirectHost: string;
  loopbackOnly: boolean;
  scopes: Scope[];
  resource: string;
  state: string | null;
  codeChallenge: string;
}

/** A request that cannot be answered by redirecting: the client or the
 * redirect URI is not trusted, so the person sees the error instead
 * (RFC 6749 §4.1.2.1). */
export class UntrustedRequestError extends Error {}

/** Where an error goes once the redirect URI is trusted: back to the client,
 * with state and iss (RFC 9207). */
export function errorRedirect(
  redirectUri: string,
  state: string | null,
  code: string,
  description: string,
) {
  const u = new URL(redirectUri);

  u.searchParams.set("error", code);
  u.searchParams.set("error_description", description);
  if (state) u.searchParams.set("state", state);
  u.searchParams.set("iss", issuer());

  return u.toString();
}

/** Validate /oauth/authorize's query. Throws UntrustedRequestError (show it)
 * or returns { request } / { redirect } (an error sent back to the client). */
export async function validateAuthorize(
  q: URLSearchParams,
): Promise<{ request: AuthorizeRequest } | { redirect: string }> {
  const clientId = q.get("client_id") ?? "";
  const redirectUri = q.get("redirect_uri") ?? "";
  let client;

  try {
    client = await getClient(clientId);
  } catch (e) {
    throw new UntrustedRequestError(
      e instanceof ClientError ? e.message : "unknown client",
    );
  }
  if (!redirectUri || !redirectMatches(client.redirectUris, redirectUri)) {
    throw new UntrustedRequestError(
      "the redirect_uri is not one this client registered",
    );
  }

  const state = q.get("state");
  const fail = (code: string, description: string) => ({
    redirect: errorRedirect(redirectUri, state, code, description),
  });

  if (q.get("response_type") !== "code")
    return fail("unsupported_response_type", "response_type must be code");
  const challenge = q.get("code_challenge") ?? "";

  if (
    q.get("code_challenge_method") !== "S256" ||
    !isValidChallenge(challenge)
  ) {
    return fail(
      "invalid_request",
      "PKCE is required: code_challenge with code_challenge_method S256",
    );
  }
  const { scopes: asked, unknown } = parseScopes(q.get("scope"));

  if (unknown.length)
    return fail("invalid_scope", `unknown scope(s): ${unknown.join(" ")}`);
  const scopes = [
    ...new Set<Scope>([
      PUBLIC_SCOPE,
      ...(asked.length ? asked : DEFAULT_SCOPES),
    ]),
  ];
  const resource = normaliseResource(q.get("resource"));

  if (!resource)
    return fail(
      "invalid_target",
      `resource must be one of: ${resources().join(" ")}`,
    );

  return {
    request: {
      clientId: client.id,
      clientName: client.name,
      clientKind: client.kind,
      clientHost:
        client.kind === "metadata_document" ? new URL(client.id).host : null,
      redirectUri,
      redirectHost: new URL(redirectUri).host,
      loopbackOnly: client.redirectUris.every(isLoopback),
      scopes,
      resource,
      state,
      codeChallenge: challenge,
    },
  };
}

/** The resource a request names, matched case-insensitively on scheme and
 * host and without a trailing slash; absent means the agent API. */
export function normaliseResource(raw: string | null): string | null {
  if (!raw) return resources()[0];
  let u: URL;

  try {
    u = new URL(raw);
  } catch {
    return null;
  }
  if (u.hash) return null;
  const canonical = `${u.protocol}//${u.host}${u.pathname}`.replace(/\/+$/, "");

  return resources().find((r) => r === canonical) ?? null;
}

export type Reach = { all: true } | { all: false; projectIds: string[] };

/** The person said yes: record the grant, issue a code, return where to send
 * the browser. Project ids must be the person's own. */
export async function approve(
  userId: string,
  req: AuthorizeRequest,
  reach: Reach,
): Promise<string> {
  let projectIds: string[] = [];

  if (!reach.all) {
    const owned = await prisma.project.findMany({
      where: { id: { in: reach.projectIds }, ownerId: userId },
      select: { id: true },
    });

    projectIds = owned.map((p) => p.id);
    if (projectIds.length === 0) {
      throw new OAuthError(
        "invalid_request",
        "choose at least one of your projects, or all of them",
      );
    }
  }
  const data = {
    scopes: req.scopes,
    allProjects: reach.all,
    projectIds,
    resource: req.resource,
    revokedAt: null,
  };
  const grant = await prisma.oAuthGrant.upsert({
    where: { userId_clientId: { userId, clientId: req.clientId } },
    create: { userId, clientId: req.clientId, ...data },
    update: data,
  });
  const code = newSecret("mec");

  await prisma.oAuthCode.create({
    data: {
      codeHash: sha256(code),
      grantId: grant.id,
      redirectUri: req.redirectUri,
      codeChallenge: req.codeChallenge,
      resource: req.resource,
      scopes: req.scopes,
      expiresAt: new Date(Date.now() + LIFETIME.codeSeconds * 1000),
    },
  });
  const u = new URL(req.redirectUri);

  u.searchParams.set("code", code);
  if (req.state) u.searchParams.set("state", req.state);
  u.searchParams.set("iss", issuer());

  return u.toString();
}

export interface TokenResponse {
  access_token: string;
  token_type: "Bearer";
  expires_in: number;
  refresh_token: string;
  scope: string;
}

async function issueTokens(
  grantId: string,
  resource: string,
  scopes: string[],
): Promise<TokenResponse> {
  const access = newSecret("mea");
  const refresh = newSecret("mer");
  const now = Date.now();

  await prisma.oAuthToken.createMany({
    data: [
      {
        grantId,
        kind: "access",
        tokenHash: sha256(access),
        resource,
        scopes,
        expiresAt: new Date(now + LIFETIME.accessSeconds * 1000),
      },
      {
        grantId,
        kind: "refresh",
        tokenHash: sha256(refresh),
        resource,
        scopes,
        expiresAt: new Date(now + LIFETIME.refreshSeconds * 1000),
      },
    ],
  });

  return {
    access_token: access,
    token_type: "Bearer",
    expires_in: LIFETIME.accessSeconds,
    refresh_token: refresh,
    scope: scopes.join(" "),
  };
}

/** Revoke a grant and everything issued under it. */
export async function revokeGrant(grantId: string) {
  const now = new Date();

  await prisma.$transaction([
    prisma.oAuthGrant.update({
      where: { id: grantId },
      data: { revokedAt: now },
    }),
    prisma.oAuthToken.updateMany({
      where: { grantId, revokedAt: null },
      data: { revokedAt: now },
    }),
  ]);
}

const invalidGrant = (why: string) => new OAuthError("invalid_grant", why);

/** grant_type=authorization_code */
export async function exchangeCode(p: {
  code: string;
  clientId: string;
  redirectUri: string;
  codeVerifier: string;
  resource: string | null;
}): Promise<TokenResponse> {
  const row = await prisma.oAuthCode.findUnique({
    where: { codeHash: sha256(p.code) },
    include: { grant: true },
  });

  if (!row) throw invalidGrant("unknown authorization code");
  if (row.grant.clientId !== p.clientId)
    throw invalidGrant("the code was issued to another client");
  if (row.usedAt) {
    // A code presented twice was intercepted, or replayed: revoke what it
    // produced (RFC 6749 §4.1.2).
    await revokeGrant(row.grantId);
    throw invalidGrant("the authorization code was already used");
  }
  if (row.expiresAt < new Date())
    throw invalidGrant("the authorization code expired");
  if (row.redirectUri !== p.redirectUri)
    throw invalidGrant("redirect_uri differs from the authorization request");
  if (!pkceMatches(p.codeVerifier, row.codeChallenge))
    throw invalidGrant("code_verifier does not match the code_challenge");
  if (p.resource !== null && normaliseResource(p.resource) !== row.resource) {
    throw new OAuthError(
      "invalid_target",
      "resource differs from the authorization request",
    );
  }
  if (row.grant.revokedAt) throw invalidGrant("the grant was revoked");
  const claimed = await prisma.oAuthCode.updateMany({
    where: { codeHash: row.codeHash, usedAt: null },
    data: { usedAt: new Date() },
  });

  if (claimed.count !== 1)
    throw invalidGrant("the authorization code was already used");

  return issueTokens(row.grantId, row.resource, row.scopes);
}

/** grant_type=refresh_token — rotates: the presented token is spent and a new
 * pair issued. A spent token presented again means it was stolen: the grant
 * is revoked (MCP authorization, token theft; OAuth 2.1 §4.3.1). */
export async function refresh(p: {
  refreshToken: string;
  clientId: string;
  scope: string | null;
  resource: string | null;
}): Promise<TokenResponse> {
  const row = await prisma.oAuthToken.findUnique({
    where: { tokenHash: sha256(p.refreshToken) },
    include: { grant: true },
  });

  if (!row || row.kind !== "refresh")
    throw invalidGrant("unknown refresh token");
  if (row.grant.clientId !== p.clientId)
    throw invalidGrant("the token was issued to another client");
  if (row.usedAt) {
    await revokeGrant(row.grantId);
    throw invalidGrant(
      "the refresh token was already used; the grant is revoked — authorize again",
    );
  }
  if (row.revokedAt || row.grant.revokedAt)
    throw invalidGrant("the grant was revoked");
  if (row.expiresAt < new Date())
    throw invalidGrant("the refresh token expired");
  if (p.resource !== null && normaliseResource(p.resource) !== row.resource) {
    throw new OAuthError("invalid_target", "resource differs from the grant");
  }
  let scopes = row.scopes;

  if (p.scope) {
    const { scopes: asked, unknown } = parseScopes(p.scope);

    if (unknown.length || asked.some((s) => !row.scopes.includes(s))) {
      throw new OAuthError(
        "invalid_scope",
        "a refresh can only narrow the granted scopes",
      );
    }
    scopes = [...new Set([PUBLIC_SCOPE, ...asked])];
  }
  const claimed = await prisma.oAuthToken.updateMany({
    where: { id: row.id, usedAt: null },
    data: { usedAt: new Date() },
  });

  if (claimed.count !== 1) {
    await revokeGrant(row.grantId);
    throw invalidGrant(
      "the refresh token was already used; the grant is revoked — authorize again",
    );
  }

  return issueTokens(row.grantId, row.resource, scopes);
}

/** RFC 7009: revoking any token of a grant ends that client's session —
 * every token of the grant. Unknown tokens are not an error. */
export async function revokeToken(token: string) {
  const row = await prisma.oAuthToken.findUnique({
    where: { tokenHash: sha256(token) },
  });

  if (!row) return;
  await prisma.oAuthToken.updateMany({
    where: { grantId: row.grantId, revokedAt: null },
    data: { revokedAt: new Date() },
  });
}

export interface VerifiedToken {
  grantId: string;
  clientId: string;
  clientName: string;
  userId: string;
  role: string;
  scopes: string[];
  projects: "all" | string[];
  expiresAt: Date;
}

/** An access token presented to `resource`, or null. The reach (projects)
 * is read from the grant, so narrowing it applies at once. */
export async function verifyAccessToken(
  token: string,
  resource: string,
): Promise<VerifiedToken | null> {
  const row = await prisma.oAuthToken.findUnique({
    where: { tokenHash: sha256(token) },
    include: {
      grant: {
        include: {
          user: { select: { role: true } },
          client: { select: { name: true } },
        },
      },
    },
  });
  const now = new Date();

  if (!row || row.kind !== "access" || row.revokedAt || row.expiresAt <= now)
    return null;
  if (row.grant.revokedAt || row.resource !== resource) return null;
  if (
    !row.grant.lastUsedAt ||
    now.getTime() - row.grant.lastUsedAt.getTime() > 5 * 60 * 1000
  ) {
    await prisma.oAuthGrant.update({
      where: { id: row.grantId },
      data: { lastUsedAt: now },
    });
  }

  return {
    grantId: row.grantId,
    clientId: row.grant.clientId,
    clientName: row.grant.client.name,
    userId: row.grant.userId,
    role: row.grant.user.role,
    scopes: row.scopes,
    projects: row.grant.allProjects ? "all" : row.grant.projectIds,
    expiresAt: row.expiresAt,
  };
}

/** A project an agent created joins a project-limited grant, so the agent
 * can go on working in it (docs/specs/MCP-PROTOCOL.md §3.2). */
export async function addProjectToGrant(grantId: string, projectId: string) {
  const grant = await prisma.oAuthGrant.findUnique({ where: { id: grantId } });

  if (!grant || grant.allProjects || grant.projectIds.includes(projectId))
    return;
  await prisma.oAuthGrant.update({
    where: { id: grantId },
    data: { projectIds: { push: projectId } },
  });
}
