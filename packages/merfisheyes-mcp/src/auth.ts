import { spawn } from "node:child_process";
import { createHash, randomBytes } from "node:crypto";
import { createServer, type Server } from "node:http";

import {
  API_BASE,
  DEFAULT_SCOPES,
  SIGN_IN_TIMEOUT_MS,
  SITE,
} from "./config.js";
import { type Saved, load, store } from "./credentials.js";

/**
 * Signing in to MERFISHEYES as an OAuth 2.1 public client, as the MCP
 * authorization spec asks of a client (docs/specs/MCP-PROTOCOL.md §6 in the
 * MERFISHEYES repository): discover the authorization server from the
 * resource's metadata, check its issuer, identify by the site's Client ID
 * Metadata Document (or register, where the server cannot fetch one), PKCE
 * S256, the resource indicator, a loopback redirect on a free port, the
 * `iss` of the answer checked, refresh tokens rotated.
 */

interface ServerMetadata {
  issuer: string;
  authorization_endpoint: string;
  token_endpoint: string;
  registration_endpoint?: string;
  revocation_endpoint?: string;
  code_challenge_methods_supported?: string[];
  client_id_metadata_document_supported?: boolean;
}

interface Discovered {
  resource: string;
  server: ServerMetadata;
}

let discovered: Promise<Discovered> | null = null;

async function getJson(url: string): Promise<Record<string, unknown>> {
  const res = await fetch(url, {
    headers: { accept: "application/json" },
    signal: AbortSignal.timeout(15_000),
  });

  if (!res.ok) throw new Error(`${url} answered ${res.status}`);

  return (await res.json()) as Record<string, unknown>;
}

/** The resource and its authorization server, checked (RFC 9728, RFC 8414). */
export function discover(): Promise<Discovered> {
  discovered ??= (async () => {
    const path = new URL(API_BASE).pathname;
    const prm = await getJson(
      `${SITE}/.well-known/oauth-protected-resource${path}`,
    );
    const resource = String(prm.resource ?? "");
    const issuer = String(
      (prm.authorization_servers as string[] | undefined)?.[0] ?? "",
    );

    if (resource !== API_BASE)
      throw new Error(
        `the site names its resource ${resource}, not ${API_BASE}`,
      );
    if (!issuer) throw new Error("the site names no authorization server");
    const server = (await getJson(
      `${issuer}/.well-known/oauth-authorization-server`,
    )) as unknown as ServerMetadata;

    // A document that names another issuer than the one asked must not be used.
    if (server.issuer !== issuer)
      throw new Error(
        `the authorization server says it is ${server.issuer}, not ${issuer}`,
      );
    if (!server.code_challenge_methods_supported?.includes("S256")) {
      throw new Error(
        "the authorization server does not support PKCE S256; refusing to sign in",
      );
    }

    return { resource, server };
  })().catch((e) => {
    discovered = null;
    throw e;
  });

  return discovered;
}

const REDIRECT_PATH = "/callback";

/** The client id: the site's metadata document over HTTPS, or a registered
 * one (a local http:// deployment cannot be fetched as a document). */
async function clientIdFor(server: ServerMetadata): Promise<string> {
  const saved = load(server.issuer);

  if (saved?.clientId) return saved.clientId;
  let clientId: string;

  if (
    server.client_id_metadata_document_supported &&
    server.issuer.startsWith("https://")
  ) {
    clientId = `${server.issuer}/oauth/merfisheyes-mcp.json`;
  } else if (server.registration_endpoint) {
    const res = await fetch(server.registration_endpoint, {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({
        client_name: "MERFISHEYES MCP",
        redirect_uris: [`http://127.0.0.1${REDIRECT_PATH}`],
        grant_types: ["authorization_code", "refresh_token"],
        token_endpoint_auth_method: "none",
        application_type: "native",
      }),
      signal: AbortSignal.timeout(15_000),
    });

    if (!res.ok)
      throw new Error(
        `registration answered ${res.status}: ${await res.text()}`,
      );
    clientId = String(((await res.json()) as { client_id: string }).client_id);
  } else {
    throw new Error(
      "the authorization server offers no way for this client to identify itself",
    );
  }
  store(server.issuer, { clientId });

  return clientId;
}

interface TokenAnswer {
  access_token: string;
  refresh_token?: string;
  expires_in?: number;
  scope?: string;
}

async function tokenRequest(
  server: ServerMetadata,
  form: Record<string, string>,
): Promise<TokenAnswer> {
  const res = await fetch(server.token_endpoint, {
    method: "POST",
    headers: {
      "content-type": "application/x-www-form-urlencoded",
      accept: "application/json",
    },
    body: new URLSearchParams(form),
    signal: AbortSignal.timeout(15_000),
  });
  const body = (await res.json().catch(() => ({}))) as Record<string, string>;

  if (!res.ok)
    throw new Error(
      body.error_description ||
        body.error ||
        `the token endpoint answered ${res.status}`,
    );

  return body as unknown as TokenAnswer;
}

function keep(
  server: ServerMetadata,
  clientId: string,
  t: TokenAnswer,
  previous?: Saved,
) {
  store(server.issuer, {
    clientId,
    accessToken: t.access_token,
    refreshToken: t.refresh_token ?? previous?.refreshToken,
    expiresAt: Date.now() + (t.expires_in ?? 3600) * 1000,
    scope: t.scope ?? previous?.scope,
  });
}

/** Open a URL in the person's browser; quietly does nothing where none can be
 * opened (the link is also returned to the agent to show). */
export function openBrowser(url: string): void {
  if (process.env.MERFISHEYES_NO_BROWSER) return;
  const [cmd, args] =
    process.platform === "darwin"
      ? ["open", [url]]
      : process.platform === "win32"
        ? ["cmd", ["/c", "start", "", url]]
        : ["xdg-open", [url]];

  try {
    spawn(cmd, args, { stdio: "ignore", detached: true })
      .on("error", () => {})
      .unref();
  } catch {
    // no browser here: the agent shows the link instead
  }
}

export interface SignIn {
  url: string;
  done: Promise<void>;
}

let pending: (SignIn & { scopes: string }) | null = null;

/**
 * Start signing in: a one-shot listener on 127.0.0.1 waits for the browser to
 * come back with a code, which is swapped for tokens and kept. Returns the
 * link to open at once; `done` settles when the person has answered (or after
 * SIGN_IN_TIMEOUT_MS). Asking again while one is open returns the same link.
 */
export async function startSignIn(
  scopes: string[] = DEFAULT_SCOPES,
): Promise<SignIn> {
  const scope = [...new Set(scopes)].join(" ");

  if (pending && pending.scopes === scope) return pending;
  const { resource, server } = await discover();
  const clientId = await clientIdFor(server);
  const verifier = randomBytes(32).toString("base64url");
  const state = randomBytes(16).toString("base64url");
  let listener: Server;
  const port = await new Promise<number>((resolve, reject) => {
    listener = createServer();
    listener.on("error", reject);
    listener.listen(0, "127.0.0.1", () =>
      resolve((listener.address() as { port: number }).port),
    );
  });
  const redirectUri = `http://127.0.0.1:${port}${REDIRECT_PATH}`;
  const url = new URL(server.authorization_endpoint);

  url.search = new URLSearchParams({
    response_type: "code",
    client_id: clientId,
    redirect_uri: redirectUri,
    code_challenge: createHash("sha256").update(verifier).digest("base64url"),
    code_challenge_method: "S256",
    scope,
    resource,
    state,
  }).toString();

  const done = new Promise<void>((resolve, reject) => {
    const finish = (error?: Error) => {
      clearTimeout(timer);
      listener.close();
      pending = null;
      if (error) reject(error);
      else resolve();
    };
    const timer = setTimeout(
      () =>
        finish(
          new Error("nobody answered the sign-in link in time; start again"),
        ),
      SIGN_IN_TIMEOUT_MS,
    );

    timer.unref();
    listener.on("request", (req, res) => {
      const answer = new URL(req.url ?? "/", redirectUri);

      if (answer.pathname !== REDIRECT_PATH) {
        res.writeHead(404).end();

        return;
      }
      const page = (status: number, text: string) =>
        res
          .writeHead(status, { "content-type": "text/html; charset=utf-8" })
          .end(
            `<!doctype html><title>MERFISHEYES</title><p style="font:16px system-ui;margin:3em">${text}</p>`,
          );

      // Another tab or site cannot finish this sign-in: the state must match,
      // and the issuer the answer names must be the one discovered (RFC 9207).
      if (answer.searchParams.get("state") !== state) {
        page(
          400,
          "This answer is not for the sign-in that is waiting. Close this tab.",
        );

        return;
      }
      if (answer.searchParams.get("iss") !== server.issuer) {
        page(
          400,
          "This answer came from an unexpected server. Close this tab.",
        );
        finish(
          new Error(
            "the sign-in answer named another issuer; nothing was kept",
          ),
        );

        return;
      }
      const error = answer.searchParams.get("error");

      if (error) {
        page(200, "Not signed in. You can close this tab.");
        finish(
          new Error(
            error === "access_denied"
              ? "you declined; nothing was shared"
              : `sign-in failed: ${error}`,
          ),
        );

        return;
      }
      tokenRequest(server, {
        grant_type: "authorization_code",
        code: answer.searchParams.get("code") ?? "",
        redirect_uri: redirectUri,
        client_id: clientId,
        code_verifier: verifier,
        resource,
      })
        .then((t) => {
          keep(server, clientId, t);
          page(
            200,
            "Signed in to MERFISHEYES. You can close this tab and go back to your agent.",
          );
          finish();
        })
        .catch((e: Error) => {
          page(500, "Signing in failed. Go back to your agent and try again.");
          finish(e);
        });
    });
  });

  done.catch(() => {}); // a sign-in nobody waits on must not crash the server
  pending = { url: url.toString(), done, scopes: scope };

  return pending;
}

/** A usable access token, refreshed when it has (nearly) expired; null when
 * signed out. `force` refreshes even an unexpired one (after a 401). */
export async function accessToken(force = false): Promise<string | null> {
  const { resource, server } = await discover();
  const saved = load(server.issuer);

  if (!saved?.accessToken) return null;
  if (!force && saved.expiresAt && saved.expiresAt - 60_000 > Date.now())
    return saved.accessToken;
  if (!saved.refreshToken) return null;
  try {
    const t = await tokenRequest(server, {
      grant_type: "refresh_token",
      refresh_token: saved.refreshToken,
      client_id: saved.clientId,
      resource,
    });

    keep(server, saved.clientId, t, saved);

    return t.access_token;
  } catch {
    // Spent, expired or revoked: forget the tokens, keep the client id.
    store(server.issuer, { clientId: saved.clientId });

    return null;
  }
}

/** Sign out: revoke the grant's tokens at the server and forget them here. */
export async function signOut(): Promise<boolean> {
  const { server } = await discover();
  const saved = load(server.issuer);

  if (!saved?.accessToken && !saved?.refreshToken) return false;
  if (server.revocation_endpoint) {
    await fetch(server.revocation_endpoint, {
      method: "POST",
      headers: { "content-type": "application/x-www-form-urlencoded" },
      body: new URLSearchParams({
        token: saved.refreshToken ?? saved.accessToken ?? "",
        client_id: saved.clientId,
      }),
      signal: AbortSignal.timeout(15_000),
    }).catch(() => {});
  }
  store(server.issuer, { clientId: saved.clientId });

  return true;
}
