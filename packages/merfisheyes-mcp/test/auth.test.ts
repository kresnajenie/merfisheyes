import { mkdtempSync, readFileSync, statSync } from "node:fs";
import { createServer, type Server } from "node:http";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { afterAll, beforeAll, describe, expect, it } from "vitest";

/**
 * Signing in as an MCP client against a stand-in authorization server:
 * discovery and its checks, registration (an http:// issuer cannot use a
 * metadata document), PKCE, the loopback answer with its state and `iss`
 * checks, refresh rotation, sign-out — and where the tokens are kept.
 */

let site: Server;
let base = "";
const seen: { authorize?: URL; tokenForms: URLSearchParams[]; revoked: string[] } = { tokenForms: [], revoked: [] };
let issuerOverride: string | null = null;

beforeAll(async () => {
  site = createServer((req, res) => {
    const url = new URL(req.url ?? "/", base);
    const json = (status: number, body: unknown) => res.writeHead(status, { "content-type": "application/json" }).end(JSON.stringify(body));
    let body = "";

    req.on("data", (c) => (body += c));
    req.on("end", () => {
      switch (url.pathname) {
        case "/.well-known/oauth-protected-resource/api/agent/v1":
          return json(200, { resource: `${base}/api/agent/v1`, authorization_servers: [base] });
        case "/.well-known/oauth-authorization-server":
          return json(200, {
            issuer: issuerOverride ?? base,
            authorization_endpoint: `${base}/oauth/authorize`,
            token_endpoint: `${base}/api/oauth/token`,
            registration_endpoint: `${base}/api/oauth/register`,
            revocation_endpoint: `${base}/api/oauth/revoke`,
            code_challenge_methods_supported: ["S256"],
            client_id_metadata_document_supported: true,
          });
        case "/api/oauth/register":
          return json(201, { client_id: "mcl_test" });
        case "/api/oauth/token": {
          const form = new URLSearchParams(body);

          seen.tokenForms.push(form);
          const n = seen.tokenForms.length;

          return json(200, { access_token: `mea_${n}`, refresh_token: `mer_${n}`, expires_in: 3600, scope: "catalog:read projects:read" });
        }
        case "/api/oauth/revoke":
          seen.revoked.push(new URLSearchParams(body).get("token") ?? "");

          return json(200, {});
        default:
          return json(404, {});
      }
    });
  });
  await new Promise<void>((r) => site.listen(0, "127.0.0.1", r));
  base = `http://127.0.0.1:${(site.address() as { port: number }).port}`;
  process.env.MERFISHEYES_URL = base;
  process.env.MERFISHEYES_CREDENTIALS = join(mkdtempSync(join(tmpdir(), "meyes-")), "sub", "credentials.json");
  process.env.MERFISHEYES_NO_BROWSER = "1";
});

afterAll(() => site.close());

const auth = () => import("../src/auth.js");

describe("signing in", () => {
  it("registers, sends PKCE and the resource, and keeps the tokens privately", async () => {
    const { startSignIn } = await auth();
    const s = await startSignIn();
    const url = new URL(s.url);

    expect(url.origin + url.pathname).toBe(`${base}/oauth/authorize`);
    expect(url.searchParams.get("client_id")).toBe("mcl_test"); // http:// issuer: no metadata document
    expect(url.searchParams.get("code_challenge_method")).toBe("S256");
    expect(url.searchParams.get("resource")).toBe(`${base}/api/agent/v1`);
    expect(url.searchParams.get("scope")).toBe("catalog:read projects:read");
    const redirect = url.searchParams.get("redirect_uri")!;

    expect(redirect).toMatch(/^http:\/\/127\.0\.0\.1:\d+\/callback$/);

    // A forged answer (wrong state) is refused and the sign-in keeps waiting.
    const forged = await fetch(`${redirect}?code=x&state=nope&iss=${encodeURIComponent(base)}`);

    expect(forged.status).toBe(400);

    const ok = await fetch(`${redirect}?code=the-code&state=${url.searchParams.get("state")}&iss=${encodeURIComponent(base)}`);

    expect(ok.status).toBe(200);
    await s.done;
    const exchange = seen.tokenForms.at(-1)!;

    expect(exchange.get("grant_type")).toBe("authorization_code");
    expect(exchange.get("code")).toBe("the-code");
    expect(exchange.get("redirect_uri")).toBe(redirect);
    expect(exchange.get("code_verifier")).toMatch(/^[A-Za-z0-9_-]{43}$/);

    const path = process.env.MERFISHEYES_CREDENTIALS!;

    expect(statSync(path).mode & 0o777).toBe(0o600);
    expect(JSON.parse(readFileSync(path, "utf8"))[base]).toMatchObject({ clientId: "mcl_test", accessToken: "mea_1", refreshToken: "mer_1" });
  });

  it("refreshes with rotation, then signs out by revoking", async () => {
    const { accessToken, signOut } = await auth();

    expect(await accessToken()).toBe("mea_1");
    expect(await accessToken(true)).toBe("mea_2");
    expect(seen.tokenForms.at(-1)!.get("refresh_token")).toBe("mer_1");
    expect(JSON.parse(readFileSync(process.env.MERFISHEYES_CREDENTIALS!, "utf8"))[base].refreshToken).toBe("mer_2");
    expect(await signOut()).toBe(true);
    expect(seen.revoked).toEqual(["mer_2"]);
    expect(await accessToken()).toBeNull();
  });

  it("refuses an answer that names another issuer", async () => {
    const { startSignIn } = await auth();
    const s = await startSignIn(["catalog:read", "projects:read", "projects:write"]);
    const url = new URL(s.url);
    const res = await fetch(`${url.searchParams.get("redirect_uri")}?code=c&state=${url.searchParams.get("state")}&iss=https%3A%2F%2Fevil.example`);

    expect(res.status).toBe(400);
    await expect(s.done).rejects.toThrow(/another issuer/);
  });
});
