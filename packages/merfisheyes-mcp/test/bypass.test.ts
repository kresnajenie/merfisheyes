import { mkdtempSync } from "node:fs";
import { createServer, type Server } from "node:http";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { afterAll, beforeAll, expect, it } from "vitest";

/**
 * A site behind Vercel's Deployment Protection (staging): with
 * MERFISHEYES_VERCEL_BYPASS set, every request carries the bypass header;
 * without it, the redirect to Vercel's login is reported as such.
 */

let site: Server;
let base = "";
const headersSeen: Record<string, string | undefined>[] = [];
let protect = false;

beforeAll(async () => {
  site = createServer((req, res) => {
    const url = new URL(req.url ?? "/", base);

    headersSeen.push({ path: url.pathname, bypass: req.headers["x-vercel-protection-bypass"] as string | undefined });
    if (protect && !req.headers["x-vercel-protection-bypass"]) {
      res.writeHead(302, { location: "https://vercel.com/sso-api?url=x" }).end();

      return;
    }
    const json = (b: unknown) => res.writeHead(200, { "content-type": "application/json" }).end(JSON.stringify(b));

    if (url.pathname === "/.well-known/oauth-protected-resource/api/agent/v1") {
      return json({ resource: `${base}/api/agent/v1`, authorization_servers: [base] });
    }
    if (url.pathname === "/.well-known/oauth-authorization-server") {
      return json({
        issuer: base,
        authorization_endpoint: `${base}/oauth/authorize`,
        token_endpoint: `${base}/api/oauth/token`,
        registration_endpoint: `${base}/api/oauth/register`,
        code_challenge_methods_supported: ["S256"],
      });
    }
    if (url.pathname === "/api/oauth/register") return json({ client_id: "mcl_bypass" });

    return json({ species: [], tissues: [], platforms: [] });
  });
  await new Promise<void>((r) => site.listen(0, "127.0.0.1", r));
  base = `http://127.0.0.1:${(site.address() as { port: number }).port}`;
  process.env.MERFISHEYES_URL = base;
  process.env.MERFISHEYES_CREDENTIALS = join(mkdtempSync(join(tmpdir(), "meyes-bp-")), "credentials.json");
  process.env.MERFISHEYES_NO_BROWSER = "1";
  process.env.MERFISHEYES_VERCEL_BYPASS = "s3cret";
});

afterAll(() => site.close());

it("sends the bypass header on discovery, registration and API calls", async () => {
  protect = true;
  const { startSignIn } = await import("../src/auth.js");
  const { apiGet } = await import("../src/api.js");

  await startSignIn();
  await apiGet("/catalog/facets", {});
  const paths = headersSeen.map((h) => h.path);

  expect(paths).toEqual(
    expect.arrayContaining([
      "/.well-known/oauth-protected-resource/api/agent/v1",
      "/.well-known/oauth-authorization-server",
      "/api/oauth/register",
      "/api/agent/v1/catalog/facets",
    ]),
  );
  expect(headersSeen.every((h) => h.bypass === "s3cret")).toBe(true);
});
