import { createHash, randomBytes } from "node:crypto";

import { type APIRequestContext, expect, test } from "@playwright/test";

import { signInDev } from "../../helpers/auth";

/**
 * The site as an MCP server (docs/plans/MCP.md pillar E): what a connector
 * does. Unsigned, `/mcp` answers 401 pointing at its own resource metadata;
 * the connector signs the person in for that resource (consent limited to
 * one project) and then calls tools over Streamable HTTP, which reach the
 * agent API in process as that person. A `/mcp` token is not an agent-API
 * token, and disconnecting stops it.
 *
 * Tagged @upload: it needs the suite's database (no S3). E2E_AGENT=1 runs it
 * against any migrated database the dev server's DATABASE_URL points at.
 */

const ENABLED = process.env.E2E_UPLOAD === "1" || process.env.E2E_AGENT === "1";
const REDIRECT = "http://127.0.0.1:5998/callback";

/** One JSON-RPC call to /mcp; the answer is JSON or one SSE event. */
async function rpc(
  request: APIRequestContext,
  token: string,
  method: string,
  params: unknown = {},
) {
  const res = await request.post("/mcp", {
    headers: {
      authorization: `Bearer ${token}`,
      "content-type": "application/json",
      accept: "application/json, text/event-stream",
    },
    data: { jsonrpc: "2.0", id: 1, method, params },
  });

  expect(res.status(), `${method}`).toBe(200);
  const text = await res.text();

  return JSON.parse(
    res.headers()["content-type"]?.includes("event-stream")
      ? /^data: (.*)$/m.exec(text)![1]
      : text,
  );
}

test.describe("@upload the site as an MCP server", () => {
  test.skip(
    !ENABLED,
    "Set E2E_UPLOAD=1 (docker infra) or E2E_AGENT=1 (any migrated database).",
  );

  test("connect, use the tools as the person, disconnect", async ({
    page,
    baseURL,
  }) => {
    const origin = new URL(baseURL!).origin;
    const resource = `${origin}/mcp`;

    // Unsigned: a 401 that names this resource's metadata.
    const anonymous = await page.request.post("/mcp", {
      headers: { "content-type": "application/json" },
      data: { jsonrpc: "2.0", id: 1, method: "tools/list" },
    });

    expect(anonymous.status()).toBe(401);
    const challenge = anonymous.headers()["www-authenticate"];

    expect(challenge).toContain(
      `resource_metadata="${origin}/.well-known/oauth-protected-resource/mcp"`,
    );
    const metadata = await (
      await page.request.get("/.well-known/oauth-protected-resource/mcp")
    ).json();

    expect(metadata.resource).toBe(resource);
    expect(metadata.authorization_servers).toEqual([origin]);

    // The person, with two projects; the connector may see one.
    await signInDev(page, `e2e-remote-${Date.now()}@example.com`);
    const make = async (title: string) =>
      (
        await (
          await page.request.post("/api/projects", { data: { title } })
        ).json()
      ).id as string;
    const chosen = await make("Connector may see this");

    await make("Connector may not see this");

    const client = await (
      await page.request.post("/api/oauth/register", {
        data: {
          client_name: "E2E Connector",
          redirect_uris: ["http://127.0.0.1/callback"],
        },
      })
    ).json();
    const verifier = randomBytes(32).toString("base64url");
    const authorize = new URL("/oauth/authorize", origin);

    authorize.search = new URLSearchParams({
      response_type: "code",
      client_id: client.client_id,
      redirect_uri: REDIRECT,
      code_challenge: createHash("sha256").update(verifier).digest("base64url"),
      code_challenge_method: "S256",
      scope: "projects:read",
      resource,
      state: "e2e-remote",
    }).toString();
    await page.goto(authorize.toString());
    await page.getByLabel("Only these projects").check();
    await page.getByLabel("Connector may see this").check();
    const [callback] = await Promise.all([
      page.waitForRequest((r) => r.url().startsWith(REDIRECT)),
      page.getByRole("button", { name: "Allow" }).click(),
    ]);
    const token = await (
      await page.request.post("/api/oauth/token", {
        form: {
          grant_type: "authorization_code",
          code: new URL(callback.url()).searchParams.get("code")!,
          redirect_uri: REDIRECT,
          client_id: client.client_id,
          code_verifier: verifier,
          resource,
        },
      })
    ).json();

    expect(token.token_type).toBe("Bearer");

    // A fresh context: no cookie, only the connector's token.
    const agent = await page
      .context()
      .browser()!
      .newContext({ baseURL: origin });
    const call = async (name: string, args: Record<string, unknown> = {}) =>
      (
        await rpc(agent.request, token.access_token, "tools/call", {
          name,
          arguments: args,
        })
      ).result;

    const init = await rpc(agent.request, token.access_token, "initialize", {
      protocolVersion: "2025-06-18",
      capabilities: {},
      clientInfo: { name: "e2e", version: "0" },
    });

    expect(init.result.serverInfo.name).toBe("merfisheyes");
    const tools = await rpc(agent.request, token.access_token, "tools/list");

    // The package's tools minus those that need the person's computer.
    expect(
      tools.result.tools.map((t: { name: string }) => t.name).sort(),
    ).toEqual([
      "meyes_find_datasets",
      "meyes_get_dataset",
      "meyes_get_processing_status",
      "meyes_get_project",
      "meyes_index_values",
      "meyes_list_facets",
      "meyes_list_genes",
      "meyes_list_my_datasets",
      "meyes_list_projects",
      "meyes_open_view",
      "meyes_search_datasets",
      "meyes_whoami",
    ]);
    expect((await call("meyes_whoami")).content[0].text).toContain(
      "Signed in to",
    );
    const projects = await call("meyes_list_projects", {
      response_format: "json",
    });

    expect(
      projects.structuredContent.items.map((p: { id: string }) => p.id),
    ).toEqual([chosen]);
    expect((await call("meyes_list_facets")).isError).toBeFalsy();

    // A /mcp token is for /mcp only.
    const direct = await agent.request.get("/api/agent/v1/projects", {
      headers: { authorization: `Bearer ${token.access_token}` },
    });

    expect(direct.status()).toBe(401);

    // Disconnect it from the account page: the connector is refused at once.
    await page.goto("/account");
    const row = page.getByRole("listitem").filter({ hasText: "E2E Connector" });

    await row.getByRole("button", { name: "Disconnect" }).click();
    await expect(row).toHaveCount(0);
    const after = await agent.request.post("/mcp", {
      headers: {
        authorization: `Bearer ${token.access_token}`,
        "content-type": "application/json",
        accept: "application/json, text/event-stream",
      },
      data: { jsonrpc: "2.0", id: 1, method: "tools/list" },
    });

    expect(after.status()).toBe(401);
    await agent.close();
  });
});
