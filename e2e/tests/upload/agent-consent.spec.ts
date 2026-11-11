import { createHash, randomBytes } from "node:crypto";

import { expect, test } from "@playwright/test";

import { signInDev } from "../../helpers/auth";

/**
 * Connecting an agent, in a browser (docs/specs/MCP-PROTOCOL.md §6; step 3
 * of docs/plans/MCP.md): an agent registers, sends the person to
 * /oauth/authorize, the consent page says who is asking and where the
 * answer goes, the person limits it to one project and allows it; the agent
 * swaps the code for a token that reaches that project only; disconnecting
 * it on /account stops the token at once.
 *
 * Tagged @upload: it needs the suite's database (no S3). E2E_AGENT=1 runs it
 * against any migrated database the dev server's DATABASE_URL points at.
 */

const ENABLED = process.env.E2E_UPLOAD === "1" || process.env.E2E_AGENT === "1";
const REDIRECT = "http://127.0.0.1:5999/callback";

test.describe("@upload agent consent", () => {
  test.skip(
    !ENABLED,
    "Set E2E_UPLOAD=1 (docker infra) or E2E_AGENT=1 (any migrated database).",
  );

  test("allow one project, use the token, disconnect", async ({
    page,
    baseURL,
  }) => {
    const origin = new URL(baseURL!).origin;
    const resource = `${origin}/api/agent/v1`;

    await signInDev(page, `e2e-agent-${Date.now()}@example.com`);
    const make = async (title: string) => {
      const res = await page.request.post("/api/projects", { data: { title } });

      expect(res.ok(), `create project ${title}`).toBeTruthy();

      return (await res.json()).id as string;
    };
    const chosen = await make("Agent may see this");
    const other = await make("Agent may not see this");

    // The agent registers itself (RFC 7591) and starts the code flow.
    const client = await (
      await page.request.post("/api/oauth/register", {
        data: {
          client_name: "E2E Agent",
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
      state: "e2e-state",
    }).toString();

    await page.goto(authorize.toString());
    await expect(
      page.getByRole("heading", { name: /Allow E2E Agent/ }),
    ).toBeVisible();
    await expect(page.getByText("127.0.0.1:5999")).toBeVisible();
    await expect(page.getByText(/app running on your computer/)).toBeVisible();
    await expect(page.getByText(/registered itself/)).toBeVisible();
    await expect(
      page.getByText("See your projects and datasets"),
    ).toBeVisible();

    await page.getByLabel("Only these projects").check();
    await page.getByLabel("Agent may see this").check();
    // The answer is a redirect to the agent's loopback listener (none runs
    // here, so the browser's request to it is what we read).
    const [callback] = await Promise.all([
      page.waitForRequest((r) => r.url().startsWith(REDIRECT)),
      page.getByRole("button", { name: "Allow" }).click(),
    ]);
    const back = new URL(callback.url());

    expect(back.searchParams.get("state")).toBe("e2e-state");
    expect(back.searchParams.get("iss")).toBe(origin);

    // The agent swaps the code for tokens.
    const token = await (
      await page.request.post("/api/oauth/token", {
        form: {
          grant_type: "authorization_code",
          code: back.searchParams.get("code")!,
          redirect_uri: REDIRECT,
          client_id: client.client_id,
          code_verifier: verifier,
          resource,
        },
      })
    ).json();

    expect(token.token_type).toBe("Bearer");

    // A fresh context: no cookie, only the bearer token.
    const agent = await page
      .context()
      .browser()!
      .newContext({ baseURL: origin });
    const asAgent = (path: string) =>
      agent.request.get(path, {
        headers: { authorization: `Bearer ${token.access_token}` },
      });
    const projects = await (await asAgent("/api/agent/v1/projects")).json();

    expect(projects.items.map((p: { id: string }) => p.id)).toEqual([chosen]);
    expect((await asAgent(`/api/agent/v1/projects/${other}`)).status()).toBe(
      404,
    );

    // Disconnect it from the account page.
    await page.goto("/account");
    const row = page.getByRole("listitem").filter({ hasText: "E2E Agent" });

    await expect(row).toContainText("Agent may see this");
    await row.getByRole("button", { name: "Disconnect" }).click();
    await expect(row).toHaveCount(0);
    expect((await asAgent("/api/agent/v1/projects")).status()).toBe(401);
    await agent.close();
  });
});
