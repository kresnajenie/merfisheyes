import { describe, expect, it } from "vitest";

import { mcpHandler } from "@/lib/mcp/server";
import type { Host } from "@/packages/merfisheyes-mcp/src/host.ts";

/**
 * The site's MCP server (`/mcp`, lib/mcp/server.ts) over HTTP, with a
 * stand-in for the agent API: what a connector sees — the tools (the
 * package's, minus the local ones), the connector's instructions, links
 * only (no "open"), and the API's answers turned into tool results. Signing
 * in and the in-process API are covered by e2e/tests/upload/agent-remote-mcp.spec.ts.
 */

const SITE = "https://merfisheyes.test";
const calls: string[] = [];
const host: Host = {
  site: SITE,
  async get<T>(path: string, params: Record<string, string | number | undefined>) {
    calls.push(`${path}?${new URLSearchParams(params as Record<string, string>)}`);
    if (path === "/me")
      return {
        signedIn: true,
        name: "Ada",
        email: "ada@example.com",
        scopes: ["catalog:read", "projects:read"],
        projects: "all",
      } as T;
    if (path === "/catalog/facets")
      return { species: ["Mouse"], tissues: ["Brain"], platforms: ["MERFISH"] } as T;
    throw new Error("Not found. Check the id.");
  },
};

async function rpc(method: string, params: unknown = {}) {
  const res = await mcpHandler(host).fetch(
    new Request(`${SITE}/mcp`, {
      method: "POST",
      headers: {
        "content-type": "application/json",
        accept: "application/json, text/event-stream",
      },
      body: JSON.stringify({ jsonrpc: "2.0", id: 1, method, params }),
    }),
  );
  const text = await res.text();
  // JSON, or one SSE event carrying it
  const json = res.headers.get("content-type")?.includes("event-stream")
    ? JSON.parse(/^data: (.*)$/m.exec(text)![1])
    : JSON.parse(text);

  return json as { result?: any; error?: any };
}

describe("the site's MCP server", () => {
  it("introduces itself with the connector's instructions", async () => {
    const r = await rpc("initialize", {
      protocolVersion: "2025-06-18",
      capabilities: {},
      clientInfo: { name: "test", version: "0" },
    });

    expect(r.result.serverInfo.name).toBe("merfisheyes");
    expect(r.result.instructions).toContain(SITE);
    expect(r.result.instructions).toContain("added this connector");
    expect(r.result.instructions).not.toContain("meyes_sign_in");
  });

  it("offers the package's tools except the ones that need the person's computer", async () => {
    const r = await rpc("tools/list");
    const names = r.result.tools.map((t: { name: string }) => t.name).sort();

    expect(names).toEqual([
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
    const view = r.result.tools.find((t: { name: string }) => t.name === "meyes_open_view");

    expect(Object.keys(view.inputSchema.properties)).not.toContain("open");
  });

  it("answers from the API, and turns its failures into tool errors", async () => {
    const me = await rpc("tools/call", { name: "meyes_whoami", arguments: {} });

    expect(me.result.content[0].text).toContain(`Signed in to ${SITE} as Ada`);
    const facets = await rpc("tools/call", { name: "meyes_list_facets", arguments: {} });

    expect(facets.result.content[0].text).toContain("MERFISH");
    const missing = await rpc("tools/call", {
      name: "meyes_get_dataset",
      arguments: { id: "nope" },
    });

    expect(missing.result.isError).toBe(true);
    expect(missing.result.content[0].text).toContain("Not found");
    expect(calls).toContain("/catalog/nope?");
  });
});
