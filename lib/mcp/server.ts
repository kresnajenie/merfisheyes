import type { Host } from "@/packages/merfisheyes-mcp/src/host.ts";

import { McpServer, createMcpHandler } from "@modelcontextprotocol/server";

import { instructions } from "@/packages/merfisheyes-mcp/src/instructions.ts";
import { registerAccountTools } from "@/packages/merfisheyes-mcp/src/tools/account.ts";
import { registerCatalogTools } from "@/packages/merfisheyes-mcp/src/tools/catalog.ts";
import { registerLibraryTools } from "@/packages/merfisheyes-mcp/src/tools/library.ts";
import { registerViewTools } from "@/packages/merfisheyes-mcp/src/tools/view.ts";

/** The version `/mcp` reports: the package's, whose tools it serves. */
export const MCP_SERVER_VERSION = "0.3.0";

/**
 * The site's MCP server for one request (MCP.md pillar E): the npm
 * package's own tools (D23) — catalog, library search, account, views — minus those that
 * need the person's computer (signing in there, opening a browser, serving
 * folders). `host` is how they reach the agent API.
 */
export function mcpHandler(host: Host) {
  return createMcpHandler(() => {
    const server = new McpServer(
      { name: "merfisheyes", version: MCP_SERVER_VERSION },
      { instructions: instructions(host.site, false) },
    );

    registerCatalogTools(server, host);
    registerLibraryTools(server);
    registerAccountTools(server, host);
    registerViewTools(server, host);

    return server;
  });
}
