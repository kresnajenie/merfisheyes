#!/usr/bin/env node
/**
 * merfisheyes-mcp — an MCP server (stdio) for MERFISHEYES.
 *
 *   merfisheyes-mcp            serve MCP on stdio (what an agent runs)
 *   merfisheyes-mcp login      sign in from a terminal, then exit
 *   merfisheyes-mcp logout     sign out
 *   merfisheyes-mcp whoami     who is signed in
 *   merfisheyes-mcp view DIR   open a folder on this computer in the viewer
 *
 * MERFISHEYES_URL points it at another deployment (default the live site).
 * Nothing is ever written to stdout but the protocol.
 */
/* eslint-disable no-console -- stderr is this server's only log: stdout carries the protocol */
import { McpServer } from "@modelcontextprotocol/server";
import { serveStdio } from "@modelcontextprotocol/server/stdio";

import { apiGet } from "./api.js";
import { openBrowser, signOut, startSignIn } from "./auth.js";
import { NAME, SITE, VERSION } from "./config.js";
import { credentialsPath } from "./credentials.js";
import { local, openLink } from "./local/serve.js";
import { registerAccountTools } from "./tools/account.js";
import { registerCatalogTools } from "./tools/catalog.js";
import { registerLocalTools } from "./tools/local.js";
import { registerViewTools } from "./tools/view.js";

const INSTRUCTIONS = `MERFISHEYES (${SITE}) is a web viewer for spatial transcriptomics: cells, molecules, cell outlines and microscopy images, at up to billions of molecules. These tools find datasets and read the person's own projects; the site does the heavy work.

- Public catalog tools (meyes_search_datasets, meyes_list_facets, meyes_get_dataset, meyes_list_genes) need no sign-in. The person's own data needs meyes_sign_in once; ask only for the scopes the task needs.
- Prefer what MERFISHEYES already records over computing it yourself: meyes_get_dataset returns the counts, cell-type levels, DE columns and QC each store's own files hold. Never fetch raw per-cell data to answer something a tool answers.
- To show the person something, call meyes_open_view: it returns a link to the viewer already set up (genes, column, cell types, 2D/3D, dot size, opacity). Give them the link. Pass open=true only when they ask to open it now, and never open more than one tab per request: for several views give links, and for two views to compare use one split screen (compare=…). Do not describe a scene you cannot see.
- Data on the person's own computer (a .zarr store, a MERFISHEYES folder, an OME-Zarr image) opens with meyes_open_local: it is served from their machine to the viewer, never uploaded. Use it, not an upload, when they ask to see their local data.
- When the person asks to sign in, or asks for their own projects or datasets while signed out, call meyes_sign_in: it opens the sign-in page and waits for them.
- Lists are paged: use limit/offset and next_offset rather than repeating a call. Calls are rate-limited; a 429 says how long to wait.
- Not every column with DE stats is a cell type (some stores carry QC columns too); group by annotation columns.`;

function createServer(): McpServer {
  const server = new McpServer(
    { name: NAME, version: VERSION },
    { instructions: INSTRUCTIONS },
  );

  registerCatalogTools(server);
  registerAccountTools(server);
  registerViewTools(server);
  registerLocalTools(server);

  return server;
}

async function cli(command: string): Promise<number> {
  switch (command) {
    case "login": {
      const s = await startSignIn();

      console.error(`Open this link to sign in to ${SITE}:\n\n  ${s.url}\n`);
      openBrowser(s.url);
      await s.done;
      console.error(`Signed in. Credentials kept in ${credentialsPath()}`);

      return 0;
    }
    case "logout":
      console.error((await signOut()) ? "Signed out." : "Was not signed in.");

      return 0;
    case "whoami": {
      const me = await apiGet<{
        signedIn: boolean;
        email?: string;
        scopes: string[];
      }>("/me", {});

      console.error(
        me.signedIn
          ? `${me.email} on ${SITE} — ${me.scopes.join(", ")}`
          : `Not signed in to ${SITE} (no saved sign-in in ${credentialsPath()}).`,
      );

      return 0;
    }
    case "view": {
      const args = process.argv.slice(3);
      const path = args.find((a) => !a.startsWith("--"));

      if (!path) {
        console.error("Usage: merfisheyes-mcp view <folder> [--no-open]");

        return 2;
      }
      const share = await local.share(path);
      const link = openLink(share);

      console.error(
        [
          `Serving ${share.path} from this computer — nothing is uploaded.`,
          "",
          `  ${link}`,
          "",
          `In Chrome, the first time: allow ${new URL(SITE).host} to reach apps on this device when it asks.`,
          "Press Ctrl+C to stop.",
        ].join("\n"),
      );
      if (!args.includes("--no-open")) openBrowser(link);

      // Serve until stopped.
      return new Promise<number>((done) =>
        process.once("SIGINT", () => local.close().then(() => done(0))),
      );
    }
    default:
      console.error(
        `Unknown command "${command}". Commands: login, logout, whoami, view; none serves MCP on stdio.`,
      );

      return 2;
  }
}

const command = process.argv[2];

if (command) {
  cli(command).then(
    (code) => process.exit(code),
    (e: Error) => {
      console.error(e.message);
      process.exit(1);
    },
  );
} else {
  serveStdio(createServer);
  console.error(`${NAME} ${VERSION} on stdio, for ${SITE}`);
}
