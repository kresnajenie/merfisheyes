#!/usr/bin/env node
/**
 * merfisheyes-mcp — an MCP server (stdio) for MERFISHEYES.
 *
 *   merfisheyes view DATA      open data on this computer in the viewer
 *                              (also: merfisheyes-mcp …; the same program)
 *   merfisheyes setup          the one-time converter setup (asks nothing)
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
import { resolve } from "node:path";

import { McpServer } from "@modelcontextprotocol/server";
import { serveStdio } from "@modelcontextprotocol/server/stdio";

import { apiGet } from "./api.ts";
import { openBrowser, signOut, startSignIn } from "./auth.ts";
import { NAME, SITE, VERSION } from "./config.ts";
import { credentialsPath } from "./credentials.ts";
import { SetupNeeded, setupPipeline, viewable } from "./local/convert.ts";
import { local, openLink } from "./local/serve.ts";
import type { Host } from "./host.ts";
import { instructions } from "./instructions.ts";
import { registerAccountTools } from "./tools/account.ts";
import { registerCatalogTools } from "./tools/catalog.ts";
import { registerLocalTools } from "./tools/local.ts";
import { registerSessionTools } from "./tools/session.ts";
import { registerSignInTools } from "./tools/signin.ts";
import { registerViewTools } from "./tools/view.ts";

/** How this package reaches the site and the person (src/host.ts). */
const host: Host = {
  site: SITE,
  get: apiGet,
  open: openBrowser,
  signInHint: `No saved sign-in in ${credentialsPath()}. Public catalog tools work; call meyes_sign_in for the person's own data.`,
};

function createServer(): McpServer {
  const server = new McpServer(
    { name: NAME, version: VERSION },
    { instructions: instructions(SITE, true) },
  );

  registerCatalogTools(server, host);
  registerSignInTools(server);
  registerAccountTools(server, host);
  registerViewTools(server, host);
  registerLocalTools(server);
  registerSessionTools(server);

  return server;
}

/** A yes/no question in a terminal; "no" when there is no terminal. */
async function confirm(question: string): Promise<boolean> {
  if (!process.stdin.isTTY) return false;
  const { createInterface } = await import("node:readline/promises");
  const rl = createInterface({ input: process.stdin, output: process.stderr });
  const answer = (await rl.question(question)).trim().toLowerCase();

  rl.close();

  return answer === "" || answer === "y" || answer === "yes";
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
        console.error(
          "Usage: merfisheyes-mcp view <folder | .h5ad | Xenium/MERSCOPE folder | molecule table> [--no-open]",
        );

        return 2;
      }
      if (!viewable(resolve(path)))
        console.error(
          `Converting ${path} with MERFISHEYES's pipeline (once; reused until it changes)…`,
        );
      const opened = await local
        .open(path, (line) => console.error(`  ${line}`))
        .catch(async (e) => {
          if (!(e instanceof SetupNeeded)) throw e;
          console.error(`\n${e.message}\n`);
          if (!(await confirm("Set it up now? [Y/n] ")))
            throw new Error(
              "Not set up. Run `merfisheyes setup` when you are ready, then `merfisheyes view` again.",
            );
          await setupPipeline((line) => console.error(`  ${line}`));

          return local.open(path, (line) => console.error(`  ${line}`));
        });
      const { share, converted } = opened;
      const link = openLink(share);

      if (converted)
        console.error(
          converted.reused
            ? `Reusing the earlier conversion in ${converted.folder}.`
            : `Converted in ${Math.round(converted.seconds)} s into ${converted.folder}.`,
        );

      console.error(
        [
          `Serving ${share.path} from this computer — nothing is uploaded.`,
          "",
          `  ${link}`,
          "",
          local.viewerDir
            ? "The viewer runs from this computer too: nothing goes to the internet."
            : `In Chrome, the first time: allow ${new URL(SITE).host} to reach apps on this device when it asks.`,
          "Press Ctrl+C to stop.",
        ].join("\n"),
      );
      if (!args.includes("--no-open")) openBrowser(link);

      // Serve until stopped.
      return new Promise<number>((done) =>
        process.once("SIGINT", () => local.close().then(() => done(0))),
      );
    }
    case "setup": {
      const python = await setupPipeline((line) => console.error(`  ${line}`));

      console.error(`MERFISHEYES's converter is ready (Python at ${python}).`);

      return 0;
    }
    default:
      console.error(
        `Unknown command "${command}". Commands: view, setup, login, logout, whoami; none serves MCP on stdio.`,
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
