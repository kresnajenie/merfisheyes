import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { apiGet } from "../api.ts";
import { openBrowser, signOut, startSignIn } from "../auth.ts";
import { DEFAULT_SCOPES, SIGN_IN_WAIT_MS, SITE } from "../config.ts";
import { text } from "../format.ts";

/**
 * Signing in and out from the person's computer (the npm package only): a
 * loopback sign-in in their browser, tokens kept in their credentials file.
 * The site's own MCP server has no such tools — a connector signs in when
 * it is added.
 */
const SCOPES = [
  "catalog:read",
  "projects:read",
  "projects:write",
  "datasets:write",
  "viewer:control",
  "export",
] as const;

export function registerSignInTools(server: McpServer) {
  server.registerTool(
    "meyes_sign_in",
    {
      title: "Sign in to MERFISHEYES",
      description: `Start signing the person in to MERFISHEYES, so tools can reach their own projects and datasets. Opens a MERFISHEYES sign-in page in the person's browser, where they sign in and choose what to allow — all their projects or chosen ones — and waits up to two minutes for them, then says who is signed in. If they take longer, or no browser could be opened, it returns the link to give them; call meyes_whoami once they are done.
Call it when a tool says sign-in is needed, or with more scopes when a tool says its permission is missing. Public catalog tools work without it.`,
      inputSchema: z.strictObject({
        scopes: z
          .array(z.enum(SCOPES))
          .optional()
          .describe(
            `What to ask for; default ${DEFAULT_SCOPES.join(", ")}. Ask only for what the task needs.`,
          ),
      }),
      annotations: {
        readOnlyHint: false,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: true,
      },
    },
    async ({ scopes }, ctx) => {
      const asked = [...new Set([...DEFAULT_SCOPES, ...(scopes ?? [])])];
      const signIn = await startSignIn(asked);
      const link = [
        signIn.url,
        ``,
        `It asks for: ${asked.join(", ")} (${DEFAULT_SCOPES.join(" and ")} are always included).`,
      ];

      // No browser opened here: hand the link over at once.
      if (!openBrowser(signIn.url))
        return text(
          [
            `Ask the person to open this link, sign in to MERFISHEYES and choose what to allow:`,
            ``,
            ...link,
            `The link works for 10 minutes and only on this computer. When they are done, call meyes_whoami.`,
          ].join("\n"),
        );

      // A sign-in page is open in their browser: wait for their answer, so
      // the agent can say who is signed in without another round.
      const outcome = await Promise.race([
        signIn.done.then(
          () => "done" as const,
          (e: Error) => e,
        ),
        new Promise<"waiting">((r) => {
          const t = setTimeout(() => r("waiting"), SIGN_IN_WAIT_MS);

          ctx.mcpReq.signal.addEventListener("abort", () => {
            clearTimeout(t);
            r("waiting");
          });
        }),
      ]);

      if (outcome === "done") {
        const me = await apiGet<{
          signedIn: boolean;
          name?: string | null;
          email?: string | null;
          scopes: string[];
          projects: "all" | string[];
        }>("/me", {}, ctx.mcpReq.signal);

        return text(
          me.signedIn
            ? `Signed in to ${SITE} as ${me.name ?? me.email}. This agent may: ${me.scopes.join(", ")}; projects: ${me.projects === "all" ? "all of theirs" : `${me.projects.length} chosen`}.`
            : `The sign-in finished but ${SITE} does not recognise it. Call meyes_sign_in again.`,
        );
      }
      if (outcome instanceof Error)
        return {
          isError: true,
          content: [
            {
              type: "text" as const,
              text: `Sign-in did not complete: ${outcome.message}. Call meyes_sign_in to try again.`,
            },
          ],
        };

      return text(
        [
          `A MERFISHEYES sign-in page is open in the person's browser and has not been answered yet. If they cannot see it, give them this link:`,
          ``,
          ...link,
          `It works for 10 minutes, only on this computer. When they are done, call meyes_whoami.`,
        ].join("\n"),
      );
    },
  );

  server.registerTool(
    "meyes_sign_out",
    {
      title: "Sign out of MERFISHEYES",
      description:
        "Disconnect this agent from the person's MERFISHEYES account: its tokens are revoked and forgotten. Public catalog tools keep working.",
      inputSchema: z.strictObject({}),
      annotations: {
        readOnlyHint: false,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: true,
      },
    },
    async () =>
      text(
        (await signOut()) ? "Signed out of MERFISHEYES." : "Was not signed in.",
      ),
  );
}
