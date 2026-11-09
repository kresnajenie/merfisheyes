import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { apiGet } from "../api.js";
import { openBrowser, signOut, startSignIn } from "../auth.js";
import { DEFAULT_SCOPES, SITE } from "../config.js";
import {
  type Page,
  join,
  many,
  n,
  one,
  paging,
  responseFormat,
  text,
} from "../format.js";

/**
 * Signing in, and the person's own datasets and projects
 * (MCP-PROTOCOL §3.1). Everything here but sign-in needs the person signed
 * in; what an agent reaches is what they chose when they allowed it.
 */

interface DatasetCard {
  id: string;
  title: string | null;
  kind: string;
  status: string;
  cells: number | null;
  genes: number | null;
  species: string | null;
  tissue: string | null;
  platform: string | null;
  viewerUrl: string | null;
  createdAt: string;
}

interface ProjectRow {
  id: string;
  title: string;
  species: string | null;
  tissue: string | null;
  platform: string | null;
  datasetIds: string[];
  updatedAt: string;
}

const datasetMd = (d: DatasetCard) =>
  [
    `## ${d.title ?? "Untitled"} (${d.id})`,
    `- ${join(d.kind, d.status, d.species, d.tissue, d.platform)}`,
    `- ${n(d.cells)} cells · ${n(d.genes)} genes · added ${d.createdAt.slice(0, 10)}`,
    ...(d.viewerUrl ? [`- Open: ${d.viewerUrl}`] : []),
  ].join("\n");

const SCOPES = [
  "catalog:read",
  "projects:read",
  "projects:write",
  "datasets:write",
  "viewer:control",
  "export",
] as const;

const READ = {
  readOnlyHint: true,
  destructiveHint: false,
  idempotentHint: true,
  openWorldHint: true,
};

export function registerAccountTools(server: McpServer) {
  server.registerTool(
    "meyes_sign_in",
    {
      title: "Sign in to MERFISHEYES",
      description: `Start signing the person in to MERFISHEYES, so tools can reach their own projects and datasets. Returns a link: show it to the person and ask them to open it (a browser tab is also opened when possible), sign in, and choose what to allow — all their projects or chosen ones. Then call meyes_whoami to confirm.
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
    async ({ scopes }) => {
      const signIn = await startSignIn([
        ...new Set([...DEFAULT_SCOPES, ...(scopes ?? [])]),
      ]);

      openBrowser(signIn.url);

      return text(
        `Ask the person to open this link, sign in to MERFISHEYES and choose what to allow:\n\n${signIn.url}\n\nThe link works for 10 minutes and only on this computer. When they are done, call meyes_whoami.`,
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

  server.registerTool(
    "meyes_whoami",
    {
      title: "Who is signed in",
      description:
        "Whether the person is signed in to MERFISHEYES, as whom, what this agent may do (scopes) and which projects it reaches.",
      inputSchema: z.strictObject({ response_format: responseFormat }),
      annotations: READ,
    },
    async ({ response_format }, ctx) => {
      const me = await apiGet<{
        signedIn: boolean;
        name?: string | null;
        email?: string | null;
        scopes: string[];
        projects: "all" | string[];
      }>("/me", {}, ctx.mcpReq.signal);

      return one(response_format, me, () =>
        me.signedIn
          ? [
              `Signed in to ${SITE} as ${me.name ?? me.email} (${me.email}).`,
              `May: ${me.scopes.join(", ")}.`,
              `Projects: ${me.projects === "all" ? "all of theirs" : `${me.projects.length} chosen (${me.projects.join(", ")})`}.`,
            ].join("\n")
          : `Not signed in to ${SITE}. Public catalog tools work; call meyes_sign_in for the person's own data.`,
      );
    },
  );

  server.registerTool(
    "meyes_list_my_datasets",
    {
      title: "List my datasets",
      description: `The signed-in person's datasets that this agent may reach, newest first, with status (QUEUED, PROCESSING, COMPLETE, FAILED) and viewer links.
Examples: "my datasets" -> {}; "which of my uploads failed" -> status="FAILED"; "my hippocampus data" -> query="hippocampus".`,
      inputSchema: z.strictObject({
        query: z.string().max(200).optional().describe("Fragment of the title"),
        status: z
          .enum(["QUEUED", "PROCESSING", "COMPLETE", "FAILED"])
          .optional()
          .describe("Only datasets in this state"),
        ...paging,
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async ({ query, status, limit, offset, response_format }, ctx) => {
      const page = await apiGet<Page<DatasetCard>>(
        "/me/datasets",
        { q: query, status, limit, offset },
        ctx.mcpReq.signal,
      );

      return many(response_format, page, {
        heading: "Your datasets",
        item: datasetMd,
        empty:
          "No dataset of yours matches (or this agent was allowed only some projects: see meyes_whoami).",
      });
    },
  );

  server.registerTool(
    "meyes_get_processing_status",
    {
      title: "Get processing status",
      description:
        "Where one of the person's datasets is in server processing (QUEUED → PROCESSING → COMPLETE or FAILED), its progress, and on failure the reason and whose fault it was. Use after an upload; meyes_list_my_datasets gives the ids.",
      inputSchema: z.strictObject({
        id: z
          .string()
          .min(1)
          .max(100)
          .describe("App dataset id, e.g. ds_a1B2c3D4e5"),
        response_format: responseFormat,
      }),
      annotations: { ...READ, idempotentHint: false },
    },
    async ({ id, response_format }, ctx) => {
      const s = await apiGet<
        DatasetCard & {
          progress: unknown;
          errorMessage: string | null;
          faultCategory: string | null;
        }
      >(`/me/datasets/${encodeURIComponent(id)}/status`, {}, ctx.mcpReq.signal);

      return one(response_format, s, () =>
        [
          `${s.title ?? s.id} (${s.id}): ${s.status}`,
          ...(s.progress ? [`Progress: ${JSON.stringify(s.progress)}`] : []),
          ...(s.errorMessage
            ? [
                `Failed: ${s.errorMessage}${s.faultCategory ? ` (fault: ${s.faultCategory})` : ""}`,
              ]
            : []),
          ...(s.viewerUrl ? [`Open: ${s.viewerUrl}`] : []),
        ].join("\n"),
      );
    },
  );

  server.registerTool(
    "meyes_list_projects",
    {
      title: "List my projects",
      description:
        "The signed-in person's projects that this agent may reach, most recently changed first, each with its dataset ids. Use meyes_get_project for one with its datasets in full.",
      inputSchema: z.strictObject({
        query: z.string().max(200).optional().describe("Fragment of the title"),
        ...paging,
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async ({ query, limit, offset, response_format }, ctx) => {
      const page = await apiGet<Page<ProjectRow>>(
        "/projects",
        { q: query, limit, offset },
        ctx.mcpReq.signal,
      );

      return many(response_format, page, {
        heading: "Your projects",
        item: (p) =>
          [
            `## ${p.title} (${p.id})`,
            `- ${join(p.species, p.tissue, p.platform) || "no metadata"}`,
            `- Datasets: ${p.datasetIds.join(", ") || "none"}`,
          ].join("\n"),
        empty:
          "No project of yours matches (or this agent was allowed only some projects: see meyes_whoami).",
      });
    },
  );

  server.registerTool(
    "meyes_get_project",
    {
      title: "Get a project",
      description:
        "One of the person's projects with its datasets (status, counts, viewer links). Ids come from meyes_list_projects.",
      inputSchema: z.strictObject({
        id: z.string().min(1).max(100).describe("Project id"),
        response_format: responseFormat,
      }),
      annotations: READ,
    },
    async ({ id, response_format }, ctx) => {
      const p = await apiGet<
        ProjectRow & { description: string | null; datasets: DatasetCard[] }
      >(`/projects/${encodeURIComponent(id)}`, {}, ctx.mcpReq.signal);

      return one(response_format, p, () =>
        [
          `# ${p.title} (${p.id})`,
          p.description ?? "",
          ...p.datasets.map(datasetMd),
        ]
          .filter(Boolean)
          .join("\n\n"),
      );
    },
  );
}
