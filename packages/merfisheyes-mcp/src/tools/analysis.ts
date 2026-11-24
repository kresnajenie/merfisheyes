import type { McpServer } from "@modelcontextprotocol/server";

import { existsSync, readFileSync } from "node:fs";
import { basename, dirname, join, resolve } from "node:path";

import * as z from "zod/v4";

import { text } from "../format.ts";
import {
  type Job,
  type Step,
  jobById,
  startStep,
  waitFor,
  workspaceFor,
} from "../local/analysis.ts";
import { SetupNeeded } from "../local/convert.ts";
import { local } from "../local/serve.ts";

/**
 * Analysis on the dataset open in the view (MCP.md pillar D, tier 2): the
 * steps an agent would otherwise script by hand — profile, a niche around
 * anchor cells, a re-clustering with its UMAP, markers, a saved copy — run
 * in MERFISHEYES's own Python as jobs (src/local/analysis.ts). New labels
 * and embeddings appear in the open view at once.
 */

const ACT = {
  readOnlyHint: false,
  destructiveHint: false,
  idempotentHint: false,
  openWorldHint: false,
};

const path = z
  .string()
  .max(4096)
  .optional()
  .describe(
    "The dataset (as given to meyes_open_local); default: the one opened last",
  );

const cellsSpec = z.strictObject({
  column: z.string().min(1).max(200),
  values: z.array(z.string().max(200)).min(1).max(500),
});

/** Compact JSON for the agent, cut to a size a client takes. */
const compact = (v: unknown, max = 20_000) => {
  const s = JSON.stringify(v);

  return s.length > max
    ? `${s.slice(0, max)}\n… (cut; ask for fewer values)`
    : s;
};

function report(job: Job, extra: string[] = []) {
  const secs = ((Date.now() - job.started) / 1000).toFixed(0);

  if (job.status === "running")
    return text(
      [
        `Still running after ${secs} s (job ${job.id}). Last: ${job.lines.at(-1) ?? "starting"}.`,
        `Call meyes_job with id="${job.id}" to wait for it; tell the person it is running.`,
      ].join("\n"),
    );
  if (job.status === "failed")
    return {
      isError: true,
      content: [
        { type: "text" as const, text: job.error ?? "The step failed." },
      ],
    };
  const result: Record<string, unknown> = { ...job.result };

  delete result.ok;
  const image =
    job.step === "markers" &&
    typeof result.heatmap === "string" &&
    existsSync(result.heatmap)
      ? [
          {
            type: "image" as const,
            data: readFileSync(result.heatmap).toString("base64"),
            mimeType: "image/png",
          },
        ]
      : [];

  return {
    content: [
      ...image,
      {
        type: "text" as const,
        text: [
          compact(result),
          ...(job.shown?.length
            ? [`In the open view now: ${job.shown.join(", ")}.`]
            : []),
          ...extra,
        ].join("\n"),
      },
    ],
  };
}

async function runStep(
  step: Step,
  args: Record<string, unknown>,
  dataset: string | undefined,
  extra: string[] = [],
) {
  try {
    const d = local.dataset(dataset);
    const job = await startStep(step, args, d);

    return report(await waitFor(job), extra);
  } catch (e) {
    if (e instanceof SetupNeeded)
      return text(
        [
          `Not run: ${e.message}`,
          "",
          "Ask the person whether that download is OK. Only if they agree, call meyes_setup_converter, then try again.",
        ].join("\n"),
      );

    return {
      isError: true,
      content: [{ type: "text" as const, text: (e as Error).message }],
    };
  }
}

export function registerAnalysisTools(server: McpServer) {
  server.registerTool(
    "meyes_profile",
    {
      title: "Profile the open dataset",
      description: `What the dataset open in the view holds, before deciding anything: cells and genes, every obs column (categorical: values with cell counts; numeric: range), embeddings, the cells' spacing (whether coordinates are in microns), which columns already have DE stats, results already in the workspace.
genes: check which marker genes the panel measures and in what fraction of cells (panels are often a few hundred genes: say which markers are missing). breakdown: counts of a column's values per section/sample/batch (e.g. which section holds the dentate gyrus cells).
Runs in MERFISHEYES's own Python (set up once, with the person's yes).`,
      inputSchema: z.strictObject({
        path,
        genes: z
          .array(z.string().max(100))
          .max(300)
          .optional()
          .describe("Marker genes to check against the panel"),
        breakdown: z
          .strictObject({
            column: z.string().min(1).max(200),
            by: z
              .string()
              .min(1)
              .max(200)
              .describe("e.g. the section or sample column"),
            values: z.array(z.string().max(200)).max(200).optional(),
          })
          .optional(),
        max_values: z.number().int().min(1).max(200).default(25),
      }),
      annotations: { ...ACT, readOnlyHint: true, idempotentHint: true },
    },
    async ({ path: p, ...a }) => runStep("profile", a, p),
  );

  server.registerTool(
    "meyes_niche",
    {
      title: "Select a niche around anchor cells",
      description: `Mark the cells within radius_um of anchor cells (cells with given labels) — a niche such as the dentate gyrus around its granule and immature neurons — as a new label column niche_<name> ("<name>" inside, "" outside). Shown in the open view at once (colour by it; hide "" to see only the niche).
within: a section/sample column, so distances are measured inside each section only. exclude: cells never in the niche (e.g. olfactory-bulb neuroblasts when the niche is the SVZ). min_anchors > 1 drops isolated anchor cells. describe_by: the niche's composition by a label column.
Check the coordinate units first (meyes_profile): radius_um assumes microns. Ask the person which anchors and radius they want; say what you chose and why.`,
      inputSchema: z.strictObject({
        path,
        name: z
          .string()
          .regex(/^[A-Za-z0-9_]+$/)
          .max(40)
          .describe('e.g. "DG", "SVZ"'),
        anchors: cellsSpec,
        radius_um: z.number().positive().max(5000),
        within: z.string().max(200).optional(),
        exclude: cellsSpec.optional(),
        min_anchors: z.number().int().min(1).max(100).default(1),
        describe_by: z.string().max(200).optional(),
      }),
      annotations: ACT,
    },
    async ({ path: p, ...a }) => runStep("niche", a, p),
  );

  server.registerTool(
    "meyes_recluster",
    {
      title: "Re-cluster a subset with its own UMAP",
      description: `For a subset of cells (a niche, some cell types): normalise (log1p of normalize_total), PCA, neighbours, UMAP and Leiden at each resolution. Writes label columns <name>_leiden_<res> (C0, C1, …; "" outside the subset) and an embedding X_umap_<name> (shown as umap_<name>), all shown in the open view at once. For the primary resolution it answers each cluster's size, top genes and — with label_column — its main existing labels, and with within its cells per section, to name the clusters.
Name clusters from that evidence (and the literature), mark uncertain names with "?", and say a cluster made of one section only may be a batch effect. Long subsets run as a job (meyes_job).`,
      inputSchema: z.strictObject({
        path,
        name: z
          .string()
          .regex(/^[A-Za-z0-9_]+$/)
          .max(40),
        cells: cellsSpec.describe(
          'The subset, e.g. {column:"niche_DG", values:["DG"]}',
        ),
        resolutions: z
          .array(z.number().positive().max(5))
          .min(1)
          .max(8)
          .default([0.5]),
        primary_resolution: z.number().positive().max(5).optional(),
        label_column: z.string().max(200).optional(),
        within: z.string().max(200).optional(),
        n_pcs: z.number().int().min(2).max(100).default(30),
        n_neighbors: z.number().int().min(2).max(200).default(15),
        top_n: z.number().int().min(1).max(30).default(8),
      }),
      annotations: ACT,
    },
    async ({ path: p, ...a }) => runStep("recluster", a, p),
  );

  server.registerTool(
    "meyes_markers",
    {
      title: "Marker genes per label",
      description: `Top genes for each value of a label column (cell types, clusters from meyes_recluster): mean log-normalised expression in the value against the other cells, and the fraction of its cells expressing. Uses the dataset's own DE stats when it has them for that column (instant), else computes them the same way.
genes + heatmap=true: a heatmap of chosen genes (e.g. a lineage's markers from the literature) across the values, returned as an image and saved — names the genes the panel lacks. These are markers, not proof of what drives a process: say so.`,
      inputSchema: z.strictObject({
        path,
        column: z.string().min(1).max(200),
        values: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe("Only these values (default: all)"),
        cells: cellsSpec.optional().describe("Compute within this subset only"),
        top_n: z.number().int().min(1).max(50).default(10),
        genes: z.array(z.string().max(100)).max(100).optional(),
        heatmap: z.boolean().default(false),
        scale: z.enum(["zscore", "minmax"]).default("zscore"),
        title: z.string().max(200).optional(),
        save_to: z
          .string()
          .max(4096)
          .optional()
          .describe("Where to save the heatmap PNG"),
      }),
      annotations: { ...ACT, readOnlyHint: true },
    },
    async ({ path: p, heatmap, save_to, ...a }) => {
      let file: string | undefined;

      if (heatmap) {
        try {
          file = save_to
            ? resolve(save_to)
            : join(
                workspaceFor(local.dataset(p)),
                "figures",
                `markers-${a.column}-${Date.now()}.png`,
              );
        } catch (e) {
          return {
            isError: true,
            content: [{ type: "text" as const, text: (e as Error).message }],
          };
        }
        const { mkdirSync } = await import("node:fs");

        mkdirSync(dirname(file), { recursive: true });
      }

      return runStep(
        "markers",
        { ...a, ...(file ? { heatmap: file } : {}) },
        p,
      );
    },
  );

  server.registerTool(
    "meyes_save_analysis",
    {
      title: "Save the results into a copy",
      description: `Write the workspace's results (labels from meyes_niche / meyes_recluster, embeddings X_umap_<name>) into a copy of the dataset — never into the original. Ask the person which format first: "zarr" (a SpatialData store, the format MERFISHEYES reads best: opens at once, and the new labels carry DE stats for the viewer's DEG panel) or "h5ad" (their original .h5ad with the new obs columns and obsm embeddings, for scanpy).
Default place: beside the original, <name>_merfisheyes.zarr / .h5ad. open=true (zarr) switches the open view to the copy.`,
      inputSchema: z.strictObject({
        path,
        format: z.enum(["zarr", "h5ad"]),
        dest: z.string().max(4096).optional(),
        columns: z
          .array(z.string().max(200))
          .max(100)
          .optional()
          .describe("Default: all the workspace's labels"),
        embeddings: z
          .array(z.string().max(200))
          .max(20)
          .optional()
          .describe("Default: all the workspace's embeddings"),
        overwrite: z.boolean().default(false),
        open: z.boolean().default(false),
      }),
      annotations: ACT,
    },
    async ({ path: p, format, dest, open, ...a }) => {
      let d;

      try {
        d = local.dataset(p);
      } catch (e) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: (e as Error).message }],
        };
      }
      const base = basename(d.input).replace(/\.(h5ad|zarr)$/i, "");
      const target = resolve(
        dest ?? join(dirname(d.input), `${base}_merfisheyes.${format}`),
      );
      const extra: string[] = [];
      const res = await runStep(
        "save",
        { ...a, format, dest: target, original: d.input },
        p,
      );

      if (
        !("isError" in res) &&
        open &&
        format === "zarr" &&
        existsSync(target)
      ) {
        try {
          const { share } = await local.open(target);
          const link = await local.navigate(share);

          extra.push(`The open view now shows the copy (${link}).`);
        } catch (e) {
          extra.push(
            `Saved; open it with meyes_open_local (${(e as Error).message}).`,
          );
        }
        res.content.push({ type: "text" as const, text: extra.join("\n") });
      }

      return res;
    },
  );

  server.registerTool(
    "meyes_job",
    {
      title: "Wait for an analysis job",
      description:
        "The result of an analysis step that was still running (meyes_profile, meyes_niche, meyes_recluster, meyes_markers, meyes_save_analysis hand back a job id when they take longer than a tool call): waits up to wait_seconds more, then answers with the result or that it is still running.",
      inputSchema: z.strictObject({
        id: z.string().min(1).max(100),
        wait_seconds: z.number().int().min(0).max(50).default(45),
      }),
      annotations: { ...ACT, readOnlyHint: true, idempotentHint: true },
    },
    async ({ id, wait_seconds }) => {
      try {
        return report(await waitFor(jobById(id), wait_seconds * 1000));
      } catch (e) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: (e as Error).message }],
        };
      }
    },
  );
}
