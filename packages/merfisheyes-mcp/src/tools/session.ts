import type { McpServer } from "@modelcontextprotocol/server";

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, join, resolve } from "node:path";

import * as z from "zod/v4";

import { text } from "../format.ts";
import { cacheRoot } from "../local/convert.ts";
import { local } from "../local/serve.ts";
import { idColumn, parseTable } from "../local/table.ts";

/**
 * The live session (MCP.md pillar D3): the agent reads, changes and
 * pictures the view open in the person's browser — the tab
 * meyes_open_local's link opened, served by this package. The page runs the
 * commands (lib/agent/view-commands.ts in the viewer).
 */

const fail = (e: unknown) => ({
  isError: true,
  content: [{ type: "text" as const, text: (e as Error).message }],
});

const ACT = {
  readOnlyHint: false,
  destructiveHint: false,
  idempotentHint: true,
  openWorldHint: false,
};

export function registerSessionTools(server: McpServer) {
  server.registerTool(
    "meyes_view_state",
    {
      title: "What the open view shows",
      description: `What the MERFISHEYES view open in the person's browser shows right now: the dataset, cells and genes, what it is coloured by (a gene or a column), which values are hidden or picked out, dot size, 2D/3D, and every column with its type and (for categorical ones) its values.
Use it before meyes_view_set to learn the exact column and value names. Needs a view opened with meyes_open_local and still open in a tab.`,
      inputSchema: z.strictObject({}),
      annotations: { ...ACT, readOnlyHint: true },
    },
    async () => {
      try {
        return text(JSON.stringify(await local.command("state", {}), null, 2));
      } catch (e) {
        return fail(e);
      }
    },
  );

  server.registerTool(
    "meyes_view_set",
    {
      title: "Change the open view",
      description: `Change the MERFISHEYES view open in the person's browser: colour by a gene or by a column (cell type, cluster, pseudotime, fov, …), show only some values of the current categorical column or hide some (e.g. a bad FOV), show all again, dot size, 2D or 3D, reset the camera. Answers with the view's state afterwards.
Names must be the dataset's own (meyes_view_state lists them); gene names are matched without regard to case.
Examples: "colour by Slc17a7" -> gene="Slc17a7"; "colour by cell type and hide fov 5" -> column="fov", hide=["5"] then column="celltype"; "only radial glia and immature neurons" -> column="celltype", show_only=["radial glia","immature"].`,
      inputSchema: z.strictObject({
        gene: z
          .string()
          .min(1)
          .max(100)
          .optional()
          .describe("Colour by this gene"),
        column: z
          .string()
          .min(1)
          .max(200)
          .optional()
          .describe(
            "Colour by this column (and the column show_only/hide refer to)",
          ),
        show_only: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe("Show only these values of the column; hide the rest"),
        hide: z
          .array(z.string().max(200))
          .max(200)
          .optional()
          .describe("Hide these values of the column"),
        show_all: z.boolean().optional().describe("Show every value again"),
        treat_as: z
          .enum(["categorical", "numerical"])
          .optional()
          .describe(
            "Treat the column as categories or as numbers (the viewer's switch); a numeric column is treated as categories automatically when you show or hide its values",
          ),
        dot_size: z
          .number()
          .min(0.1)
          .max(10)
          .optional()
          .describe("Dot size multiplier: 1 is the default"),
        view: z.enum(["2D", "3D"]).optional(),
        reset_camera: z
          .boolean()
          .optional()
          .describe("Fit the data in view again"),
      }),
      annotations: ACT,
    },
    async (a) => {
      try {
        return text(JSON.stringify(await local.command("set", a), null, 2));
      } catch (e) {
        return fail(e);
      }
    },
  );

  server.registerTool(
    "meyes_view_snapshot",
    {
      title: "Picture of the open view",
      description: `A PNG of the MERFISHEYES view open in the person's browser, as drawn, with a legend beside it (the visible categories and their colours, or the gene and its colour scale). Returned as an image you can look at, and saved to a file — save_to (a .png path, e.g. in the folder of a slide deck being made) or MERFISHEYES's snapshots folder. Use it to check a view before showing it, and for figures.`,
      inputSchema: z.strictObject({
        save_to: z
          .string()
          .max(4096)
          .optional()
          .describe(
            "Where to save the PNG; default: MERFISHEYES's snapshots folder",
          ),
        width: z
          .number()
          .int()
          .min(400)
          .max(3000)
          .optional()
          .describe(
            "Width of the scene in pixels (the legend is added beside it)",
          ),
      }),
      annotations: ACT,
    },
    async ({ save_to, width }) => {
      try {
        const shot = await local.command<{
          png_base64: string;
          width: number;
          height: number;
          legend?: { value: string; color: string }[];
          gene?: string;
          column?: string;
        }>("snapshot", { width }, 60_000);
        const file = save_to
          ? resolve(save_to)
          : join(
              dirname(cacheRoot()),
              "snapshots",
              `snapshot-${new Date().toISOString().replace(/[:.]/g, "-")}.png`,
            );

        mkdirSync(dirname(file), { recursive: true });
        writeFileSync(file, Buffer.from(shot.png_base64, "base64"));

        return {
          content: [
            {
              type: "image" as const,
              data: shot.png_base64,
              mimeType: "image/png",
            },
            {
              type: "text" as const,
              text: [
                `Saved ${shot.width}×${shot.height} PNG to ${file}.`,
                shot.gene
                  ? `Coloured by gene ${shot.gene}.`
                  : `Coloured by ${shot.column}: ${(shot.legend ?? []).map((l) => l.value).join(", ")}.`,
              ].join("\n"),
            },
          ],
        };
      } catch (e) {
        return fail(e);
      }
    },
  );
}

export function registerAddColumnTool(server: McpServer) {
  server.registerTool(
    "meyes_view_add_column",
    {
      title: "Add your results to the open view",
      description: `Add columns you computed (pseudotime, new clusters, labels such as "radial glia" / "immature", scores) to the MERFISHEYES view open in the person's browser, without changing their data: they then behave like the dataset's own columns — colour by them, show only or hide values, snapshot.
Write the results to a CSV or TSV first: one row per cell, a header, an id column that matches the dataset's cell ids (pandas' df.to_csv() writes the index as the first column, which is used) and one column per result. Without an id column, give exactly one row per cell in the dataset's cell order. Cells missing from the file get an empty value.
Example: df[["pseudotime","fate"]].to_csv("/tmp/results.csv") then file="/tmp/results.csv", color_by=true.`,
      inputSchema: z.strictObject({
        file: z.string().min(1).max(4096).describe("The CSV or TSV"),
        id_column: z
          .string()
          .max(200)
          .optional()
          .describe(
            "The column with cell ids (default: an unnamed first column, or cell_id / obs_names / barcode / index)",
          ),
        columns: z
          .array(z.string().max(200))
          .max(50)
          .optional()
          .describe("Which columns to add (default: all but the id column)"),
        color_by: z
          .boolean()
          .default(false)
          .describe("Colour by the first added column"),
        replace: z
          .boolean()
          .default(false)
          .describe("Replace a column of the same name already in the view"),
      }),
      annotations: ACT,
    },
    async ({ file, id_column, columns, color_by, replace }) => {
      try {
        const path = resolve(file);
        const { header, rows } = parseTable(readFileSync(path, "utf8"));
        const id = idColumn(header, id_column);
        const wanted = (columns ?? header.filter((_, i) => i !== id)).map(
          (name) => {
            const i = header.indexOf(name);

            if (i < 0)
              throw new Error(
                `No column "${name}" in ${path}; its columns: ${header.join(", ")}.`,
              );

            return { name, i };
          },
        );

        if (!wanted.length)
          throw new Error(`${path} has no column to add besides the ids.`);
        const result = await local.command<{
          added: { name: string; type: string; n_values: number }[];
          cells: number;
          matched: number;
          unmatched: number;
          cells_without_a_value: number;
        }>(
          "add_column",
          {
            ids: id === null ? undefined : rows.map((r) => r[id] ?? ""),
            columns: wanted.map(({ name, i }) => ({
              name,
              values: rows.map((r) => r[i] ?? ""),
            })),
            color_by,
            replace,
          },
          120_000,
        );

        return text(
          [
            `Added to the open view: ${result.added.map((c) => `${c.name} (${c.type}, ${c.n_values} distinct values)`).join("; ")}.`,
            id === null
              ? `Values taken in cell order (${result.cells} cells).`
              : `Matched ${result.matched} of ${rows.length} rows to the dataset's ${result.cells} cells by "${header[id] || "the first column"}"${result.unmatched ? `; ${result.unmatched} ids matched no cell` : ""}${result.cells_without_a_value ? `; ${result.cells_without_a_value} cells have no value` : ""}.`,
            color_by
              ? `Coloured by ${result.added[0].name}.`
              : "Colour by them with meyes_view_set.",
          ].join("\n"),
        );
      } catch (e) {
        return fail(e);
      }
    },
  );
}
