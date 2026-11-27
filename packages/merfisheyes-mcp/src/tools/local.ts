import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { openBrowser } from "../auth.ts";
import { SITE } from "../config.ts";
import { text } from "../format.ts";
import { SetupNeeded, setupPipeline } from "../local/convert.ts";
import { local, openLink } from "../local/serve.ts";

/**
 * meyes_open_local: a folder on the person's computer in the MERFISHEYES
 * viewer (MCP.md §0, pillar D1). The package serves it on 127.0.0.1 for as
 * long as this server runs; the browser reads it from there.
 */
export function registerLocalTools(server: McpServer) {
  server.registerTool(
    "meyes_open_local",
    {
      title: "Open a folder from this computer in the viewer",
      description: `Open data on the person's own computer in the MERFISHEYES viewer, without uploading it: the folder is served on 127.0.0.1 (only to the MERFISHEYES site, only while this server runs) and the viewer reads it from there.
Opens a SpatialData store (.zarr), a folder holding one, a MERFISHEYES-processed folder (manifest.json) or an OME-Zarr image as they are; an .h5ad (with obsm['spatial']), a Xenium export, a MERSCOPE region, a Visium output or a molecule table (.parquet/.csv with gene, x, y) is first converted once by MERFISHEYES's pipeline into a store in a cache folder (the input is never changed; reused until the input changes). The first conversion needs a one-time setup (Python and the pipeline's packages, about 1.3 GB): this tool then answers that setup is needed, without converting — ask the person, and only if they agree call meyes_setup_converter, then call this again. A large input takes minutes to convert.
Returns a link; set open=true only when the person asked to open it now. The viewer itself is served from their computer too when the standalone viewer is installed; otherwise the site's viewer reads the folder, and in Chrome the person must allow ${new URL(SITE).host} to reach apps on this device the first time.
Examples: "show my sample in MERFISHEYES" -> path="/data/run1/sample.zarr"; "open the output folder" -> path="./out".`,
      inputSchema: z.strictObject({
        path: z
          .string()
          .min(1)
          .max(4096)
          .describe(
            "The folder or file: absolute, or relative to where this server was started",
          ),
        open: z
          .boolean()
          .default(false)
          .describe("Open it in a browser tab now; only when the person asked"),
      }),
      annotations: {
        readOnlyHint: true,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: false,
      },
    },
    async ({ path, open }) => {
      try {
        const { share, converted } = await local.open(path);
        const link = openLink(share);

        if (open) openBrowser(link);

        return text(
          [
            ...(converted
              ? [
                  converted.reused
                    ? `${path} was converted before and has not changed since: reusing ${converted.folder}.`
                    : `Converted ${path} with MERFISHEYES's pipeline in ${Math.round(converted.seconds)} s into ${converted.folder} (the input is unchanged).`,
                ]
              : []),
            `Serving ${share.path} from this computer (nothing is uploaded; only while this MCP server runs).`,
            "",
            link,
            "",
            open
              ? "Opened in a browser tab. If none opened, give the person the link."
              : "Give the person this link; it was not opened.",
            local.viewerDir
              ? "The viewer runs from this computer too (the standalone viewer): nothing goes to the internet."
              : `In Chrome, the first time: allow ${new URL(SITE).host} to reach apps on this device when it asks, or the viewer cannot read the folder.`,
          ].join("\n"),
        );
      } catch (e) {
        if (e instanceof SetupNeeded)
          return text(
            [
              `Not converted yet. ${e.message}`,
              "",
              "Ask the person whether that download is OK. Only if they agree, call meyes_setup_converter, then call meyes_open_local again with the same path.",
            ].join("\n"),
          );

        return {
          isError: true,
          content: [{ type: "text" as const, text: (e as Error).message }],
        };
      }
    },
  );

  server.registerTool(
    "meyes_setup_converter",
    {
      title: "Set up MERFISHEYES's converter (one-time download)",
      description: `The one-time setup that lets meyes_open_local convert raw data (.h5ad, Xenium, MERSCOPE, molecule tables) and runs the analysis tools (meyes_profile, meyes_niche, meyes_recluster, meyes_markers, meyes_save_analysis) and meyes_raw_check: downloads Python 3.12 and MERFISHEYES's pipeline packages (about 1.3 GB) into MERFISHEYES's own folder, plus uv (about 15 MB, from its official release, checksum verified) if the computer has none. Nothing is installed system-wide. Takes a minute or two.
Call it ONLY after the person has agreed to the download (meyes_open_local says when it is needed and what it downloads). Does nothing if the setup is already done.`,
      inputSchema: z.strictObject({}),
      annotations: {
        readOnlyHint: false,
        destructiveHint: false,
        idempotentHint: true,
        openWorldHint: true,
      },
    },
    async () => {
      try {
        const lines: string[] = [];
        const started = Date.now();
        const python = await setupPipeline((l) => lines.push(l));

        return text(
          `MERFISHEYES's converter is ready (${Math.round((Date.now() - started) / 1000)} s; Python at ${python}). Call meyes_open_local again to convert and open the data.`,
        );
      } catch (e) {
        return {
          isError: true,
          content: [{ type: "text" as const, text: (e as Error).message }],
        };
      }
    },
  );
}
