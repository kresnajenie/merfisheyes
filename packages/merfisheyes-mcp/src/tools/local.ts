import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { openBrowser } from "../auth.ts";
import { SITE } from "../config.ts";
import { text } from "../format.ts";
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
Opens a SpatialData store (.zarr), a folder holding one, a MERFISHEYES-processed folder (manifest.json) or an OME-Zarr image as they are; an .h5ad (with obsm['spatial']), a Xenium export, a MERSCOPE region, a Visium output or a molecule table (.parquet/.csv with gene, x, y) is first converted once by MERFISHEYES's pipeline into a store in a cache folder (the input is never changed; reused until the input changes). Converting needs Python with the pipeline's packages; a large input takes minutes.
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
        return {
          isError: true,
          content: [{ type: "text" as const, text: (e as Error).message }],
        };
      }
    },
  );
}
