import type { McpServer } from "@modelcontextprotocol/server";

import * as z from "zod/v4";

import { openBrowser } from "../auth.js";
import { SITE } from "../config.js";
import { text } from "../format.js";
import { local, openLink } from "../local/serve.js";

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
Opens what the viewer reads: a SpatialData store (.zarr, Zarr v3 with consolidated metadata), a folder holding one, a MERFISHEYES-processed folder (manifest.json), or an OME-Zarr image. Not yet: an .h5ad file or raw MERSCOPE/Xenium output — convert those first (scripts/prepare.py in the MERFISHEYES repository writes a SpatialData store).
Returns a link; set open=true only when the person asked to open it now. In Chrome the first time, the person must allow ${new URL(SITE).host} to reach apps on this device.
Examples: "show my sample in MERFISHEYES" -> path="/data/run1/sample.zarr"; "open the output folder" -> path="./out".`,
      inputSchema: z.strictObject({
        path: z
          .string()
          .min(1)
          .max(4096)
          .describe(
            "The folder: absolute, or relative to where this server was started",
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
        const share = await local.share(path);
        const link = openLink(share);

        if (open) openBrowser(link);

        return text(
          [
            `Serving ${share.path} from this computer (nothing is uploaded; only while this MCP server runs).`,
            "",
            link,
            "",
            open
              ? "Opened in a browser tab. If none opened, give the person the link."
              : "Give the person this link; it was not opened.",
            `In Chrome, the first time: allow ${new URL(SITE).host} to reach apps on this device when it asks, or the viewer cannot read the folder.`,
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
