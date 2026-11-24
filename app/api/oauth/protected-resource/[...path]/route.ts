import { preflight } from "@/lib/oauth/http";
import { protectedResourceMetadata } from "@/lib/oauth/protected-resource";

// GET /.well-known/oauth-protected-resource/<path> (rewritten here, the path
// kept as segments): the metadata of the resource at <path> — `/mcp`, or
// `/api/agent/v1`. (A `?path=` rewrite lost the path and answered every
// path with the agent API's metadata.)
export async function GET(
  _req: Request,
  ctx: { params: Promise<{ path: string[] }> },
) {
  return protectedResourceMetadata((await ctx.params).path.join("/"));
}

export const OPTIONS = preflight;
export const dynamic = "force-dynamic";
