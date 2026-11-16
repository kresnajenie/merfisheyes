import { preflight } from "@/lib/oauth/http";
import { protectedResourceMetadata } from "@/lib/oauth/protected-resource";

// GET /.well-known/oauth-protected-resource (rewritten here): the agent
// API's metadata. A resource at a path is served by [...path]/route.ts.
export function GET() {
  return protectedResourceMetadata(null);
}

export const OPTIONS = preflight;
export const dynamic = "force-dynamic";
