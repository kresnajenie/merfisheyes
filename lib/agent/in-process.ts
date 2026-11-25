import { type Caller, trustCaller } from "@/lib/agent/caller";
import { issuer } from "@/lib/oauth/config";
import { explainStatus } from "@/packages/merfisheyes-mcp/src/errors.ts";

/**
 * The agent API called from inside the site, as a caller the site's own MCP
 * server (`/mcp`) has already verified: the same route handlers — visibility,
 * scopes, rate limits, answers — without an HTTP round trip (which, behind
 * Vercel's Deployment Protection on staging, would not get through). One
 * implementation of the API; `/mcp` is a second door to it.
 */

/** A route's GET: each declares its own params type, so they are taken as
 * a function of a request and a params promise. */
type Handler = (
  req: Request,

  ctx?: { params: Promise<any> },
) => Promise<Response>;

/** Agent API paths (MCP-PROTOCOL §4) → their handler and parameter names.
 * Fixed segments before `:id` (`/catalog/search` is not a catalog id). */
const ROUTES: [RegExp, string[], () => Promise<{ GET: Handler }>][] = [
  [
    /^\/catalog\/search$/,
    [],
    () => import("@/app/api/agent/v1/catalog/search/route"),
  ],
  [
    /^\/catalog\/facets$/,
    [],
    () => import("@/app/api/agent/v1/catalog/facets/route"),
  ],
  [
    /^\/catalog\/([^/]+)\/genes$/,
    ["id"],
    () => import("@/app/api/agent/v1/catalog/[id]/genes/route"),
  ],
  [
    /^\/catalog\/([^/]+)$/,
    ["id"],
    () => import("@/app/api/agent/v1/catalog/[id]/route"),
  ],
  [
    /^\/datasets\/([^/]+)\/genes$/,
    ["id"],
    () => import("@/app/api/agent/v1/datasets/[id]/genes/route"),
  ],
  [
    /^\/datasets\/([^/]+)$/,
    ["id"],
    () => import("@/app/api/agent/v1/datasets/[id]/route"),
  ],
  [/^\/me$/, [], () => import("@/app/api/agent/v1/me/route")],
  [
    /^\/me\/datasets$/,
    [],
    () => import("@/app/api/agent/v1/me/datasets/route"),
  ],
  [
    /^\/me\/datasets\/([^/]+)\/status$/,
    ["id"],
    () => import("@/app/api/agent/v1/me/datasets/[id]/status/route"),
  ],
  [/^\/projects$/, [], () => import("@/app/api/agent/v1/projects/route")],
  [
    /^\/projects\/([^/]+)$/,
    ["id"],
    () => import("@/app/api/agent/v1/projects/[id]/route"),
  ],
];

export async function agentGet<T>(
  caller: Caller,
  path: string,
  params: Record<string, string | number | undefined>,
  signal?: AbortSignal,
): Promise<T> {
  const route = ROUTES.find(([re]) => re.test(path));

  if (!route) throw new Error(`No agent API route for ${path}.`);
  const [re, names, load] = route;
  const match = re.exec(path)!;
  const routeParams = Object.fromEntries(
    names.map((n, i) => [n, decodeURIComponent(match[i + 1])]),
  );
  const url = new URL(`${issuer()}/api/agent/v1${path}`);

  for (const [k, v] of Object.entries(params))
    if (v !== undefined && v !== "") url.searchParams.set(k, String(v));
  const req = trustCaller(new Request(url, { signal }), caller);
  const res = await (
    await load()
  ).GET(req, {
    params: Promise.resolve(routeParams),
  });

  if (!res.ok) {
    const body = (await res.json().catch(() => ({}))) as {
      error?: string;
      message?: string;
    };

    throw new Error(
      explainStatus(res.status, body, res.headers, {
        signedIn: caller.kind === "user",
        remote: true,
      }),
    );
  }

  return (await res.json()) as T;
}
