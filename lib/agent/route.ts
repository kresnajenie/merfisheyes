import type { Scope } from "@/lib/oauth/config";

import { type Caller, hasScope, resolveCaller } from "@/lib/agent/caller";
import {
  insufficientScope,
  problem,
  tooManyRequests,
  unauthorized,
} from "@/lib/agent/respond";
import { LIMITS, checkRateLimits, clientIp } from "@/lib/rate-limit";

/**
 * Every /api/agent/v1 handler: resolve the caller, count the call against
 * its rate limits (docs/plans/MCP.md D14), check the scopes the call needs
 * (all of them named in one challenge), run, and turn an unexpected failure
 * into an answer that says what to do — never a stack trace. `expensive`
 * marks a call that reads S3 or asks AWS Batch: it has a tighter limit.
 */
export function agentRoute<
  P extends Record<string, string> = Record<string, never>,
>(
  scopes: Scope[],
  handler: (ctx: {
    req: Request;
    url: URL;
    caller: Caller;
    params: P;
  }) => Promise<Response>,
  opts: { expensive?: boolean } = {},
) {
  return async (req: Request, ctx?: { params: Promise<P> }) => {
    const caller = await resolveCaller(req);

    if (caller instanceof Response) return caller;
    const subject =
      caller.kind === "user" ? `user:${caller.userId}` : `ip:${clientIp(req)}`;
    const limited = await checkRateLimits([
      ...(caller.kind === "user"
        ? LIMITS.agentUser(caller)
        : LIMITS.agentAnonymous(clientIp(req))),
      ...(opts.expensive ? LIMITS.agentExpensive(subject) : []),
    ]);

    if (!limited.ok) return tooManyRequests(limited.retryAfter);
    if (!scopes.every((s) => hasScope(caller, s))) {
      return caller.kind === "anonymous"
        ? unauthorized(null, "sign in: this needs an access token", scopes)
        : insufficientScope(scopes);
    }
    try {
      return await handler({
        req,
        url: new URL(req.url),
        caller,
        params: ctx ? await ctx.params : ({} as P),
      });
    } catch (e) {
      console.error(`agent API ${new URL(req.url).pathname}:`, e);

      return problem(
        500,
        "server_error",
        "the server failed on this request; try again, and report it if it keeps failing",
      );
    }
  };
}
