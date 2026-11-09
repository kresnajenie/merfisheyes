import { NextResponse } from "next/server";

import {
  DEFAULT_SCOPES,
  type Scope,
  resourceMetadataUrl,
} from "@/lib/oauth/config";

/**
 * Answers of the agent API (docs/specs/MCP-PROTOCOL.md §4, §6.5): errors as
 * { error, message } where message says how to recover, the OAuth challenges
 * on 401 and 403, and the paging every list shares.
 */

export function problem(
  status: number,
  error: string,
  message: string,
  headers?: HeadersInit,
) {
  return NextResponse.json({ error, message }, { status, headers });
}

function challenge(params: Record<string, string>) {
  return (
    "Bearer " +
    Object.entries(params)
      .map(([k, v]) => `${k}="${v.replace(/"/g, "'")}"`)
      .join(", ")
  );
}

/** 401: no usable token. */
export function unauthorized(
  error: string | null,
  description: string,
  scopes: Scope[] = DEFAULT_SCOPES,
) {
  const params: Record<string, string> = {
    resource_metadata: resourceMetadataUrl(),
  };

  if (error) params.error = error;
  params.scope = scopes.join(" ");
  if (error) params.error_description = description;

  return problem(401, error ?? "unauthorized", description, {
    "WWW-Authenticate": challenge(params),
  });
}

/** 403: the token is valid but lacks scopes; `scopes` is everything this
 * call needs, in one challenge. */
export function insufficientScope(scopes: Scope[]) {
  const description = `this needs the scope(s) ${scopes.join(", ")}: authorize again asking for them`;

  return problem(403, "insufficient_scope", description, {
    "WWW-Authenticate": challenge({
      error: "insufficient_scope",
      scope: scopes.join(" "),
      resource_metadata: resourceMetadataUrl(),
      error_description: description,
    }),
  });
}

/** 429: over a rate limit; says how long to wait (Retry-After, seconds). */
export function tooManyRequests(retryAfter: number) {
  return problem(
    429,
    "rate_limited",
    `too many requests: wait ${retryAfter} s before trying again, and make fewer calls (page through results instead of repeating them)`,
    { "Retry-After": String(retryAfter) },
  );
}

export function notFound(what: string) {
  return problem(404, "not_found", `no ${what} with that id is visible to you`);
}

export interface Paging {
  limit: number;
  offset: number;
}

export function paging(url: URL, defaultLimit = 20, maxLimit = 100): Paging {
  const n = (k: string, d: number) => {
    const v = Number(url.searchParams.get(k));

    return Number.isFinite(v) && url.searchParams.has(k) ? Math.floor(v) : d;
  };

  return {
    limit: Math.min(maxLimit, Math.max(1, n("limit", defaultLimit))),
    offset: Math.max(0, n("offset", 0)),
  };
}

export function page<T>(items: T[], total: number, p: Paging) {
  const next = p.offset + items.length;
  const hasMore = total > next;

  return {
    total,
    count: items.length,
    offset: p.offset,
    has_more: hasMore,
    ...(hasMore ? { next_offset: next } : {}),
    items,
  };
}
