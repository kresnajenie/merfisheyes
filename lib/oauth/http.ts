import { NextResponse } from "next/server";

/**
 * HTTP for the OAuth endpoints: their errors are RFC 6749 §5.2 JSON, and
 * they answer any origin (a browser-based agent client calls them; nothing
 * here reads a cookie, so no credentials are allowed).
 */

export const CORS = {
  "Access-Control-Allow-Origin": "*",
  "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
  "Access-Control-Allow-Headers":
    "Content-Type, Authorization, MCP-Protocol-Version",
  "Access-Control-Max-Age": "86400",
};

export function preflight() {
  return new NextResponse(null, { status: 204, headers: CORS });
}

export function oauthJson(body: unknown, status = 200) {
  return NextResponse.json(body, {
    status,
    headers: { ...CORS, "Cache-Control": "no-store", Pragma: "no-cache" },
  });
}

export function oauthError(error: string, description: string, status = 400) {
  return oauthJson({ error, error_description: description }, status);
}

/** A token or revocation request's parameters: form-encoded as the RFCs say,
 * JSON accepted too. */
export async function formParams(req: Request): Promise<URLSearchParams> {
  const type = req.headers.get("content-type") ?? "";

  if (type.includes("application/json")) {
    const body = (await req.json().catch(() => ({}))) as Record<
      string,
      unknown
    >;

    return new URLSearchParams(
      Object.entries(body)
        .filter(([, v]) => typeof v === "string")
        .map(([k, v]) => [k, v as string]),
    );
  }

  return new URLSearchParams(await req.text());
}
