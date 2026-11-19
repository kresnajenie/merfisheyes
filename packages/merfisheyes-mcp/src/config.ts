/**
 * Settings, from the environment. MERFISHEYES_URL points the server at
 * another deployment (staging, or a local dev server).
 */

export const SITE = (
  process.env.MERFISHEYES_URL || "https://www.merfisheyes.com"
).replace(/\/+$/, "");

/** The agent API, and the protected resource tokens are bound to. */
export const API_BASE = `${SITE}/api/agent/v1`;

/**
 * The published query index (MCP.md pillar B): `catalog.parquet` and the
 * files beside it. A folder path works too (tests, a local build).
 */
export const INDEX_URL = (
  process.env.MERFISHEYES_INDEX_URL ||
  "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/bil-psc-zarr/.index/latest"
).replace(/\/+$/, "");

export const NAME = "merfisheyes-mcp";
export const VERSION = "0.2.0";

/** What sign-in asks for when nothing more is needed (MCP-PROTOCOL §6.3). */
export const DEFAULT_SCOPES = ["catalog:read", "projects:read"];

/** A tool result stays under this many characters by returning fewer items. */
export const CHARACTER_LIMIT = 25_000;

/** How long a sign-in link stays open for its answer. */
export const SIGN_IN_TIMEOUT_MS = 10 * 60 * 1000;

/** How long meyes_sign_in waits for the person before handing back the link
 * (an MCP client may give up on a longer tool call). */
export const SIGN_IN_WAIT_MS = Number(
  process.env.MERFISHEYES_SIGN_IN_WAIT_MS || 120_000,
);

/**
 * A Vercel "Protection Bypass for Automation" secret, for a deployment kept
 * private by Vercel's Deployment Protection (the staging site): without it
 * every request is redirected to a Vercel login a program cannot pass. Never
 * needed against the live site. Keep it on your own machine.
 */
export const VERCEL_BYPASS = process.env.MERFISHEYES_VERCEL_BYPASS || "";

/** Headers every request to the site carries. */
export function siteHeaders(): Record<string, string> {
  return VERCEL_BYPASS ? { "x-vercel-protection-bypass": VERCEL_BYPASS } : {};
}
