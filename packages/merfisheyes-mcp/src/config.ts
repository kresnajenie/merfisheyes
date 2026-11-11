/**
 * Settings, from the environment. MERFISHEYES_URL points the server at
 * another deployment (staging, or a local dev server).
 */

export const SITE = (
  process.env.MERFISHEYES_URL || "https://www.merfisheyes.com"
).replace(/\/+$/, "");

/** The agent API, and the protected resource tokens are bound to. */
export const API_BASE = `${SITE}/api/agent/v1`;

export const NAME = "merfisheyes-mcp";
export const VERSION = "0.1.0";

/** What sign-in asks for when nothing more is needed (MCP-PROTOCOL §6.3). */
export const DEFAULT_SCOPES = ["catalog:read", "projects:read"];

/** A tool result stays under this many characters by returning fewer items. */
export const CHARACTER_LIMIT = 25_000;

/** How long a sign-in link stays open for its answer. */
export const SIGN_IN_TIMEOUT_MS = 10 * 60 * 1000;
