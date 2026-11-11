import { mkdirSync, readFileSync, renameSync, writeFileSync } from "node:fs";
import { homedir } from "node:os";
import { dirname, join } from "node:path";

/**
 * Tokens and the client id, kept per authorization server (keyed by its
 * issuer: credentials are never reused against another one, MCP
 * authorization §"Authorization Server Binding"). One JSON file readable by
 * the person only (0600, in a 0700 folder); written to a temporary file and
 * renamed, so a crash never leaves half a file — a rotated refresh token that
 * is lost revokes the whole grant.
 */

export interface Saved {
  clientId: string;
  accessToken?: string;
  refreshToken?: string;
  /** ms since epoch */
  expiresAt?: number;
  scope?: string;
}

export function credentialsPath(): string {
  if (process.env.MERFISHEYES_CREDENTIALS)
    return process.env.MERFISHEYES_CREDENTIALS;
  const base =
    process.env.XDG_CONFIG_HOME ||
    (process.platform === "win32"
      ? process.env.APPDATA || join(homedir(), "AppData", "Roaming")
      : join(homedir(), ".config"));

  return join(base, "merfisheyes-mcp", "credentials.json");
}

function readAll(): Record<string, Saved> {
  try {
    const parsed = JSON.parse(readFileSync(credentialsPath(), "utf8"));

    return parsed && typeof parsed === "object" ? parsed : {};
  } catch {
    return {};
  }
}

export function load(issuer: string): Saved | undefined {
  return readAll()[issuer];
}

/** Replace what is kept for an issuer; undefined forgets it entirely. */
export function store(issuer: string, saved: Saved | undefined): void {
  const all = readAll();

  if (saved) all[issuer] = saved;
  else delete all[issuer];
  const path = credentialsPath();

  mkdirSync(dirname(path), { recursive: true, mode: 0o700 });
  const tmp = `${path}.${process.pid}.tmp`;

  writeFileSync(tmp, JSON.stringify(all, null, 2), { mode: 0o600 });
  renameSync(tmp, path);
}
