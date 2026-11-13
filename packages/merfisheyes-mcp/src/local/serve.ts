import {
  createReadStream,
  existsSync,
  readdirSync,
  realpathSync,
  statSync,
} from "node:fs";
import { createServer, type Server } from "node:http";
import { randomBytes } from "node:crypto";
import { basename, join, resolve, sep } from "node:path";

import { SITE } from "../config.js";

/**
 * MERFISHEYES on your own data (MCP.md §0, pillar D1): a folder on this
 * computer served on 127.0.0.1 for the site's viewer to read — the browser
 * fetches the files from here; nothing is uploaded.
 *
 * What it allows, and nothing more:
 * - listens on 127.0.0.1 only; refuses a request whose Host is not that
 *   address (a page cannot reach it through a DNS name that rebinds to it);
 * - each folder is reachable only under a random token in the path, and only
 *   inside that folder (symlinks that lead out of it are refused);
 * - GET and HEAD of files, byte ranges (zarr chunks are read in ranges); no
 *   listing of folders;
 * - CORS for the MERFISHEYES site the package talks to (`SITE`) only, with the
 *   answers Chrome's private/local network preflight asks for.
 */

export interface Share {
  /** The folder on disk. */
  path: string;
  /** Its address for the viewer: http://127.0.0.1:<port>/<token>/<name> */
  url: string;
}

const TYPES: Record<string, string> = {
  json: "application/json",
  parquet: "application/vnd.apache.parquet",
  csv: "text/csv",
  png: "image/png",
};

/** A folder's store when the folder only holds one: `sample/` with a single
 * `*.zarr` inside, no manifest and no zarr.json of its own. Otherwise the
 * folder itself (the viewer finds `mapping.json` / `sample/sample.zarr`). */
export function storeIn(path: string): string {
  const has = (f: string) => existsSync(join(path, f));

  if (
    has("zarr.json") ||
    has(".zattrs") ||
    has("manifest.json") ||
    has("manifest.json.gz") ||
    has("mapping.json")
  )
    return path;
  const zarrs = readdirSync(path, { withFileTypes: true }).filter(
    (d) => d.isDirectory() && d.name.endsWith(".zarr"),
  );

  return zarrs.length === 1 ? join(path, zarrs[0].name) : path;
}

/** The first byte range of a `Range` header within `size`, or "bad". */
export function parseRange(
  header: string | undefined,
  size: number,
): { start: number; end: number } | "bad" | null {
  if (!header) return null;
  const m = /^bytes=(\d*)-(\d*)$/.exec(header.split(",")[0].trim());

  if (!m || (m[1] === "" && m[2] === "")) return "bad";
  let start: number;
  let end: number;

  if (m[1] === "") {
    // suffix: the last N bytes
    start = Math.max(0, size - Number(m[2]));
    end = size - 1;
  } else {
    start = Number(m[1]);
    end = m[2] === "" ? size - 1 : Math.min(Number(m[2]), size - 1);
  }
  if (start > end || start >= size) return "bad";

  return { start, end };
}

export class LocalServer {
  private server: Server | null = null;
  private port = 0;
  private readonly roots = new Map<string, { root: string; name: string }>();
  private readonly origin = new URL(SITE).origin;

  /** Serve a folder (starting the server on first use); the same folder
   * twice gives the same address. */
  async share(path: string): Promise<Share> {
    const abs = resolve(path);

    if (!existsSync(abs) || !statSync(abs).isDirectory())
      throw new Error(`${abs} is not a folder on this computer.`);
    const folder = storeIn(abs);
    const root = realpathSync(folder);

    await this.start();
    for (const [token, r] of this.roots)
      if (r.root === root)
        return { path: folder, url: this.address(token, r.name) };
    const token = randomBytes(12).toString("base64url");
    const name = basename(root);

    this.roots.set(token, { root, name });

    return { path: folder, url: this.address(token, name) };
  }

  shared(): Share[] {
    return [...this.roots].map(([token, r]) => ({
      path: r.root,
      url: this.address(token, r.name),
    }));
  }

  async close(): Promise<void> {
    const s = this.server;

    this.server = null;
    if (s) await new Promise<void>((r) => s.close(() => r()));
  }

  private address(token: string, name: string) {
    return `http://127.0.0.1:${this.port}/${token}/${encodeURIComponent(name)}`;
  }

  private async start(port = Number(process.env.MERFISHEYES_LOCAL_PORT || 0)) {
    if (this.server) return;
    const server = createServer((req, res) => this.handle(req, res));

    await new Promise<void>((ok, fail) => {
      server.once("error", fail);
      server.listen(port, "127.0.0.1", () => ok());
    });
    this.server = server;
    this.port = (server.address() as { port: number }).port;
  }

  private handle(
    req: import("node:http").IncomingMessage,
    res: import("node:http").ServerResponse,
  ) {
    const send = (status: number, body = "") => {
      res.writeHead(status, { "content-type": "text/plain" }).end(body);
    };

    // Only this address: a DNS name rebound to 127.0.0.1 is refused.
    if (req.headers.host !== `127.0.0.1:${this.port}`) return send(403);

    res.setHeader("Access-Control-Allow-Origin", this.origin);
    res.setHeader("Vary", "Origin");
    res.setHeader(
      "Access-Control-Expose-Headers",
      "Content-Length, Content-Range, Accept-Ranges",
    );
    res.setHeader("Cache-Control", "no-cache");
    if (req.method === "OPTIONS") {
      res.setHeader("Access-Control-Allow-Methods", "GET, HEAD");
      res.setHeader("Access-Control-Allow-Headers", "range");
      // Chrome's private network preflight
      if (req.headers["access-control-request-private-network"])
        res.setHeader("Access-Control-Allow-Private-Network", "true");

      return res.writeHead(204).end();
    }
    if (req.method !== "GET" && req.method !== "HEAD") return send(405);

    let parts: string[];

    try {
      parts = new URL(req.url ?? "/", "http://x").pathname
        .split("/")
        .filter(Boolean)
        .map(decodeURIComponent);
    } catch {
      return send(400);
    }
    const [token, name, ...rest] = parts;
    const share = token ? this.roots.get(token) : undefined;

    if (
      !share ||
      name !== share.name ||
      rest.some((p) => p === ".." || /[\\/]/.test(p))
    )
      return send(404);
    let file: string;

    try {
      file = realpathSync(join(share.root, ...rest));
    } catch {
      return send(404);
    }
    // Inside the folder only, symlinks included.
    if (file !== share.root && !file.startsWith(share.root + sep))
      return send(404);
    const st = statSync(file);

    if (!st.isFile()) return send(404);
    const range = parseRange(req.headers.range, st.size);

    if (range === "bad") {
      res.setHeader("Content-Range", `bytes */${st.size}`);

      return send(416);
    }
    const ext = file.slice(file.lastIndexOf(".") + 1).toLowerCase();

    res.setHeader("Content-Type", TYPES[ext] ?? "application/octet-stream");
    res.setHeader("Accept-Ranges", "bytes");
    if (range) {
      res.setHeader(
        "Content-Range",
        `bytes ${range.start}-${range.end}/${st.size}`,
      );
      res.setHeader("Content-Length", range.end - range.start + 1);
      res.writeHead(206);
    } else {
      res.setHeader("Content-Length", st.size);
      res.writeHead(200);
    }
    if (req.method === "HEAD") return res.end();
    createReadStream(file, range ?? undefined).pipe(res);
  }
}

/** The link that opens a shared folder in the site's viewer. */
export const openLink = (share: Share) =>
  `${SITE}/open?url=${encodeURIComponent(share.url)}`;

/** One server per process: the CLI's, or the MCP server's for its tools. */
export const local = new LocalServer();
