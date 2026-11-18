import {
  createReadStream,
  existsSync,
  readdirSync,
  realpathSync,
  statSync,
} from "node:fs";
import {
  createServer,
  type IncomingMessage,
  type Server,
  type ServerResponse,
} from "node:http";
import { randomBytes } from "node:crypto";
import { createRequire } from "node:module";
import { basename, dirname, join, resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";

import { SITE } from "../config.ts";
import { type Converted, convert, viewable } from "./convert.ts";

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
  html: "text/html; charset=utf-8",
  js: "text/javascript; charset=utf-8",
  css: "text/css; charset=utf-8",
  txt: "text/plain; charset=utf-8",
  svg: "image/svg+xml",
  ico: "image/x-icon",
  woff: "font/woff",
  woff2: "font/woff2",
  map: "application/json",
};

/**
 * The standalone viewer's files (pillar D2, `packages/merfisheyes-viewer`,
 * built with `next build`): MERFISHEYES_VIEWER_DIR, else the installed
 * `@merfisheyes/viewer` package, else the copy built in this repository.
 * Null when none is there — then links open the site's viewer instead.
 * MERFISHEYES_VIEWER=site always uses the site.
 */
export function findViewer(): string | null {
  if (process.env.MERFISHEYES_VIEWER === "site") return null;
  const ok = (dir: string) => existsSync(join(dir, "open", "index.html"));
  const candidates: string[] = [];

  if (process.env.MERFISHEYES_VIEWER_DIR)
    candidates.push(resolve(process.env.MERFISHEYES_VIEWER_DIR));
  try {
    candidates.push(
      join(
        dirname(
          createRequire(import.meta.url).resolve(
            "@merfisheyes/viewer/package.json",
          ),
        ),
        "out",
      ),
    );
  } catch {
    // not installed
  }
  // src/local or dist/local → packages/merfisheyes-viewer/out
  candidates.push(
    fileURLToPath(new URL("../../../merfisheyes-viewer/out", import.meta.url)),
  );

  return candidates.find(ok) ?? null;
}

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
  /** The standalone viewer served beside the data, or null (the site's). */
  readonly viewerDir: string | null;
  /** The secret a page needs for the live session (pillar D3): only a page
   * this server opened has it (it is in the link). */
  readonly sessionToken = randomBytes(16).toString("base64url");
  private readonly views: ServerResponse[] = [];
  private readonly pending = new Map<
    number,
    { ok: (v: unknown) => void; fail: (e: Error) => void }
  >();
  private nextId = 1;

  constructor(opts: { viewerDir?: string | null } = {}) {
    this.viewerDir =
      opts.viewerDir === undefined ? findViewer() : opts.viewerDir;
  }

  /** The link that opens a shared folder: the standalone viewer on this
   * same address when there is one, else the site's `/open`. */
  link(share: Share): string {
    const url = encodeURIComponent(share.url);

    return this.viewerDir
      ? `http://127.0.0.1:${this.port}/open/?url=${url}&session=${this.sessionToken}`
      : `${SITE}/open?url=${url}`;
  }

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

  /** Open anything: a folder the viewer reads is shared as it is; anything
   * else (an .h5ad, a Xenium export, a MERSCOPE region, a molecule table)
   * is converted once into the cache (convert.ts) and that is shared. */
  async open(
    path: string,
    onLine?: (line: string) => void,
  ): Promise<{ share: Share; converted: Converted | null }> {
    const abs = resolve(path);

    if (viewable(abs)) return { share: await this.share(abs), converted: null };
    const converted = await convert(abs, { onLine });

    return { share: await this.share(converted.folder), converted };
  }

  /** Whether a page of this server is open and listening for commands. */
  get viewOpen(): boolean {
    return this.views.length > 0;
  }

  /**
   * Send a command to the open view (the most recently opened tab) and wait
   * for its answer: `state`, `set`, `snapshot` (lib/agent/view-commands.ts
   * in the viewer). Throws when no view is open, the view answers with an
   * error, or nothing comes back in time.
   */
  command<T>(command: string, args: unknown, timeoutMs = 30_000): Promise<T> {
    const view = this.views.at(-1);

    if (!view)
      return Promise.reject(
        new Error(
          "No MERFISHEYES view is open: open the data with meyes_open_local, have the person open the link, and keep that tab open.",
        ),
      );
    const id = this.nextId++;

    return new Promise<T>((ok, fail) => {
      const timer = setTimeout(() => {
        this.pending.delete(id);
        fail(
          new Error(`The view did not answer within ${timeoutMs / 1000} s.`),
        );
      }, timeoutMs);

      this.pending.set(id, {
        ok: (v) => {
          clearTimeout(timer);
          ok(v as T);
        },
        fail: (e) => {
          clearTimeout(timer);
          fail(e);
        },
      });
      view.write(`data: ${JSON.stringify({ id, command, args })}\n\n`);
    });
  }

  /** `/session/events` (the page listens) and `/session/reply` (it
   * answers): this server's own pages only — the session token, no CORS. */
  private session(
    path: string,
    url: URL,
    req: IncomingMessage,
    res: ServerResponse,
  ) {
    res.removeHeader("Access-Control-Allow-Origin");
    if (url.searchParams.get("token") !== this.sessionToken)
      return res.writeHead(403).end();
    if (path === "events" && req.method === "GET") {
      res.writeHead(200, {
        "content-type": "text/event-stream",
        "cache-control": "no-cache",
        connection: "keep-alive",
      });
      res.write(": connected\n\n");
      this.views.push(res);
      const beat = setInterval(() => res.write(": .\n\n"), 15_000);

      req.on("close", () => {
        clearInterval(beat);
        const i = this.views.indexOf(res);

        if (i >= 0) this.views.splice(i, 1);
      });

      return;
    }
    if (path === "reply" && req.method === "POST") {
      const chunks: Buffer[] = [];
      let size = 0;

      req.on("data", (c: Buffer) => {
        size += c.length;
        // a snapshot is the largest answer: a few MB of PNG
        if (size > 64 * 1024 * 1024) req.destroy();
        else chunks.push(c);
      });
      req.on("end", () => {
        try {
          const m = JSON.parse(Buffer.concat(chunks).toString("utf8"));
          const waiting = this.pending.get(m.id);

          if (waiting) {
            this.pending.delete(m.id);
            if (m.ok) waiting.ok(m.result);
            else waiting.fail(new Error(m.error ?? "the view failed"));
          }
          res.writeHead(204).end();
        } catch {
          res.writeHead(400).end();
        }
      });

      return;
    }
    res.writeHead(404).end();
  }

  shared(): Share[] {
    return [...this.roots].map(([token, r]) => ({
      path: r.root,
      url: this.address(token, r.name),
    }));
  }

  async close(): Promise<void> {
    const s = this.server;

    for (const v of this.views.splice(0)) v.end();

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

    // Only this computer's own names for itself (what a browser, or an
    // editor forwarding the port, sends): a DNS name rebound to 127.0.0.1
    // is refused.
    const hosts = [
      `127.0.0.1:${this.port}`,
      `localhost:${this.port}`,
      `[::1]:${this.port}`,
    ];

    if (!hosts.includes(req.headers.host ?? "")) return send(403);

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
    {
      const url = new URL(req.url ?? "/", "http://x");
      const m = /^\/session\/(events|reply)$/.exec(url.pathname);

      if (m) return this.session(m[1], url, req, res);
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

    if (!share) return this.viewerFile(parts, req, res);
    if (name !== share.name || rest.some((p) => p === ".." || /[\\/]/.test(p)))
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
    this.sendFile(file, st.size, req, res);
  }

  /** A file of the standalone viewer: `/open/` → `open/index.html`, `/` →
   * the open page; nothing outside the viewer's folder. */
  private viewerFile(
    parts: string[],
    req: import("node:http").IncomingMessage,
    res: import("node:http").ServerResponse,
  ) {
    const notFound = () =>
      res.writeHead(404, { "content-type": "text/plain" }).end();

    if (!this.viewerDir || parts.some((p) => p === ".." || /[\\/]/.test(p)))
      return notFound();
    if (parts.length === 0)
      return res.writeHead(302, { location: "/open/" }).end();
    let file = join(this.viewerDir, ...parts);

    if (existsSync(file) && statSync(file).isDirectory())
      file = join(file, "index.html");
    if (!existsSync(file) || !statSync(file).isFile()) return notFound();
    this.sendFile(file, statSync(file).size, req, res);
  }

  private sendFile(
    file: string,
    size: number,
    req: import("node:http").IncomingMessage,
    res: import("node:http").ServerResponse,
  ) {
    const st = { size };
    const send = (status: number) =>
      res.writeHead(status, { "content-type": "text/plain" }).end();
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

/** The link that opens a shared folder (see LocalServer.link). */
export const openLink = (share: Share) => local.link(share);

/** One server per process: the CLI's, or the MCP server's for its tools. */
export const local = new LocalServer();
