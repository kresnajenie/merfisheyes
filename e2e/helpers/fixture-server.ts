import fs from "node:fs";
import http from "node:http";
import { type AddressInfo } from "node:net";

/**
 * What a path resolves to: a file on disk, a body made up for the test
 * (text, or bytes — optionally trickled out in pieces, so a test can act
 * while a download is still in progress), or nothing (404).
 */
export type Resolved =
  | { file: string }
  | { body: string }
  | { bytes: Buffer; throttle?: { chunkBytes: number; delayMs: number } }
  | undefined;

export interface ServedRequest {
  path: string;
  range?: string;
  status: number;
  /** Bytes of body sent. */
  bytes: number;
  /** Wall time from request to the last byte, ms. */
  ms: number;
  /** performance.now() of the server when the request arrived. */
  at: number;
}

export interface FixtureServer {
  /** Origin of the server, without a trailing slash. */
  url: string;
  /** Every GET served so far, in order. */
  requests: ServedRequest[];
  /** Bytes sent for the requests whose path starts with `prefix`. */
  bytesSince(index: number, prefix?: string): number;
  close(): Promise<void>;
}

export interface FixtureServerOptions {
  /**
   * Cap the rate bodies are sent at, so "a 100 Mbit/s connection" is
   * reproducible in CI: `12_000_000` is about that. Bodies go out in 64 KB
   * pieces paced to the rate; ranged and small responses pay the same rate.
   */
  bytesPerSecond?: number;
}

const PACE_CHUNK = 64 * 1024;

/**
 * A stand-in for a public bucket: a read-only HTTP server on a free local
 * port that allows any origin, answers the preflight a browser sends before a
 * suffix `Range` request, serves byte ranges (which sharded zarr arrays are
 * read with), and tells the browser it may expose timing and sizes for
 * these responses (`Timing-Allow-Origin`), so a test can read bytes per
 * element kind from `performance.getEntriesByType("resource")`. `resolve`
 * maps a URL path to what to serve.
 */
export async function startFixtureServer(
  resolve: (path: string) => Resolved,
  options: FixtureServerOptions = {},
): Promise<FixtureServer> {
  const requests: ServedRequest[] = [];
  const server = http.createServer((req, res) => {
    const startedAt = performance.now();
    const path = decodeURIComponent(new URL(req.url!, "http://x").pathname);

    res.setHeader("Access-Control-Allow-Origin", "*");
    res.setHeader("Timing-Allow-Origin", "*");
    res.setHeader(
      "Access-Control-Expose-Headers",
      "Content-Length, Content-Range",
    );
    if (req.method === "OPTIONS") {
      res.setHeader("Access-Control-Allow-Methods", "GET, HEAD");
      res.setHeader("Access-Control-Allow-Headers", "range");
      res.writeHead(204).end();

      return;
    }

    const found = resolve(path);
    let bytes: Buffer | undefined;
    let throttle: { chunkBytes: number; delayMs: number } | undefined;

    if (found && "body" in found) bytes = Buffer.from(found.body);
    else if (found && "bytes" in found) {
      bytes = found.bytes;
      throttle = found.throttle;
    } else if (
      found &&
      fs.existsSync(found.file) &&
      fs.statSync(found.file).isFile()
    ) {
      bytes = fs.readFileSync(found.file);
    }

    const range = req.headers.range;
    let status = 200;
    let body: Buffer | undefined;

    if (!bytes) {
      status = 404;
    } else if (range) {
      // bytes=a-b, bytes=a- or bytes=-n
      const [, from, to] = /^bytes=(\d*)-(\d*)$/.exec(range)!;
      const start = from ? Number(from) : bytes.length - Number(to);
      const end = from && to ? Number(to) : bytes.length - 1;

      status = 206;
      res.setHeader("Content-Range", `bytes ${start}-${end}/${bytes.length}`);
      body = bytes.subarray(start, end + 1);
    } else {
      body = bytes;
    }

    const record: ServedRequest = {
      path,
      range,
      status,
      bytes: body?.length ?? 0,
      ms: 0,
      at: startedAt,
    };

    if (req.method === "GET") requests.push(record);
    res.on("finish", () => {
      record.ms = performance.now() - startedAt;
    });

    if (!body) {
      res.writeHead(status).end();

      return;
    }
    if (req.method === "HEAD") {
      res.setHeader("Content-Length", body.length);
      res.writeHead(status).end();

      return;
    }

    // Trickle: a fixed piece size with a fixed pause (a test's own timing),
    // or paced to a bandwidth (the same for every response).
    const pace = throttle
      ? throttle
      : options.bytesPerSecond
        ? {
            chunkBytes: PACE_CHUNK,
            delayMs: (PACE_CHUNK / options.bytesPerSecond) * 1000,
          }
        : null;

    if (!pace) {
      res.writeHead(status).end(body);

      return;
    }

    const { chunkBytes, delayMs } = pace;
    const whole = body;
    let offset = 0;

    res.writeHead(status);
    const next = () => {
      if (res.destroyed) return;
      res.write(whole.subarray(offset, offset + chunkBytes));
      offset += chunkBytes;
      if (offset >= whole.length) res.end();
      else setTimeout(next, delayMs);
    };

    next();
  });

  await new Promise<void>((done) => server.listen(0, "127.0.0.1", done));

  return {
    url: `http://127.0.0.1:${(server.address() as AddressInfo).port}`,
    requests,
    bytesSince: (index, prefix = "") =>
      requests
        .slice(index)
        .filter((r) => r.path.startsWith(prefix))
        .reduce((sum, r) => sum + r.bytes, 0),
    close: () => new Promise((done) => server.close(() => done())),
  };
}
