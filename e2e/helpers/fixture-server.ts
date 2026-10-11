import fs from "node:fs";
import http from "node:http";
import { type AddressInfo } from "node:net";

/**
 * What a path resolves to: a file on disk, a body made up for the test, or
 * nothing (404).
 */
export type Resolved = { file: string } | { body: string } | undefined;

export interface FixtureServer {
  /** Origin of the server, without a trailing slash. */
  url: string;
  /** Every GET served so far, in order. */
  requests: { path: string; range?: string; status: number }[];
  close(): Promise<void>;
}

/**
 * A stand-in for a public bucket: a read-only HTTP server on a free local
 * port that allows any origin, answers the preflight a browser sends before a
 * suffix `Range` request, and serves byte ranges (which sharded zarr arrays
 * are read with). `resolve` maps a URL path to what to serve.
 */
export async function startFixtureServer(
  resolve: (path: string) => Resolved,
): Promise<FixtureServer> {
  const requests: FixtureServer["requests"] = [];
  const server = http.createServer((req, res) => {
    const path = decodeURIComponent(new URL(req.url!, "http://x").pathname);

    res.setHeader("Access-Control-Allow-Origin", "*");
    if (req.method === "OPTIONS") {
      res.setHeader("Access-Control-Allow-Methods", "GET, HEAD");
      res.setHeader("Access-Control-Allow-Headers", "range");
      res.writeHead(204).end();

      return;
    }

    const found = resolve(path);
    let bytes: Buffer | undefined;

    if (found && "body" in found) bytes = Buffer.from(found.body);
    else if (
      found &&
      fs.existsSync(found.file) &&
      fs.statSync(found.file).isFile()
    ) {
      bytes = fs.readFileSync(found.file);
    }

    const range = req.headers.range;
    let status = 200;

    if (!bytes) {
      status = 404;
      res.writeHead(status).end();
    } else if (range) {
      // bytes=a-b, bytes=a- or bytes=-n
      const [, from, to] = /^bytes=(\d*)-(\d*)$/.exec(range)!;
      const start = from ? Number(from) : bytes.length - Number(to);
      const end = from && to ? Number(to) : bytes.length - 1;

      status = 206;
      res.writeHead(status, {
        "Content-Range": `bytes ${start}-${end}/${bytes.length}`,
      });
      res.end(bytes.subarray(start, end + 1));
    } else {
      res.writeHead(status).end(bytes);
    }
    if (req.method === "GET") requests.push({ path, range, status });
  });

  await new Promise<void>((done) => server.listen(0, "127.0.0.1", done));

  return {
    url: `http://127.0.0.1:${(server.address() as AddressInfo).port}`,
    requests,
    close: () => new Promise((done) => server.close(() => done())),
  };
}
