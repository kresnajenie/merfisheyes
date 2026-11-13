import { mkdirSync, mkdtempSync, symlinkSync, writeFileSync } from "node:fs";
import { request } from "node:http";
import { tmpdir } from "node:os";
import { join } from "node:path";

import { afterAll, beforeAll, describe, expect, it } from "vitest";

process.env.MERFISHEYES_URL = "https://www.merfisheyes.com";

const { LocalServer, parseRange, storeIn, openLink } = await import(
  "../src/local/serve.js"
);

/** A raw request, so the Host header can be set (fetch would not allow it). */
function get(
  url: string,
  headers: Record<string, string> = {},
  method = "GET",
): Promise<{
  status: number;
  headers: Record<string, string | string[] | undefined>;
  body: string;
}> {
  const u = new URL(url);

  return new Promise((ok, fail) => {
    const req = request(
      { host: "127.0.0.1", port: u.port, path: u.pathname, method, headers },
      (res) => {
        let body = "";

        res.on("data", (c) => (body += c));
        res.on("end", () =>
          ok({ status: res.statusCode!, headers: res.headers, body }),
        );
      },
    );

    req.on("error", fail);
    req.end();
  });
}

let dir: string;
let secret: string;
const srv = new LocalServer();
let base = "";

beforeAll(async () => {
  dir = mkdtempSync(join(tmpdir(), "meyes-local-"));
  const store = join(dir, "sample", "sample.zarr");

  mkdirSync(join(store, "table"), { recursive: true });
  writeFileSync(
    join(store, "zarr.json"),
    '{"attributes":{"spatialdata_attrs":{}}}',
  );
  writeFileSync(join(store, "table", "chunk"), "0123456789");
  secret = join(dir, "secret.txt");
  writeFileSync(secret, "do not serve");
  symlinkSync(secret, join(store, "escape.txt"));
  base = (await srv.share(join(dir, "sample"))).url;
});

afterAll(() => srv.close());

describe("serving a folder from this computer", () => {
  it("points at the one store a folder holds, and links /open on the site", () => {
    expect(base).toMatch(
      /^http:\/\/127\.0\.0\.1:\d+\/[\w-]{16}\/sample\.zarr$/,
    );
    expect(storeIn(join(dir, "sample"))).toBe(
      join(dir, "sample", "sample.zarr"),
    );
    expect(openLink({ path: "", url: base })).toBe(
      `https://www.merfisheyes.com/open?url=${encodeURIComponent(base)}`,
    );
  });

  it("serves files, with CORS for the MERFISHEYES site only", async () => {
    const r = await get(`${base}/zarr.json`, { host: new URL(base).host });

    expect(r.status).toBe(200);
    expect(r.body).toContain("spatialdata_attrs");
    expect(r.headers["access-control-allow-origin"]).toBe(
      "https://www.merfisheyes.com",
    );
    expect(r.headers["content-type"]).toBe("application/json");
  });

  it("answers byte ranges, suffix ranges and a bad range", async () => {
    const h = { host: new URL(base).host };
    const part = await get(`${base}/table/chunk`, { ...h, range: "bytes=2-4" });

    expect(part.status).toBe(206);
    expect(part.body).toBe("234");
    expect(part.headers["content-range"]).toBe("bytes 2-4/10");
    expect(
      (await get(`${base}/table/chunk`, { ...h, range: "bytes=-3" })).body,
    ).toBe("789");
    expect(
      (await get(`${base}/table/chunk`, { ...h, range: "bytes=20-" })).status,
    ).toBe(416);
  });

  it("answers the browser's preflight, private network included", async () => {
    const r = await get(
      `${base}/zarr.json`,
      {
        host: new URL(base).host,
        origin: "https://www.merfisheyes.com",
        "access-control-request-method": "GET",
        "access-control-request-headers": "range",
        "access-control-request-private-network": "true",
      },
      "OPTIONS",
    );

    expect(r.status).toBe(204);
    expect(r.headers["access-control-allow-headers"]).toBe("range");
    expect(r.headers["access-control-allow-private-network"]).toBe("true");
  });

  it("refuses everything outside the folder", async () => {
    const h = { host: new URL(base).host };
    const token = new URL(base).pathname.split("/")[1];
    const origin = new URL(base).origin;

    // no listing, wrong token, wrong name, climbing out, a symlink out
    expect((await get(`${base}/table`, h)).status).toBe(404);
    expect(
      (await get(`${origin}/wrongtoken/sample.zarr/zarr.json`, h)).status,
    ).toBe(404);
    expect((await get(`${origin}/${token}/other/zarr.json`, h)).status).toBe(
      404,
    );
    expect((await get(`${base}/../../secret.txt`, h)).status).toBe(404);
    expect((await get(`${base}/%2E%2E%2F%2E%2E%2Fsecret.txt`, h)).status).toBe(
      404,
    );
    expect((await get(`${base}/escape.txt`, h)).status).toBe(404);
    // writes
    expect((await get(`${base}/zarr.json`, h, "PUT")).status).toBe(405);
  });

  it("refuses a request addressed to another host name (DNS rebinding)", async () => {
    const r = await get(`${base}/zarr.json`, { host: "evil.example.com" });

    expect(r.status).toBe(403);
  });

  it("parses ranges", () => {
    expect(parseRange(undefined, 10)).toBeNull();
    expect(parseRange("bytes=0-", 10)).toEqual({ start: 0, end: 9 });
    expect(parseRange("bytes=5-100", 10)).toEqual({ start: 5, end: 9 });
    expect(parseRange("bytes=-100", 10)).toEqual({ start: 0, end: 9 });
    expect(parseRange("items=0-1", 10)).toBe("bad");
    expect(parseRange("bytes=-", 10)).toBe("bad");
  });
});
