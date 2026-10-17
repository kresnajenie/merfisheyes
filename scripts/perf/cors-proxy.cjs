#!/usr/bin/env node
// Read-only proxy in front of a public bucket that adds the CORS header the
// bucket grants only to the standard dev ports, and logs every response's
// bytes and timing — so the budgets (docs/specs/BUDGETS.md) can be measured
// on the real stores from any local dev server (scripts/perf/real-store.ts).
//
//   node scripts/perf/cors-proxy.cjs [port=3977] [bucket=https://merfisheyes-bil.s3.us-west-2.amazonaws.com]
//
// Range requests pass through. A mapping.json has its absolute bucket URLs
// rewritten to this proxy so a linked molecule folder is reached the same
// way. GET /__log returns the requests served so far as JSON; GET /__reset
// clears it.
const http = require("node:http");
const https = require("node:https");

const PORT = Number(process.argv[2] || 3977);
const BUCKET = (
  process.argv[3] || "https://merfisheyes-bil.s3.us-west-2.amazonaws.com"
).replace(/\/+$/, "");

let log = [];

const server = http.createServer((req, res) => {
  res.setHeader("Access-Control-Allow-Origin", "*");
  res.setHeader("Timing-Allow-Origin", "*");
  res.setHeader("Access-Control-Expose-Headers", "Content-Range, Content-Length");
  if (req.method === "OPTIONS") {
    res.setHeader("Access-Control-Allow-Methods", "GET, HEAD");
    res.setHeader("Access-Control-Allow-Headers", "range");
    res.writeHead(204).end();
    return;
  }
  if (req.url === "/__log") {
    res.setHeader("Content-Type", "application/json");
    res.writeHead(200).end(JSON.stringify(log));
    return;
  }
  if (req.url === "/__reset") {
    log = [];
    res.writeHead(204).end();
    return;
  }
  const startedAt = Date.now();
  const entry = { path: decodeURIComponent(req.url), range: req.headers.range, status: 0, bytes: 0, ms: 0, at: startedAt };
  const headers = {};
  if (req.headers.range) headers.Range = req.headers.range;
  const upstream = https.request(BUCKET + req.url, { method: req.method, headers }, (up) => {
    const passthrough = ["content-type", "content-length", "content-range", "accept-ranges", "last-modified", "etag"];
    const isMapping = req.url.endsWith("/mapping.json") && up.statusCode === 200;
    for (const h of passthrough) {
      if (up.headers[h] && !(isMapping && h === "content-length")) res.setHeader(h, up.headers[h]);
    }
    entry.status = up.statusCode;
    if (!isMapping) {
      res.writeHead(up.statusCode);
      up.on("data", (c) => (entry.bytes += c.length));
      up.on("end", () => {
        entry.ms = Date.now() - startedAt;
        log.push(entry);
      });
      up.pipe(res);
      return;
    }
    let body = "";
    up.setEncoding("utf8");
    up.on("data", (c) => (body += c));
    up.on("end", () => {
      const rewritten = body.split(BUCKET).join(`http://127.0.0.1:${PORT}`);
      entry.bytes = Buffer.byteLength(rewritten);
      entry.ms = Date.now() - startedAt;
      log.push(entry);
      res.writeHead(200);
      res.end(rewritten);
    });
  });
  upstream.on("error", (e) => {
    entry.status = 502;
    entry.ms = Date.now() - startedAt;
    log.push(entry);
    res.writeHead(502).end(String(e));
  });
  upstream.end();
});

server.listen(PORT, "127.0.0.1", () => {
  console.log(`cors proxy for ${BUCKET} on http://127.0.0.1:${PORT}`);
});
