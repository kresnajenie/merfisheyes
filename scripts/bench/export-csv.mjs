// Flatten bench results into one CSV row per dataset, with viewer links.
//
//   node scripts/bench/export-csv.mjs [--env linux-x64] [--base https://dev.merfisheyes.com]
//        [--out bench/results/<env>/results.csv]
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

const args = process.argv.slice(2);
const flag = (n, d = null) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : d;
};
const ENV_NAME = flag("env", process.env.PERF_ENV || `${os.platform()}-${os.arch()}`);
const BASE = flag("base", "https://dev.merfisheyes.com").replace(/\/$/, "");
const DIR = path.join("bench", "results", ENV_NAME);
const OUT = flag("out", path.join(DIR, "results.csv"));

const rows = new Map();
for (const f of fs.readdirSync(DIR).filter((n) => n.endsWith(".json"))) {
  const side = f.startsWith("local-") ? "local" : f.startsWith("server-") ? "server" : null;
  if (!side) continue;
  const r = JSON.parse(fs.readFileSync(path.join(DIR, f)));
  if (!rows.has(r.datasetId)) rows.set(r.datasetId, { id: r.datasetId, format: r.format, sizeGB: r.sizeGB, relevantGB: r.relevantGB, shape: r.shape ?? {} });
  rows.get(r.datasetId)[side] = r;
}

const esc = (v) => {
  if (v == null) return "";
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};
const secs = (ms) => (ms == null ? "" : (ms / 1000).toFixed(1));

const header = [
  "dataset", "format", "file_gb", "pipeline_reads_gb", "cells", "genes", "molecules",
  "browser_outcome", "browser_process_s", "browser_peak_gb", "browser_genes_loaded", "browser_error",
  "server_outcome", "server_process_s", "server_peak_gb", "server_tier",
  "escalated_tier", "escalated_process_s", "escalated_peak_gb",
  "upload_s", "queue_s", "fargate_usd", "output_gb", "s3_usd_month",
  "merfisheyes_link",
];
const lines = [header.join(",")];
for (const r of [...rows.values()].sort((a, b) => a.format.localeCompare(b.format) || a.id.localeCompare(b.id))) {
  const l = r.local, s = r.server, e = s?.escalation, c = s?.cost;
  const serverWorked = s && (s.outcome === "ok" || e?.outcome === "ok");
  const link = serverWorked && s.datasetRowId
    ? `${BASE}/${r.format === "single-molecule" ? "sm-viewer" : "viewer"}/${s.datasetRowId}`
    : "";
  lines.push([
    r.id, r.format, r.sizeGB, r.relevantGB ?? "",
    r.shape.cells ?? "", r.shape.genes ?? "", r.shape.molecules ?? "",
    l?.outcome ?? "", secs(l?.processMs), l?.peakRssGB ?? "",
    l?.stats?.exprGeneCount ?? l?.stats?.geneCount ?? "", (l?.error ?? "").slice(0, 120),
    s?.outcome ?? "", secs(s?.processMs), s?.peakRssGB ?? "", s ? "16GB" : "",
    e ? `${Math.round(Number(e.tier) / 1024)}GB` : "", secs(e?.processMs), e?.peakRssGB ?? "",
    secs(s?.clientUploadMs), secs(s?.queueMs),
    c?.totalUsd ?? "", c?.outputBytes ? (c.outputBytes / 1024 ** 3).toFixed(3) : "", c?.outputUsdMonth ?? "",
    link,
  ].map(esc).join(","));
}
fs.writeFileSync(OUT, lines.join("\n") + "\n");
console.log(`${rows.size} datasets -> ${OUT}`);
