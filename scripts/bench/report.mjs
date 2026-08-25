// Turn bench results into the comparison table and the claim sentence.
//
//   node scripts/bench/report.mjs [--env linux-x64] [--md bench/RESULTS.md]
//
// Joins the local (browser) and server (Batch) results per dataset and reports
// them side by side. `process` means the same thing on both sides — parse +
// build, excluding byte transfer and Fargate provisioning — so the columns are
// comparable; the excluded phases are printed separately so nothing is hidden.
import fs from "node:fs";
import os from "node:os";
import path from "node:path";

const args = process.argv.slice(2);
const flag = (n, d = null) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : d;
};
const ENV_NAME = flag("env", process.env.PERF_ENV || `${os.platform()}-${os.arch()}`);
const DIR = path.join("bench", "results", ENV_NAME);
const MD = flag("md");

if (!fs.existsSync(DIR)) {
  console.error(`No results in ${DIR}`);
  process.exit(2);
}

const rows = new Map();
for (const f of fs.readdirSync(DIR).filter((n) => n.endsWith(".json"))) {
  const r = JSON.parse(fs.readFileSync(path.join(DIR, f), "utf8"));
  const id = r.datasetId;
  if (!id) continue;
  if (!rows.has(id)) rows.set(id, { id, format: r.format, sizeGB: r.sizeGB, relevantGB: r.relevantGB, shape: r.shape });
  // The record's own "path" key is the dataset path (the runner's spread
  // overwrites its side marker) — the filename prefix is the reliable side.
  rows.get(id)[f.startsWith("server-") ? "server" : "local"] = r;
  if (r.machine?.cpu) rows.get(id).machine = r.machine;
}

const list = [...rows.values()].sort(
  (a, b) => (a.relevantGB || a.sizeGB) - (b.relevantGB || b.sizeGB),
);
const secs = (ms) => (ms == null ? "–" : ms >= 60_000 ? `${(ms / 60000).toFixed(1)}m` : `${(ms / 1000).toFixed(1)}s`);
const mem = (g) => (g == null ? "–" : `${g} GB`);
// A single-cell run that finished with zero genes did not do the same work:
// the browser's Xenium adapter reads no expression when the export ships it as
// cell_feature_matrix.h5 / .tar.gz, so it "succeeds" without any genes. That is
// a degraded result, not a comparable one.
const degraded = (x) =>
  x && x.outcome === "ok" &&
  ((x.stats && x.stats.dataType === "single_cell" && (x.stats.exprGeneCount ?? x.stats.geneCount) === 0) ||
    (x.numGenes != null && x.numGenes <= 1)); // server writes a 1-gene placeholder

const shapeOf = (s) =>
  s?.cells ? `${s.cells.toLocaleString()} cells × ${s.genes ?? "?"}`
  : s?.molecules ? `${s.molecules.toLocaleString()} molecules` : "–";

const machine = [...rows.values()].find((r) => r.machine)?.machine;
const lines = [];
const out = (s = "") => { lines.push(s); console.log(s); };

out(`# Benchmark — ${ENV_NAME}`);
out();
if (machine) {
  out(`**Machine:** ${machine.cpu} · ${machine.cores} cores · ${machine.totalRamGB} GB RAM`
    + (machine.chromium ? ` · Chromium ${machine.chromium}` : ""));
  out();
}
out("`process` = parse + build the chunked layout. Byte transfer and Fargate");
out("provisioning are excluded from it and listed separately.");
out();
out("| dataset | format | size | shape | local | local peak | server | server peak |");
out("|---|---|---:|---|---|---:|---|---:|");
for (const r of list) {
  const l = r.local, s = r.server;
  const tierGB = (t) => Math.round(Number(t) / 1024);
  const esc = (x) => !x?.escalation ? ""
    : x.escalation.outcome === "ok"
      ? ` → ok @${tierGB(x.escalation.tier)} GB ${secs(x.escalation.processMs)}${x.escalation.peakRssGB ? ` (peak ${x.escalation.peakRssGB} GB)` : ""}`
      : ` → oom @${tierGB(x.escalation.tier)} GB`;
  const cell = (x) => !x ? "–"
    : degraded(x) ? `${secs(x.processMs)} ⚠️`
    : x.outcome === "ok" ? secs(x.processMs)
    : `**${x.outcome}**` + esc(x);
  out(`| ${r.id.replace(/^[a-z-]+__/, "")} | ${r.format} | ${(r.relevantGB || r.sizeGB).toFixed(2)} GB | ${shapeOf(r.shape)} `
    + `| ${cell(l)} | ${mem(l?.peakRssGB)} | ${cell(s)} | ${mem(s?.peakRssGB)} |`);
}
out();

const localOk = list.filter((r) => r.local?.outcome === "ok" && !degraded(r.local));
const localBad = list.filter((r) => r.local && r.local.outcome !== "ok");
const serverOk = list.filter((r) => r.server?.outcome === "ok");
const size = (r) => r.relevantGB || r.sizeGB;

if (list.some((r) => degraded(r.local) || degraded(r.server))) {
  out("⚠️ = finished but loaded **no real gene expression**. Neither path reads");
  out("`cell_feature_matrix.h5`, which is how these Xenium exports ship it — the");
  out("browser loads zero genes and the server publishes a 1-gene placeholder.");
  out("Xenium rows are cells-only on both paths, not full comparisons.");
  out();
}

out("## Ceiling");
out();
if (localOk.length) {
  const big = localOk.reduce((a, b) => (size(a) > size(b) ? a : b));
  out(`- Largest processed **in the browser**: \`${big.id}\` — ${size(big).toFixed(2)} GB, `
    + `${shapeOf(big.shape)}, ${secs(big.local.processMs)}, peak ${mem(big.local.peakRssGB)}`);
}
const crashed = localBad.filter((r) => r.local.outcome === "oom");
if (crashed.length) {
  const small = crashed.reduce((a, b) => (size(a) < size(b) ? a : b));
  const peaks = crashed.map((r) => r.local.peakRssGB).filter(Boolean);
  out(`- Smallest **browser crash**: \`${small.id}\` — ${size(small).toFixed(2)} GB, ${shapeOf(small.shape)}`);
  if (peaks.length) {
    out(`- Browser crashes peaked at ${Math.min(...peaks)}–${Math.max(...peaks)} GB`
      + (machine ? ` on a machine with ${machine.totalRamGB} GB of RAM` : "")
      + " — the cap is per-renderer, not the hardware.");
  }
}
if (serverOk.length) {
  const big = serverOk.reduce((a, b) => (size(a) > size(b) ? a : b));
  out(`- Largest processed **on the server**: \`${big.id}\` — ${size(big).toFixed(2)} GB, `
    + `${shapeOf(big.shape)}, ${secs(big.server.processMs)}, peak ${mem(big.server.peakRssGB)}`);
}
out();

const both = list.filter((r) => r.local?.outcome === "ok" && r.server?.outcome === "ok");
if (both.length) {
  out("## Where the server wins (same datasets, both paths succeeded)");
  out();
  out("| dataset | local process | server process | local peak | server peak |");
  out("|---|---:|---:|---:|---:|");
  for (const r of both) {
    out(`| ${r.id.replace(/^[a-z-]+__/, "")} | ${secs(r.local.processMs)} | ${secs(r.server.processMs)} `
      + `| ${mem(r.local.peakRssGB)} | ${mem(r.server.peakRssGB)} |`);
  }
  out();
}

const withUpload = list.filter((r) => r.server?.clientUploadMs);
if (withUpload.length) {
  out("## Excluded from `process` (reported, not hidden)");
  out();
  out("| dataset | upload | MB/s | queue (Fargate) |");
  out("|---|---:|---:|---:|");
  for (const r of withUpload) {
    out(`| ${r.id.replace(/^[a-z-]+__/, "")} | ${secs(r.server.clientUploadMs)} | ${r.server.uploadMBps ?? "–"} | ${secs(r.server.queueMs)} |`);
  }
  out();
  out("Upload time is network-bound and varies by more than 10× between machines,");
  out("which is why it is never part of the headline number.");
  out();
}

if (MD) {
  fs.writeFileSync(MD, lines.join("\n") + "\n");
  console.error(`\nwrote ${MD}`);
}
