// Bundle bench results into the JSON blob the plots page embeds.
//
//   node scripts/bench/benchdata.mjs [--env linux-x64] [--out <file>]
//
// Joins local-*.json / server-*.json in bench/results/<env>/ per dataset
// (archived/ is ignored) and attaches the machine + browser limits from
// limits.json. The output replaces the `const DATA = …;` line in the plots
// page.
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
const OUT = flag("out", path.join(DIR, "benchdata.json"));

const trunc = (s) => (s ? String(s).slice(0, 120) : "");
const side = (r, isServer) => ({
  outcome: r.outcome,
  processMs: r.processMs ?? null,
  peakRssGB: r.peakRssGB ?? null,
  genes: isServer ? r.numGenes ?? null : r.stats?.exprGeneCount ?? r.stats?.geneCount ?? null,
  cells: isServer ? null : (r.stats?.dataType === "single_cell" ? r.stats?.pointCount ?? null : null),
  uploadMs: isServer ? r.clientUploadMs ?? null : null,
  queueMs: isServer ? r.queueMs ?? null : null,
  err: trunc(r.error),
  ...(isServer && r.escalation ? { escalation: r.escalation } : {}),
  ...(isServer && r.cost ? { cost: r.cost } : {}),
});

const rows = new Map();
for (const f of fs.readdirSync(DIR).filter((n) => n.endsWith(".json"))) {
  const isLocal = f.startsWith("local-");
  const isServer = f.startsWith("server-");
  if (!isLocal && !isServer) continue;
  const r = JSON.parse(fs.readFileSync(path.join(DIR, f), "utf8"));
  const id = r.datasetId;
  if (!rows.has(id)) {
    rows.set(id, {
      id, format: r.format, sizeGB: r.sizeGB,
      relevantGB: r.relevantGB ?? null, shape: r.shape ?? null,
    });
  }
  rows.get(id)[isServer ? "server" : "local"] = side(r, isServer);
}

const limits = JSON.parse(fs.readFileSync(path.join(DIR, "limits.json"), "utf8"));
const data = {
  rows: [...rows.values()].sort((a, b) => a.format.localeCompare(b.format) || a.sizeGB - b.sizeGB),
  limits: { maxArrayBufferGB: limits.maxArrayBufferGB, jsHeapLimitGB: limits.jsHeapLimitGB },
  machine: limits.machine,
};
fs.writeFileSync(OUT, JSON.stringify(data));
console.log(`${data.rows.length} rows -> ${OUT}`);
for (const r of data.rows) {
  console.log(`  ${r.format.padEnd(16)} ${r.id.slice(0, 58).padEnd(58)} L:${r.local?.outcome ?? "-"} S:${r.server?.outcome ?? "-"}`);
}
