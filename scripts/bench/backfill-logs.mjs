// Re-parse CloudWatch worker logs for already-recorded server bench results.
//
//   node scripts/bench/backfill-logs.mjs [--env linux-x64]
//
// Useful when the log parser has been improved after a ladder ran (the runner
// process keeps the code it started with): each server-*.json that has a
// batchJobId but a null processMs gets its log re-read and its metrics
// updated in place. Upload/queue numbers are left untouched.
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { execFileSync } from "node:child_process";

const args = process.argv.slice(2);
const flag = (n, d = null) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : d;
};
const ENV_NAME = flag("env", process.env.PERF_ENV || `${os.platform()}-${os.arch()}`);
const DIR = path.join("bench", "results", ENV_NAME);
const REGION = process.env.AWS_REGION || "us-west-2";

const aws = (argv) => {
  try {
    return JSON.parse(execFileSync("aws", [...argv, "--region", REGION, "--output", "json"],
      { encoding: "utf8", maxBuffer: 64 * 1024 * 1024 }));
  } catch {
    return null;
  }
};

// Same parser as scripts/bench/server-run.mjs (with the SM elapsed fallback).
function parseWorkerLog(stream) {
  const events = [];
  let token = null;
  for (let i = 0; i < 200; i++) {
    const argv = ["logs", "get-log-events", "--log-group-name", "/aws/batch/job",
      "--log-stream-name", stream, "--start-from-head"];
    if (token) argv.push("--next-token", token);
    const page = aws(argv);
    if (!page) break;
    events.push(...(page.events ?? []).map((e) => e.message));
    if (!page.nextForwardToken || page.nextForwardToken === token) break;
    token = page.nextForwardToken;
    if ((page.events ?? []).length === 0) break;
  }
  const text = events.join("\n");
  const num = (re) => { const m = text.match(re); return m ? m[1] : null; };
  const elapsedToMs = (s) => {
    if (!s) return null;
    const m = s.match(/(?:(\d+)m\s*)?([\d.]+)s/);
    return m ? Math.round(((Number(m[1] ?? 0) * 60) + Number(m[2])) * 1000) : null;
  };
  let lastElapsed = null;
  for (const m of text.matchAll(/\[((?:\d+m )?[\d.]+s)\] /g)) {
    const ms = elapsedToMs(m[1]);
    if (ms != null && (lastElapsed == null || ms > lastElapsed)) lastElapsed = ms;
  }
  const stages = {};
  for (const m of text.matchAll(/\[([\dhms. ]+)\]\s*=== STEP ([\w.]+): (.+?) ===/g)) {
    stages[`step${m[2]}`] = { at: elapsedToMs(m[1]), label: m[3].split("(")[0].trim() };
  }
  return {
    processMs: elapsedToMs(num(/Total time: ([\dhms. ]+)/)) ?? lastElapsed,
    peakRssGB: Number(num(/Peak RSS: ([\d.]+) GB/)) || null,
    chunkWriteMs: elapsedToMs(num(/Chunk writing done \(([\dhms. ]+)\)/)),
    stages,
    oom: /Out of memory|OutOfMemory|oom_kill/i.test(text),
  };
}

let touched = 0;
for (const f of fs.readdirSync(DIR).filter((n) => n.startsWith("server-") && n.endsWith(".json"))) {
  const p = path.join(DIR, f);
  const r = JSON.parse(fs.readFileSync(p, "utf8"));
  if (!r.batchJobId || r.processMs != null) continue;
  const job = aws(["batch", "describe-jobs", "--jobs", r.batchJobId])?.jobs?.[0];
  const stream = job?.container?.logStreamName ?? job?.attempts?.at(-1)?.container?.logStreamName;
  if (!stream) { console.log(`skip ${r.datasetId}: no log stream`); continue; }
  const parsed = parseWorkerLog(stream);
  Object.assign(r, parsed);
  fs.writeFileSync(p, JSON.stringify(r, null, 1) + "\n");
  touched++;
  console.log(`${r.datasetId}: processMs=${r.processMs} peak=${r.peakRssGB ?? "–"} oom=${parsed.oom}`);
}
console.log(`${touched} result file(s) backfilled`);
