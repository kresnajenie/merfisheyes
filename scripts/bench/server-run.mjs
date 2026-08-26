// Benchmark the SERVER path: upload the raw dataset, let AWS Batch process it,
// and read the worker's own stage timings back out of CloudWatch.
//
//   node scripts/bench/server-run.mjs --datasets bench/datasets.json \
//        --base http://localhost:3000 --email you@example.com \
//        [--only <id> ...] [--format h5ad ...] [--max-gb 5] [--conc 8]
//
// The app it talks to must run with ALLOW_DEV_EMAIL_LOGIN=true and AWS Batch
// configured. Uploads mirror lib/ingest/classify-folder.ts exactly: a 67 GB
// Xenium export sends only the handful of files the processor reads, which is
// also what the browser would send.
//
// The headline metric is processMs — the processor's own "Total time", i.e. the
// same parse+build work the browser does. Byte transfer, Fargate provisioning
// and the worker's S3 round trips are recorded separately and never folded in.
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { execFileSync } from "node:child_process";

const args = process.argv.slice(2);
const flag = (n, d = null) => {
  const i = args.indexOf(`--${n}`);
  return i >= 0 && args[i + 1] && !args[i + 1].startsWith("--") ? args[i + 1] : d;
};
const multi = (n) => {
  const out = [];
  for (let i = 0; i < args.length; i++)
    if (args[i] === `--${n}`)
      for (let j = i + 1; j < args.length && !args[j].startsWith("--"); j++) out.push(args[j]);
  return out;
};

const DATASETS = flag("datasets", "bench/datasets.json");
const BASE = flag("base", "http://localhost:3000");
const EMAIL = flag("email", process.env.BENCH_EMAIL);
const ONLY = multi("only");
const FORMATS = multi("format");
const MAX_GB = Number(flag("max-gb", "0")) || Infinity;
const CONC = Number(flag("conc", "8"));
const REGION = process.env.AWS_REGION || "us-west-2";
const ENV_NAME = process.env.PERF_ENV || `${os.platform()}-${os.arch()}`;
const OUT_DIR = path.join("bench", "results", ENV_NAME);

if (!EMAIL) {
  console.error("--email (or BENCH_EMAIL) is required: the owner account to upload as");
  process.exit(2);
}

// ---- the browser's folder filter, mirrored (lib/ingest/classify-folder.ts) ---
const XENIUM_EXPRESSION = new Set([
  "cell_by_gene.csv", "cell_feature_matrix.h5", "features.tsv", "features.tsv.gz",
  "barcodes.tsv", "barcodes.tsv.gz", "matrix.mtx", "matrix.mtx.gz",
]);
const isTop = (key) => !key.includes("/");
const base = (key) => key.split("/").pop().toLowerCase();

function walk(dir, prefix = "") {
  const out = [];
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    if (e.name.startsWith(".")) continue;
    const abs = path.join(dir, e.name);
    const key = prefix ? `${prefix}/${e.name}` : e.name;
    let st;
    try { st = fs.statSync(abs); } catch { continue; }
    if (st.isDirectory()) out.push(...walk(abs, key));
    else out.push({ key, abs, size: st.size });
  }
  return out;
}

function filesToUpload(dataset, rootDir) {
  const abs = path.join(rootDir, dataset.datasetRoot ?? dataset.path);
  const fmt = dataset.declaredFormat;
  if (fs.statSync(abs).isFile()) {
    return [{ key: path.basename(abs), abs, size: fs.statSync(abs).size }];
  }
  const all = walk(abs);
  if (fmt === "xenium") {
    const cells = all.filter((f) => isTop(f.key) && ["cells.csv", "cells.csv.gz"].includes(base(f.key)));
    const expr = all.filter(
      (f) => XENIUM_EXPRESSION.has(base(f.key)) &&
        (isTop(f.key) || f.key.toLowerCase().startsWith("cell_feature_matrix/")),
    );
    return [...cells, ...expr];
  }
  if (fmt === "merscope") {
    return all.filter((f) => {
      const n = base(f.key);
      if (!isTop(f.key) || !n.endsWith(".csv")) return false;
      const meta = n.includes("metadata") && !n.includes("gene");
      const cbg = n.includes("cell") && n.includes("gene");
      return meta || cbg;
    });
  }
  return all;
}

// ---- app session ------------------------------------------------------------
const jar = new Map();
const setCookies = (res) => {
  for (const c of res.headers.getSetCookie?.() ?? []) {
    const [kv] = c.split(";");
    const i = kv.indexOf("=");
    jar.set(kv.slice(0, i).trim(), kv.slice(i + 1));
  }
};
async function api(url, init = {}) {
  const res = await fetch(url, {
    ...init,
    headers: { ...(init.headers || {}), cookie: [...jar].map(([k, v]) => `${k}=${v}`).join("; ") },
    redirect: "manual",
  });
  setCookies(res);
  return res;
}
const asJson = async (res) => {
  const t = await res.text();
  try { return JSON.parse(t); } catch { return { raw: t.slice(0, 300) }; }
};

async function signIn() {
  const { csrfToken } = await asJson(await api(`${BASE}/api/auth/csrf`));
  await api(`${BASE}/api/auth/callback/dev-email`, {
    method: "POST",
    headers: { "content-type": "application/x-www-form-urlencoded" },
    body: new URLSearchParams({ csrfToken, email: EMAIL, callbackUrl: `${BASE}/` }),
  });
  const s = await asJson(await api(`${BASE}/api/auth/session`));
  if (!s?.user?.id) throw new Error(`sign-in failed: ${JSON.stringify(s).slice(0, 200)}`);
  return s.user;
}

// ---- upload -----------------------------------------------------------------
async function putWithRetry(url, body, contentType, attempt = 1) {
  try {
    const res = await fetch(url, { method: "PUT", headers: contentType ? { "content-type": contentType } : {}, body });
    if (!res.ok) throw new Error(`PUT ${res.status}`);
    return res.headers.get("etag");
  } catch (e) {
    if (attempt >= 4) throw e;
    await new Promise((r) => setTimeout(r, 2000 * attempt));
    return putWithRetry(url, body, contentType, attempt + 1);
  }
}

async function uploadOne(file, target, datasetId, uploadId) {
  const fh = await fs.promises.open(file.abs, "r");
  try {
    if (target.mode === "single") {
      const buf = Buffer.alloc(file.size);
      if (file.size) await fh.read(buf, 0, file.size, 0);
      await putWithRetry(target.url, buf, "application/octet-stream");
      return { key: file.key, multipart: undefined };
    }
    const { partSize, parts, s3UploadId } = target;
    const etags = new Array(parts.length);
    let next = 0;
    const worker = async () => {
      while (next < parts.length) {
        const p = next++;
        const start = p * partSize;
        const len = Math.min(partSize, file.size - start);
        const buf = Buffer.alloc(len);
        await fh.read(buf, 0, len, start);
        const etag = await putWithRetry(parts[p].url, buf, null);
        if (!etag) throw new Error(`missing ETag for part ${p + 1} of ${file.key}`);
        etags[p] = { partNumber: parts[p].partNumber, etag };
      }
    };
    await Promise.all(Array.from({ length: Math.min(CONC, parts.length) }, worker));
    return { key: file.key, multipart: { s3UploadId, parts: etags } };
  } finally {
    await fh.close();
  }
}

// ---- AWS --------------------------------------------------------------------
const aws = (argv) => {
  try {
    return JSON.parse(execFileSync("aws", [...argv, "--region", REGION, "--output", "json"], {
      encoding: "utf8", maxBuffer: 64 * 1024 * 1024,
    }));
  } catch {
    return null;
  }
};

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
  const stages = {};
  for (const m of text.matchAll(/\[([\dhms. ]+)\]\s*=== STEP ([\w.]+): (.+?) ===/g)) {
    stages[`step${m[2]}`] = { at: elapsedToMs(m[1]), label: m[3].split("(")[0].trim() };
  }
  // process_single_molecule.py prints no "Total time" summary — fall back to
  // the largest [elapsed] marker its log lines carry.
  let lastElapsed = null;
  for (const m of text.matchAll(/\[((?:\d+m )?[\d.]+s)\] /g)) {
    const ms = elapsedToMs(m[1]);
    if (ms != null && (lastElapsed == null || ms > lastElapsed)) lastElapsed = ms;
  }
  return {
    processMs: elapsedToMs(num(/Total time: ([\dhms. ]+)/)) ?? lastElapsed,
    peakRssGB: Number(num(/Peak RSS: ([\d.]+) GB/)) || null,
    chunkWriteMs: elapsedToMs(num(/Chunk writing done \(([\dhms. ]+)\)/)),
    stages,
    oom: /Out of memory|OutOfMemory|oom_kill/i.test(text),
    tail: events.slice(-3).map((l) => l.slice(0, 200)),
  };
}

async function waitForJob(jobId) {
  for (;;) {
    const d = aws(["batch", "describe-jobs", "--jobs", jobId]);
    const job = d?.jobs?.[0];
    if (!job) throw new Error(`job ${jobId} not found`);
    if (job.status === "SUCCEEDED" || job.status === "FAILED") return job;
    await new Promise((r) => setTimeout(r, 20_000));
  }
}

// ---- one dataset ------------------------------------------------------------
async function runOne(dataset, rootDir) {
  const kind = dataset.declaredFormat === "single-molecule" ? "single_molecule" : "single_cell";
  const files = filesToUpload(dataset, rootDir);
  const uploadBytes = files.reduce((a, f) => a + f.size, 0);
  const result = {
    datasetId: dataset.id, format: dataset.declaredFormat, path: dataset.path,
    sizeGB: dataset.sizeGB, relevantGB: dataset.relevantGB ?? null, shape: dataset.shape,
    uploadFiles: files.length, uploadGB: +(uploadBytes / 1024 ** 3).toFixed(3),
    outcome: "ok", processMs: null, peakRssGB: null,
    clientUploadMs: null, uploadMBps: null, queueMs: null, jobTotalMs: null,
    computeTier: null, batchJobId: null, error: null,
  };

  const init = await asJson(await api(`${BASE}/api/ingest/initiate`, {
    method: "POST",
    headers: { "content-type": "application/json" },
    body: JSON.stringify({
      metadata: { title: `bench ${dataset.id}` },
      processingParams: { kind, stages: { chunk: { chunkSize: 1 } } },
      files: files.map((f) => ({ key: f.key, size: f.size, contentType: "application/octet-stream" })),
    }),
  }));
  if (!init?.datasetId) throw new Error(`initiate failed: ${JSON.stringify(init).slice(0, 200)}`);
  const { datasetId, uploadId, uploadUrls } = init;
  result.datasetRowId = datasetId;

  const t0 = Date.now();
  for (const f of files) {
    const target = uploadUrls[f.key];
    if (!target) throw new Error(`no upload URL for ${f.key}`);
    const done = await uploadOne(f, target, datasetId, uploadId);
    const res = await api(`${BASE}/api/ingest/${datasetId}/files/${encodeURIComponent(f.key)}/complete`, {
      method: "POST", headers: { "content-type": "application/json" },
      body: JSON.stringify(done.multipart ? { uploadId, ...done.multipart } : { uploadId }),
    });
    if (!res.ok) throw new Error(`file complete failed for ${f.key}: ${res.status}`);
  }
  result.clientUploadMs = Date.now() - t0;
  result.uploadMBps = +((uploadBytes / 1e6) / (result.clientUploadMs / 1000)).toFixed(1);

  const comp = await asJson(await api(`${BASE}/api/ingest/${datasetId}/complete`, {
    method: "POST", headers: { "content-type": "application/json" }, body: JSON.stringify({ uploadId }),
  }));
  if (comp?.submitError || !comp?.batchJobId) throw new Error(`submit failed: ${comp?.submitError ?? "no jobId"}`);
  result.batchJobId = comp.batchJobId;

  const job = await waitForJob(comp.batchJobId);
  result.queueMs = job.startedAt && job.createdAt ? job.startedAt - job.createdAt : null;
  result.jobTotalMs = job.stoppedAt && job.startedAt ? job.stoppedAt - job.startedAt : null;
  result.computeTier = job.container?.resourceRequirements?.find((r) => r.type === "MEMORY")?.value ?? null;

  const stream = job.container?.logStreamName ?? job.attempts?.at(-1)?.container?.logStreamName;
  if (stream) Object.assign(result, parseWorkerLog(stream));
  if (job.status === "FAILED") {
    result.outcome = result.oom ? "oom" : "error";
    result.error = job.container?.reason ?? job.statusReason ?? "job failed";
  }
  return result;
}

async function main() {
  const db = JSON.parse(fs.readFileSync(DATASETS, "utf8"));
  let list = db.datasets;
  if (ONLY.length) list = list.filter((d) => ONLY.includes(d.id));
  if (FORMATS.length) list = list.filter((d) => FORMATS.includes(d.declaredFormat));
  // Same skip as the local runner: a folder with no cells and no molecules is
  // one the processor would reject, so uploading it proves nothing.
  const usable = (d) => (d.shape?.cells ?? 0) > 0 || (d.shape?.molecules ?? 0) > 0;
  for (const d of list.filter((x) => !usable(x))) console.log(`skipped ${d.id} — no cells/molecules detected`);
  list = list.filter(usable);
  const rung = (d) => (d.relevantGB ? d.relevantGB : d.sizeGB);
  list = list.filter((d) => rung(d) <= MAX_GB);
  list.sort((a, b) => rung(a) - rung(b));
  if (!list.length) { console.error("No datasets selected."); process.exit(2); }

  const user = await signIn();
  console.log(`signed in as ${user.email}`);
  console.log(`${"dataset".padEnd(46)} ${"upGB".padStart(6)} ${"outcome".padEnd(8)} ${"process".padStart(9)} ${"peakRSS".padStart(8)} ${"upload".padStart(9)} ${"queue".padStart(7)}`);
  console.log("-".repeat(104));

  fs.mkdirSync(OUT_DIR, { recursive: true });
  const machine = { env: ENV_NAME, base: BASE, region: REGION, uploadedFrom: ENV_NAME };
  for (const d of list) {
    let r;
    try {
      r = await runOne({ ...d }, db.root);
    } catch (e) {
      r = { datasetId: d.id, format: d.declaredFormat, outcome: "error", error: String(e?.message ?? e).slice(0, 300) };
    }
    const s = (ms) => (ms != null ? `${(ms / 1000).toFixed(1)}s` : "-");
    console.log(
      `${r.datasetId.slice(0, 46).padEnd(46)} ${String(r.uploadGB ?? "-").padStart(6)} ${r.outcome.padEnd(8)} ` +
      `${s(r.processMs).padStart(9)} ${String(r.peakRssGB ?? "-").padStart(8)} ${s(r.clientUploadMs).padStart(9)} ${s(r.queueMs).padStart(7)}` +
      (r.error ? `\n${" ".repeat(46)}   ! ${r.error}` : ""),
    );
    fs.writeFileSync(path.join(OUT_DIR, `server-${r.datasetId}.json`),
      JSON.stringify({ machine, path: "server", ...r }, null, 1) + "\n");
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
