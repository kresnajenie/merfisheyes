import { createHash, randomBytes } from "node:crypto";
import { existsSync, mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { basename, dirname, join } from "node:path";

import { parquetReadObjects } from "hyparquet";
import { compressors } from "hyparquet-compressors";

import { cacheRoot, findPrepare, pipelinePython, run } from "./convert.ts";
import { local, type OpenedDataset } from "./serve.ts";

/**
 * Analysis steps on an open dataset (MCP.md pillar D, tier 2): each runs
 * scripts/meyes_analysis.py in MERFISHEYES's own Python (the one the
 * converter set up, with the person's yes) as a background job, so a long
 * step outlives a tool call: a tool waits a while, then hands back the job
 * for meyes_job. Results go into a workspace per dataset — labels and
 * embeddings in the table's cell order — and into the open view at once; a
 * save writes them into a copy. The person's data is never changed.
 */

export type Step =
  | "profile"
  | "niche"
  | "recluster"
  | "markers"
  | "save"
  | "raw_check";

export interface Job {
  id: string;
  step: Step;
  /** The open dataset the step works on; none for a raw run's check. */
  dataset?: OpenedDataset;
  status: "running" | "done" | "failed";
  started: number;
  lines: string[];
  result?: Record<string, unknown>;
  error?: string;
  /** What was shown in the open view afterwards (or why not). */
  shown?: string[];
  /** The plots to hand back as images (raw check); null: all. */
  want?: string[] | null;
  done: Promise<void>;
}

const jobs = new Map<string, Job>();

/** How long a tool waits for its step before handing back the job (a
 * client may cut a tool call off after a few minutes). */
export const WAIT_MS = Number(process.env.MERFISHEYES_JOB_WAIT_MS || 45_000);

/** The workspace of a dataset: beside the cache, one per store. */
export function workspaceFor(d: OpenedDataset): string {
  const id = createHash("sha256").update(d.store).digest("hex").slice(0, 10);

  return join(
    dirname(cacheRoot()),
    "analysis",
    `${basename(d.store).replace(/\.zarr$/, "")}-${id}`,
  );
}

/** A script of MERFISHEYES's pipeline, beside prepare.py. */
function script(name = "meyes_analysis.py"): string {
  const prepare = findPrepare();
  const s = prepare ? join(dirname(prepare), name) : null;

  if (!s || !existsSync(s))
    throw new Error(
      `MERFISHEYES's ${name} is not available here (it belongs beside prepare.py).`,
    );

  return s;
}

/** Start a step on the open dataset; the job is returned at once, its
 * `done` settles later. */
export async function startStep(
  step: Step,
  args: Record<string, unknown>,
  dataset: OpenedDataset,
): Promise<Job> {
  const py = script();
  const work = workspaceFor(dataset);

  return startJob(step, args, work, dataset, (a, o) => [
    py,
    step,
    "--store",
    dataset.store,
    "--work",
    work,
    "--args",
    a,
    "--out",
    o,
  ]);
}

/**
 * Check a raw MERFISH run folder (RAW-MERFISH.md step 1): summary, flags and
 * plots from scripts/raw_merfish.py, written to a workspace of its own. The
 * run is only read.
 */
export async function startRawCheck(
  folder: string,
  args: Record<string, unknown>,
): Promise<Job> {
  const py = script("raw_merfish.py");
  const id = createHash("sha256").update(folder).digest("hex").slice(0, 10);
  const work = join(
    dirname(cacheRoot()),
    "raw-check",
    `${basename(folder.replace(/[\\/]+$/, "")) || "run"}-${id}`,
  );

  return startJob("raw_check", args, work, undefined, (a, o) => [
    py,
    "check",
    "--run",
    folder,
    "--work",
    work,
    "--args",
    a,
    "--out",
    o,
  ]);
}

async function startJob(
  step: Step,
  args: Record<string, unknown>,
  work: string,
  dataset: OpenedDataset | undefined,
  argv: (argsFile: string, outFile: string) => string[],
): Promise<Job> {
  // throws SetupNeeded before anything starts when Python is not set up
  const { python } = await pipelinePython(
    join(dirname(script()), "requirements.txt"),
  );
  const id = `${step}-${randomBytes(4).toString("hex")}`;

  mkdirSync(join(work, "jobs"), { recursive: true });
  const argsFile = join(work, "jobs", `${id}.args.json`);
  const outFile = join(work, "jobs", `${id}.out.json`);

  writeFileSync(argsFile, JSON.stringify(args));
  const job = {
    id,
    step,
    dataset,
    status: "running",
    started: Date.now(),
    lines: [],
  } as unknown as Job;

  job.done = run(python, argv(argsFile, outFile), (l) => {
    if (/warning|warn\(/i.test(l)) return;
    job.lines.push(l);
    if (job.lines.length > 200) job.lines.splice(0, 100);
  }).then(
    async ({ code }) => {
      if (code === 0 && existsSync(outFile)) {
        job.result = JSON.parse(readFileSync(outFile, "utf8"));
        job.status = "done";
        job.shown = job.dataset
          ? await showInView(job).catch((e) => [
              `Not shown in the open view: ${(e as Error).message}`,
            ])
          : [];
      } else {
        job.status = "failed";
        job.error =
          job.lines
            .filter((l) => l.startsWith("ERROR: "))
            .map((l) => l.slice(7))
            .join(" ") ||
          `The step stopped (exit ${code}): ${job.lines.slice(-6).join(" | ")}`;
      }
    },
    (e) => {
      job.status = "failed";
      job.error = (e as Error).message;
    },
  );
  jobs.set(id, job);

  return job;
}

/** Wait for a job up to `ms`; returns it either way. */
export async function waitFor(job: Job, ms = WAIT_MS): Promise<Job> {
  await Promise.race([job.done, new Promise((r) => setTimeout(r, ms))]);

  return job;
}

export function jobById(id: string): Job {
  const job = jobs.get(id);

  if (!job)
    throw new Error(
      `No job "${id}" in this session (jobs last as long as the MERFISHEYES server).`,
    );

  return job;
}

/** A float32 .npy (as numpy saves it) as rows of [x, y], NaN kept. */
export function readNpy2(file: string): ([number, number] | null)[] {
  const b = readFileSync(file);
  const major = b[6];
  const headerLen = major === 1 ? b.readUInt16LE(8) : b.readUInt32LE(8);
  const start = (major === 1 ? 10 : 12) + headerLen;
  const header = b.subarray(major === 1 ? 10 : 12, start).toString("latin1");

  if (
    !/'descr':\s*'<f4'/.test(header) ||
    /'fortran_order':\s*True/.test(header)
  )
    throw new Error(`Unexpected array layout in ${file}: ${header}`);
  const f = new Float32Array(
    b.buffer.slice(b.byteOffset + start, b.byteOffset + b.length),
  );
  const out: ([number, number] | null)[] = [];

  for (let i = 0; i + 1 < f.length; i += 2)
    out.push(
      Number.isFinite(f[i]) && Number.isFinite(f[i + 1])
        ? [f[i], f[i + 1]]
        : null,
    );

  return out;
}

async function readObs(file: string, name: string): Promise<string[]> {
  const b = readFileSync(file);
  const rows = (await parquetReadObjects({
    file: b.buffer.slice(b.byteOffset, b.byteOffset + b.length),
    columns: [name],
    compressors,
  })) as Record<string, unknown>[];

  return rows.map((r) => String(r[name] ?? ""));
}

/**
 * Put a step's new labels and embedding into the open view, when the view
 * shows this dataset: the person sees them at once, and the agent can
 * colour by them and picture them. Returns what was shown.
 */
async function showInView(job: Job): Promise<string[]> {
  const r = job.result ?? {};

  if (!job.dataset) return [];
  const work = workspaceFor(job.dataset);
  const obs: string[] = [];
  const obsm: string[] = [];

  if (job.step === "niche") obs.push(String(r.column));
  if (job.step === "recluster") {
    for (const c of (r.clusterings ?? []) as { column: string }[])
      obs.push(c.column);
    obsm.push(String(r.embedding));
  }
  if (!obs.length && !obsm.length) return [];
  if (!local.viewOpen)
    return [
      "No view is open: the results are in the workspace and go into a saved copy",
    ];
  const state = await local.command<{ opened?: string }>("state", {}, 120_000);

  if (state.opened !== job.dataset.share.url)
    return [
      "The open view shows another dataset: open this one to see the results.",
    ];
  const shown: string[] = [];

  if (obs.length) {
    const columns = [];

    for (const name of obs)
      columns.push({
        name,
        values: await readObs(join(work, `obs__${name}.parquet`), name),
      });
    await local.command(
      "add_column",
      { columns, replace: true, color_by: job.step === "niche" },
      120_000,
    );
    shown.push(...obs);
  }
  for (const name of obsm) {
    const short = name.replace(/^X_/, "");

    await local.command(
      "add_embedding",
      {
        name: short,
        coords: readNpy2(join(work, `obsm__${name}.npy`)),
        replace: true,
      },
      120_000,
    );
    shown.push(short);
  }

  return shown;
}
