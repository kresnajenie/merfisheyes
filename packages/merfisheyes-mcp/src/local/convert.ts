import { spawn } from "node:child_process";
import { createHash } from "node:crypto";
import {
  existsSync,
  mkdirSync,
  readFileSync,
  readdirSync,
  renameSync,
  rmSync,
  statSync,
  writeFileSync,
} from "node:fs";
import { homedir, release } from "node:os";
import { basename, dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

/**
 * D1b, first part (MCP.md pillar D): data the viewer cannot read as it is —
 * an .h5ad, a Xenium export, a MERSCOPE region, a molecule table — is
 * converted once by MERFISHEYES's own pipeline (`scripts/prepare.py`, the
 * same the server runs) into a SpatialData store in a cache folder, and that
 * store is opened. The input is never changed. A conversion is reused until
 * the input changes.
 */

/** Can the viewer open this folder as it is? (A store, a folder holding
 * one, a MERFISHEYES folder.) */
export function viewable(path: string): boolean {
  if (!existsSync(path) || !statSync(path).isDirectory()) return false;
  const has = (f: string) => existsSync(join(path, f));

  if (
    has("zarr.json") ||
    has(".zattrs") ||
    has("manifest.json") ||
    has("manifest.json.gz") ||
    has("mapping.json")
  )
    return true;

  return readdirSync(path, { withFileTypes: true }).some(
    (d) =>
      d.isDirectory() &&
      d.name.endsWith(".zarr") &&
      existsSync(join(path, d.name, "zarr.json")),
  );
}

/** Where converted stores are kept: MERFISHEYES_CACHE, else the platform's
 * cache folder. */
export function cacheRoot(): string {
  if (process.env.MERFISHEYES_CACHE)
    return resolve(process.env.MERFISHEYES_CACHE);
  const base =
    process.platform === "win32"
      ? process.env.LOCALAPPDATA || join(homedir(), "AppData", "Local")
      : process.platform === "darwin"
        ? join(homedir(), "Library", "Caches")
        : process.env.XDG_CACHE_HOME || join(homedir(), ".cache");

  return join(base, "merfisheyes", "converted");
}

/** The pipeline's scripts: MERFISHEYES_SCRIPTS_DIR, the copy the npm
 * package carries (`python/scripts`), or the repository's `scripts/`. */
export function findPrepare(): string | null {
  const candidates = [
    process.env.MERFISHEYES_SCRIPTS_DIR &&
      resolve(process.env.MERFISHEYES_SCRIPTS_DIR),
    // src/local or dist/local → the package's python/scripts
    fileURLToPath(new URL("../../python/scripts", import.meta.url)),
    // → the repository's scripts/
    fileURLToPath(new URL("../../../../scripts", import.meta.url)),
  ].filter(Boolean) as string[];
  const dir = candidates.find((d) => existsSync(join(d, "prepare.py")));

  return dir ? join(dir, "prepare.py") : null;
}

/** What identifies the input's current contents: size and modification
 * time of the file, or of a folder and everything one level inside it. */
function fingerprint(path: string): string {
  const st = statSync(path);
  const parts = [`${st.size}:${st.mtimeMs}`];

  if (st.isDirectory())
    for (const name of readdirSync(path).sort()) {
      const s = statSync(join(path, name));

      parts.push(`${name}:${s.size}:${s.mtimeMs}`);
    }

  return createHash("sha256").update(parts.join("|")).digest("hex");
}

/** Run a program; every line it prints goes to `onLine`. Resolves with its
 * exit code and output; rejects only when it cannot be started. */
export function run(
  cmd: string,
  args: string[],
  onLine?: (line: string) => void,
): Promise<{ code: number; output: string[] }> {
  return new Promise((done, fail) => {
    const output: string[] = [];
    const child = spawn(cmd, args, { stdio: ["ignore", "pipe", "pipe"] });
    const take = (b: Buffer) =>
      b
        .toString()
        .split("\n")
        .filter(Boolean)
        .forEach((l) => {
          output.push(l);
          onLine?.(l);
        });

    child.stdout.on("data", take);
    child.stderr.on("data", take);
    child.on("error", fail);
    child.on("close", (c) => done({ code: c ?? 1, output }));
  });
}

/** The uv MERFISHEYES downloads when there is none (pinned; checksummed). */
export const UV_VERSION = "0.13.9";

/** The setup is not done and needs the person's yes: what it would
 * download, and where. */
export class SetupNeeded extends Error {
  readonly needsSetup = true;
}

/** Where MERFISHEYES keeps what it sets up for itself, beside the cache of
 * converted stores: its uv and its Python environment. */
const toolsRoot = () => dirname(cacheRoot());

/** uv: MERFISHEYES_UV, the one MERFISHEYES downloaded, `uv` on the PATH, or
 * where uv's installers put it — an app started from the Dock or Start menu
 * (Claude Desktop) does not see the PATH a terminal has. */
export async function findUv(): Promise<string | null> {
  const home = homedir();
  const exe = process.platform === "win32" ? "uv.exe" : "uv";
  const candidates = [
    process.env.MERFISHEYES_UV,
    join(toolsRoot(), `uv-${UV_VERSION}`, exe),
    "uv",
    join(home, ".local", "bin", exe),
    join(home, ".cargo", "bin", exe),
    "/opt/homebrew/bin/uv",
    "/usr/local/bin/uv",
  ].filter(Boolean) as string[];

  for (const uv of candidates) {
    if (/[\\/]/.test(uv) && !existsSync(uv)) continue;
    const ok = await run(uv, ["--version"]).then(
      (r) => r.code === 0,
      () => false,
    );

    if (ok) return uv;
  }

  return null;
}

/** The uv release asset for this computer. */
export function uvAsset(
  platform = process.platform,
  arch = process.arch,
): string | null {
  const cpu = arch === "arm64" ? "aarch64" : arch === "x64" ? "x86_64" : null;

  if (!cpu) return null;
  if (platform === "darwin") return `uv-${cpu}-apple-darwin.tar.gz`;
  if (platform === "linux") return `uv-${cpu}-unknown-linux-gnu.tar.gz`;
  if (platform === "win32") return `uv-${cpu}-pc-windows-msvc.zip`;

  return null;
}

/** Download uv (UV_VERSION, from its GitHub release, checksum verified)
 * into MERFISHEYES's own folder; nothing is installed system-wide. */
async function downloadUv(onLine?: (line: string) => void): Promise<string> {
  const asset = uvAsset();

  if (!asset)
    throw new Error(
      `MERFISHEYES cannot set up uv on ${process.platform}/${process.arch}; install uv yourself (https://docs.astral.sh/uv/).`,
    );
  // (MERFISHEYES_UV_RELEASE_URL: a stand-in release, for tests)
  const release =
    process.env.MERFISHEYES_UV_RELEASE_URL ||
    `https://github.com/astral-sh/uv/releases/download/${UV_VERSION}`;
  const base = `${release}/${asset}`;
  const dir = join(toolsRoot(), `uv-${UV_VERSION}`);
  const exe = join(dir, process.platform === "win32" ? "uv.exe" : "uv");

  onLine?.(`Downloading uv ${UV_VERSION} (about 15 MB)…`);
  const [archive, sum] = await Promise.all([
    fetch(base).then((r) => {
      if (!r.ok) throw new Error(`uv download answered HTTP ${r.status}`);

      return r.arrayBuffer();
    }),
    fetch(`${base}.sha256`).then((r) => r.text()),
  ]);
  const expected = sum.trim().split(/\s+/)[0];
  const got = createHash("sha256").update(Buffer.from(archive)).digest("hex");

  if (got !== expected)
    throw new Error(
      `The uv download did not match its published checksum; nothing was installed.`,
    );
  const work = `${dir}.tmp-${process.pid}`;

  rmSync(work, { recursive: true, force: true });
  mkdirSync(work, { recursive: true });
  const file = join(work, asset);

  writeFileSync(file, Buffer.from(archive));
  // tar reads .tar.gz everywhere and .zip on Windows 10+
  const untar = await run("tar", ["-xf", file, "-C", work]);

  if (untar.code !== 0)
    throw new Error(`Could not unpack uv: ${untar.output.slice(-3).join(" ")}`);
  const inner = join(work, asset.replace(/\.(tar\.gz|zip)$/, ""));
  const found = [join(inner, basename(exe)), join(work, basename(exe))].find(
    (f) => existsSync(f),
  );

  if (!found) throw new Error("The uv archive did not hold a uv program.");
  rmSync(dir, { recursive: true, force: true });
  mkdirSync(dir, { recursive: true });
  renameSync(found, exe);
  rmSync(work, { recursive: true, force: true });

  return exe;
}

/**
 * The pins this computer cannot install as they are, replaced: the
 * pipeline's requirements are locked on Linux (the server's), and a few
 * packages have no ready-made build for some Macs (llvmlite and numba for
 * Intel Macs; pyproj before macOS 15/14). `overrides/<platform>-<arch>.txt`,
 * applied with `uv pip install --override`; null where none is needed.
 */
export function overridesFile(
  platform: string = process.platform,
  arch: string = process.arch,
): string | null {
  const file = fileURLToPath(
    new URL(`../../overrides/${platform}-${arch}.txt`, import.meta.url),
  );

  return existsSync(file) ? file : null;
}

/** Below macOS 13 the pipeline's packages have no ready-made builds. */
export function unsupportedMac(
  platform: string = process.platform,
  darwinMajor = Number(release().split(".")[0]),
): boolean {
  // Darwin 22 = macOS 13 (Ventura)
  return platform === "darwin" && darwinMajor < 22;
}

/** MERFISHEYES's own Python environment for these requirements (and this
 * computer's overrides: a change to either builds a new one). */
function managedEnv(requirements: string) {
  const overrides = overridesFile();
  const pins = createHash("sha256")
    .update(readFileSync(requirements))
    .update(overrides ? readFileSync(overrides) : "")
    .digest("hex")
    .slice(0, 10);
  const env = join(toolsRoot(), `python-${pins}`);
  const python =
    process.platform === "win32"
      ? join(env, "Scripts", "python.exe")
      : join(env, "bin", "python");

  return { env, python, ready: existsSync(join(env, ".merfisheyes-ready")) };
}

/** What the one-time setup downloads, in words for the person. */
export function setupDescription(requirements: string): string {
  return `MERFISHEYES needs a one-time setup of its own Python, to convert raw data (.h5ad, Xenium, MERSCOPE, molecule tables) and to run its analysis tools (meyes_profile … meyes_save_analysis, meyes_raw_check): it downloads Python 3.12 and its pipeline's packages (about 1.3 GB) into ${managedEnv(requirements).env}, plus uv (about 15 MB, the tool that installs them) if this computer has none. Nothing is installed system-wide; viewing data that is already a .zarr store needs none of it. Ask the person, then call meyes_setup_converter.`;
}

/**
 * The Python that runs the pipeline: MERFISHEYES_PYTHON when set (an
 * extension's unset setting can arrive as its "${…}" placeholder); else
 * MERFISHEYES's own environment once it is set up. Otherwise SetupNeeded —
 * the setup downloads a lot, so it runs only when the person agrees
 * (`setupPipeline`; the meyes_setup_converter tool; `merfisheyes setup`).
 */
export async function pipelinePython(
  requirements: string,
): Promise<{ python: string; managed: boolean }> {
  const chosen = process.env.MERFISHEYES_PYTHON;

  if (chosen && !chosen.startsWith("${"))
    return { python: chosen, managed: false };
  const managed = managedEnv(requirements);

  if (managed.ready) return { python: managed.python, managed: true };
  throw new SetupNeeded(setupDescription(requirements));
}

/** The requirements file of the pipeline this package runs. */
export function pipelineRequirements(): string {
  const prepare = findPrepare();

  if (!prepare)
    throw new Error(
      "MERFISHEYES's pipeline (prepare.py) is not available here. Set MERFISHEYES_SCRIPTS_DIR to the scripts/ folder of a MERFISHEYES checkout.",
    );

  return join(prepare, "..", "requirements.txt");
}

/**
 * The one-time setup, once the person agreed: uv (downloaded into
 * MERFISHEYES's folder when this computer has none), then Python 3.12 and
 * the pipeline's pinned packages in MERFISHEYES's own environment. Returns
 * that environment's Python.
 */
export async function setupPipeline(
  onLine?: (line: string) => void,
): Promise<string> {
  const requirements = pipelineRequirements();
  const managed = managedEnv(requirements);

  if (managed.ready) return managed.python;
  if (unsupportedMac())
    throw new Error(
      "Converting raw data needs macOS 13 (Ventura) or later: the pipeline's packages have no ready-made builds for older macOS. Opening data that is already a .zarr store still works.",
    );
  const uv = (await findUv()) ?? (await downloadUv(onLine));
  const overrides = overridesFile();

  onLine?.(
    "Setting up MERFISHEYES's Python (Python 3.12 and the pipeline's packages, about 1.3 GB; a minute or two)…",
  );
  rmSync(managed.env, { recursive: true, force: true });
  mkdirSync(toolsRoot(), { recursive: true });
  const venv = await run(
    uv,
    ["venv", "--quiet", "--python", "3.12", managed.env],
    onLine,
  );

  if (venv.code !== 0)
    throw new Error(
      `Could not create MERFISHEYES's Python environment with uv:\n${venv.output.slice(-6).join("\n")}`,
    );
  const pip = await run(
    uv,
    [
      "pip",
      "install",
      "--quiet",
      "--python",
      managed.python,
      "-r",
      requirements,
      ...(overrides ? ["--override", overrides] : []),
    ],
    onLine,
  );

  if (pip.code !== 0) {
    rmSync(managed.env, { recursive: true, force: true });
    throw new Error(
      `Could not install the pipeline's packages with uv:\n${pip.output.slice(-8).join("\n")}`,
    );
  }
  writeFileSync(join(managed.env, ".merfisheyes-ready"), requirements);

  return managed.python;
}

export interface Converted {
  /** The folder to open (it holds `mapping.json` and the store). */
  folder: string;
  /** True when an earlier conversion of the same input was reused. */
  reused: boolean;
  seconds: number;
}

const EXIT: Record<number, string> = {
  2: "MERFISHEYES does not recognise this as spatial data: it opens a SpatialData store, an .h5ad (with obsm['spatial']), a Xenium export, a MERSCOPE region, a Visium output, a raw MERFISH run (final_spots/ + final_mats/) or a molecule table (gene, x, y columns)",
  3: "the input is incomplete or unusable",
  4: "the converted store did not pass validation",
};

export interface ConvertOptions {
  onLine?: (line: string) => void;
  /** Where the converted folder goes (default: the cache). Must be new or
   * empty, or hold an earlier conversion of the same input. */
  out?: string;
  /** A raw MERFISH run: only these sets. */
  sets?: string[];
  /** A raw MERFISH run: an .h5ad whose labels join by cell name. */
  labelsFrom?: string;
}

/** A raw MERFISH run (docs/specs/raw-merfish-profile.md): its conversion
 * can resume where an interrupted one stopped. */
function isRawRun(path: string): boolean {
  return existsSync(join(path, "final_spots"));
}

/**
 * Convert `input` into a store (in the cache unless `out` says where), or
 * reuse the last conversion of the same, unchanged input with the same
 * options. A raw MERFISH run that was interrupted (the server stopped)
 * resumes from its `.partial` folder. Throws an Error saying what to do when
 * Python or its packages are missing, or the input is not spatial data.
 */
export async function convert(
  input: string,
  opts: ConvertOptions = {},
): Promise<Converted> {
  const abs = resolve(input);

  if (!existsSync(abs)) throw new Error(`${abs} does not exist.`);
  const options = {
    sets: opts.sets?.length ? opts.sets.map(String) : undefined,
    labelsFrom: opts.labelsFrom ? resolve(opts.labelsFrom) : undefined,
  };
  const keyed = options.sets || options.labelsFrom;
  const id = createHash("sha256")
    .update(keyed ? `${abs}|${JSON.stringify(options)}` : abs)
    .digest("hex")
    .slice(0, 10);
  const name = basename(abs).replace(/\.(h5ad|parquet|csv|tsv|zarr)$/i, "");
  const folder = opts.out
    ? resolve(opts.out)
    : join(cacheRoot(), `${name}-${id}`);
  const stamp = join(folder, ".merfisheyes-source.json");
  const print = fingerprint(abs);
  const same = (file: string) => {
    try {
      const was = JSON.parse(readFileSync(file, "utf8"));

      return (
        was.input === abs &&
        was.fingerprint === print &&
        JSON.stringify(was.options ?? {}) === JSON.stringify(options)
      );
    } catch {
      return false;
    }
  };

  if (existsSync(stamp) && same(stamp))
    return { folder, reused: true, seconds: 0 };
  if (
    opts.out &&
    existsSync(folder) &&
    readdirSync(folder).length &&
    !existsSync(stamp)
  )
    throw new Error(
      `${folder} is not empty and holds no earlier conversion of ${abs}: give a new or empty folder.`,
    );

  const prepare = findPrepare();

  if (!prepare)
    throw new Error(
      "MERFISHEYES's pipeline (prepare.py) is not available here. Set MERFISHEYES_SCRIPTS_DIR to the scripts/ folder of a MERFISHEYES checkout.",
    );
  const requirements = join(prepare, "..", "requirements.txt");
  const { python, managed } = await pipelinePython(requirements);
  // One work folder per output, kept when the server stops mid-way: a raw
  // run resumes there; anything else starts again
  const work = `${folder}.partial`;
  const marker = `${work}.json`; // beside it: prepare.py wants the folder empty
  const resume = isRawRun(abs) && existsSync(marker) && same(marker);

  if (!resume) rmSync(work, { recursive: true, force: true });
  mkdirSync(dirname(work), { recursive: true });
  writeFileSync(
    marker,
    JSON.stringify({ input: abs, fingerprint: print, options }, null, 2),
  );
  if (resume) opts.onLine?.(`Resuming the conversion in ${work}`);
  const started = Date.now();
  const args = [
    prepare,
    abs,
    work,
    ...(options.sets ? ["--sets", options.sets.join(",")] : []),
    ...(options.labelsFrom ? ["--labels-from", options.labelsFrom] : []),
    ...(resume ? ["--resume"] : []),
  ];
  const { code, output } = await run(python, args, opts.onLine).catch(
    (e: NodeJS.ErrnoException) => {
      throw e.code === "ENOENT"
        ? new Error(
            `No Python found at "${python}" (MERFISHEYES_PYTHON). Unset it to let MERFISHEYES set up its own, or point it at a Python 3.12 with the pipeline's packages (pip install -r ${requirements}).`,
          )
        : e;
    },
  );

  if (code !== 0) {
    rmSync(work, { recursive: true, force: true });
    rmSync(marker, { force: true });
    const missing = output
      .map((l) => /No module named '([^']+)'/.exec(l)?.[1])
      .find(Boolean);

    if (missing)
      throw new Error(
        managed
          ? `MERFISHEYES's own Python lacks "${missing}" — delete ${dirname(python)}/.. and try again to rebuild it.`
          : `The Python at "${python}" (MERFISHEYES_PYTHON) lacks the pipeline's packages (no module "${missing}"). Unset MERFISHEYES_PYTHON to let MERFISHEYES set up its own, or install them there: pip install -r ${requirements}.`,
      );
    const reason = output
      .filter((l) => /^error: /.test(l))
      .map((l) => l.slice(7))
      .join(" ");

    throw new Error(
      `Could not convert ${abs}: ${reason || EXIT[code] || `prepare.py stopped (exit ${code})`}.\n${output.slice(-8).join("\n")}`,
    );
  }
  rmSync(marker, { force: true });
  writeFileSync(
    join(work, ".merfisheyes-source.json"),
    JSON.stringify(
      {
        input: abs,
        fingerprint: print,
        options,
        converted_at: new Date().toISOString(),
      },
      null,
      2,
    ),
  );
  rmSync(folder, { recursive: true, force: true });
  renameSync(work, folder);

  return { folder, reused: false, seconds: (Date.now() - started) / 1000 };
}
