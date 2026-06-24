# Testing — MERFISHEYES

This repo has three layers of automated testing:

| Layer | Tool | What it covers | Where |
| --- | --- | --- | --- |
| Unit | Vitest | Pure helpers (coordinates, cluster utils) | `tests/unit/` |
| E2E | Playwright (TypeScript) | Real browser: local viewing + upload/save across every format | `e2e/` |
| Performance | Playwright + baseline scripts | Load / gene-selection / upload timings vs a stored baseline | `e2e/perf/`, `perf/`, `scripts/perf/` |

The E2E + performance pipeline loads **fixed test datasets**, reports clear
pass/fail with the exact failing step, and gates on performance regressions
against a committed baseline.

---

## Quick start

```bash
npm install
npm run e2e:install         # one-time: download the Playwright browser
npm run test:e2e:quick      # fast smoke suite over the tiny datasets
```

The suite launches its own dev server on **port 3100** (with the E2E hooks
enabled) so it never collides with a dev server you have running on 3000.

---

## Commands

| Command | Purpose |
| --- | --- |
| `npm test` | Unit tests (Vitest) |
| `npm run test:e2e:quick` | Fast `@quick` suite — local viewing, tiny datasets, headless |
| `npm run test:e2e` | Full E2E suite |
| `npm run test:e2e:headed` | Run with a visible browser (debugging) |
| `npm run test:e2e:ci` | Headless quick suite with CI reporters (`CI=1`) |
| `npm run test:upload` | Upload→save→reload suite (needs Docker infra, see below) |
| `npm run test:perf` | Performance suite — writes `perf/results/latest.json` |
| `npm run test:perf:compare` | Compare latest run vs baseline; non-zero exit on regression |
| `npm run test:perf:update` | Bless the latest run as the new baseline |
| `npm run test:infra:up` / `:down` | Start / stop Postgres + MinIO for upload tests |
| `npm run test:data:make-tiny` | Regenerate the committed tiny fixtures |
| `npm run test:data:fetch` | Materialize medium/large datasets on demand |

Useful env vars: `E2E_HEADED=1`, `E2E_PORT=3100`, `E2E_REUSE=1` (reuse a server
you started with `NEXT_PUBLIC_E2E=1`), `E2E_FIREFOX=1`, `PERF_SAMPLES=3`,
`PERF_FAIL_PCT=20`, `PERF_WARN_PCT=10`.

---

## The dataset databank

All test datasets are described in **`tests/data/manifest.json`** — the single
source of truth for each dataset's path, format, and *expected* counts / genes /
fields. Tests assert against these expected values.

- **Tiny datasets are committed** (`tests/data/tiny/…`, ~330 KB total) and cover
  every supported format: `h5ad`, `xenium`, `merscope`, `chunked`, and
  `single-molecule` (parquet). These drive the `@quick` CI suite.
- **Medium / large datasets are not committed.** They are produced on demand
  into the git-ignored `tests/data/generated/` directory, either generated
  (`source.kind: "generate"`) or downloaded with sha256 verification
  (`source.kind: "download"`).

Regenerate the committed tiny fixtures (deterministic, seeded):

```bash
npm run test:data:make-tiny
```

Materialize larger datasets before running perf on them:

```bash
npm run test:data:fetch -- h5ad-medium     # one dataset
npm run test:data:fetch -- --tier medium   # a whole tier
```

The generator (`scripts/testdata/generate.py`) can also be called directly:

```bash
python scripts/testdata/generate.py --format all --cells 300 --genes 40 --out tests/data/tiny
```

### Format notes (important)

- **Raw `.h5ad` fixtures are written with a _dense_ `X` matrix.** The app's local
  `H5adAdapter` (`lib/adapters/H5adAdapter.ts`) reads gene expression only from a
  dense `X` — it has no CSR/CSC path — so a *sparse* local h5ad renders spatially
  but fails gene selection with `No adapter available for gene expression data
  access`. (This only affects locally-loaded raw h5ad; uploaded/chunked datasets
  use `ChunkedDataAdapter`, which is fine.) Keep h5ad fixtures dense.
- **Large single-cell uses the `chunked` format, not raw h5ad.** A dense 500k-cell
  h5ad is a ~1 GB in-memory matrix that overwhelms the browser — exactly the case
  the app's chunked/pre-processed pipeline exists for. `chunked-large` is generated
  via `--format chunked` (which builds an intermediate h5ad and runs
  `process_spatial_data.py`), loads lazily, and supports gene selection.

### Adding a dataset

1. Add a generator/download entry (or commit a small fixture) and a
   corresponding object in `tests/data/manifest.json` with its `expected` block.
2. That's it — the viewer and perf suites iterate the manifest automatically.

---

## How render timing works (app hooks)

The app exposes a tiny, production-safe instrumentation object,
`window.__merfish` (see `lib/utils/test-hooks.ts`), **only when
`NEXT_PUBLIC_E2E=1`**. It records `performance.now()` timestamps for scene-ready
and render-complete, plus authoritative dataset stats (cell/molecule/gene
counts, dimensions). This makes timing deterministic instead of sleep-based and
lets tests verify counts exactly. With the hook disabled, tests fall back to
canvas presence (less precise). A few `data-testid` attributes
(`sc-scene-canvas`, `sm-scene-canvas`, `upload-progress`, `selected-gene-badge`)
provide stable selectors.

---

## Local viewing tests

`e2e/tests/viewer/local-viewer.spec.ts` runs one test per tiny dataset and, in
labelled `test.step()`s:

- loads the dataset locally and measures **time-to-initial-render** and
  **time-to-interactive**;
- verifies cell/molecule/gene counts and dimensions against the manifest;
- verifies the canvas actually rendered (non-blank);
- selects 2–3 known genes and measures **gene-selection latency**;
- exercises **drag-and-drop** loading for single-file formats.

---

## Upload & save tests

`e2e/tests/upload/upload-save.spec.ts` (tag `@upload`) drives the real
upload→save→reload flow for the single-cell and single-molecule paths. These
need a database and an S3 endpoint, provided hermetically by Docker:

```bash
npm run test:infra:up      # Postgres (:5433) + MinIO (:9000, console :9001)
npm run test:upload        # sets E2E_UPLOAD=1; loads .env.test automatically
npm run test:infra:down    # tear down + remove volumes
```

`docker-compose.test.yml` creates the `merfisheyes-test` bucket; `.env.test`
points the app at the local Postgres and MinIO (via the `AWS_S3_ENDPOINT`
support added to `lib/s3.ts`). `e2e/global-setup.ts` applies the Prisma schema
and ensures the bucket before the tests run. The tests verify the saved dataset
is `COMPLETE` via the API and that reloading it from S3 matches the expected
metadata and counts.

> Auth is **not** required — the upload API routes are open; only DB + S3 are
> needed.

---

## Performance regression

```bash
npm run test:perf            # measure (median of PERF_SAMPLES runs) -> perf/results/latest.json
npm run test:perf:compare    # diff vs perf/baselines/<platform>.baseline.json
npm run test:perf:update     # accept current timings as the new baseline (commit it)
```

`compare` prints a table of **absolute ms + % change** per metric and exits
non-zero if any metric regressed beyond the fail threshold (default **20%**;
10–20% is a warning). Baselines are namespaced per `platform` (and browser) so
timings from different machines are never compared. Commit
`perf/baselines/*.baseline.json`; `perf/results/` is git-ignored.

Perf is dev-server-based by default (a warmup load absorbs first-compile cost).
For paper-grade numbers, run against a production build locally.

---

## Failure reporting

- Every assertion includes **dataset, format, step, expected, actual** (and
  timing for perf) — see `e2e/helpers/assert.ts`.
- On failure Playwright captures a **screenshot, video, and trace**
  (`playwright-report/`, `test-results/`). Open a trace with
  `npx playwright show-trace <trace.zip>`.
- Tests are annotated with their dataset/format so the HTML report groups them.

---

## CI

`.github/workflows/e2e.yml`:

- **e2e-quick** — runs the `@quick` suite on every PR (tiny committed data, no
  Docker).
- **e2e-upload** — spins up Postgres + MinIO via compose and runs the `@upload`
  suite on every PR.
- **perf** — nightly (and manual `workflow_dispatch`) performance run +
  baseline comparison; uploads results and the HTML report as artifacts.

The existing `test.yml` (unit tests) is unchanged.

---

## Legacy benchmarks

The Python scripts under `scripts/benchmark/` and the chart generators in
`benchmarks/` are retained for paper figures. They are **not** part of this
pipeline and are not maintained as tests.
