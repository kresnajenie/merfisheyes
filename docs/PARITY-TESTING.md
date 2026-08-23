# Pipeline parity testing

Goal: run **both** pipelines (browser JS, server Python) on the same inputs and
diff the outputs, per [docs/PIPELINE-SPEC.md](PIPELINE-SPEC.md). Branch:
`feat/pipeline-parity` (rebased onto `develop` after memopt merged).

## Run it

```bash
npm run parity                       # generate inputs → python → node → diff; exit 1 on any difference
python3 scripts/parity/run.py --only h5ad_dense sm_csv_merscope   # a subset
python3 scripts/parity/run.py --skip-generate                     # reuse .parity/inputs
python3 scripts/parity/diff.py .parity/py/<case> .parity/js/<case> --kind sc|sm
PARITY_DEBUG=1 npx vite-node -c vitest.config.ts scripts/parity/run-js.mts .parity/inputs .parity/js <case>
```

Needs: Node (vite-node ships with vitest), Python with anndata/h5py/scipy/
pandas/pyarrow. Everything lands under `.parity/` (gitignored):
`inputs/` (+ `cases.json`), `py/<case>`, `js/<case>`, `logs/`.

## Pieces (`scripts/parity/`)

- `generate-inputs.py` — 9 deterministic cases (seeded, ~1.4 MB total):
  h5ad dense f64 X / CSR f32 X with 3D spatial / coordinates only in obs;
  Xenium folder (`cells.csv`, real 3-column `features.tsv`, `matrix.mtx.gz`);
  MERSCOPE folder (`cell_metadata.csv`, `cell_by_gene.csv`); SM parquet
  (Xenium 3D, Xenium 2D, MERSCOPE with `cell_id`) and SM CSV (MERSCOPE with
  `cell_id`). Every case carries control-probe genes, categorical obs with
  missing labels, numerical/bool obs, an all-empty column, `uns[<col>_colors]`,
  a 50-D `X_pca`; SM inputs carry NaN/Inf coordinates, blank gene names, a gene
  whose molecules are all unassigned.
- `run.py` — orchestrates; invokes the Python scripts the way
  `worker/entrypoint.py` does (`--dataset-type` + `--cell-id-col` for SM).
- `run-js.mts` — the browser pipeline headlessly in Node via **vite-node**
  (the lib code imports ESM-only packages, so tsx/CJS won't do). Single cell
  goes through `workerApi.parseH5ad/parseXenium/parseMerscope` (exported from
  `lib/workers/standardized-dataset.worker.ts`) → `StandardizedDataset.
  fromSerializedData` → `GeneChunkProcessor` → `createManifest` /
  `prepareFilesForUpload` (extracted from the modal into
  `lib/utils/sc-upload-files.ts`). Single molecule: `SingleMoleculeDataset.
  fromParquet/fromCSV` → `SingleMoleculeProcessor`. Shims: a minimal
  `FileReader` for PapaParse; unique `File.name` per case (h5wasm's virtual FS).
- `diff.py` — decompresses `.gz`; JSON compared structurally after normalising
  the documented expected differences (`manifest` name / created_at /
  dataset_id / version / processing.*, JS-only `files.*` keys, JS per-gene
  stats in `expr/index.json`, `statistics.cluster_count`); expr chunks compared
  per gene block (gene-table `uncompressedSize` tolerated); `de/*.bin`
  decoded and compared per value; everything else byte-for-byte.

## Status

`npm run parity` → **9/9 cases identical** (2026-08-23). Getting there fixed
these real divergences (browser side unless noted):

- gene filter not applied for h5ad and MERSCOPE; Xenium `features.tsv`
  header sniff dropped the first gene of a real (headerless) file
- all-empty obs columns exported; obs coordinate columns (`center_x` …)
  exported when spatial came from obs; bool obs labelled `0/1` (now
  `True`/`False` like pandas); embeddings truncated to 2 dims (now ≤ 3) and
  not rounded; 3-D spatial written as 2-D; **CSR/CSC h5ad unsupported**
  (`X` group) — now expanded to the dense float32 matrix; float64 X kept as
  float64 (DE means differed by 1 ulp; now float32 like the server)
- Xenium/MERSCOPE column typing ran on the *unique* values (unique ratio
  always 1 → integer columns came out numerical); `""` now counts as missing
- coordinates rounded with `Math.round` on float32-truncated values: both
  pipelines now round the native value in float64, round-half-even
  (`round2` ≡ `np.round(x, 2)`), then narrow to float32
- SM: empty-string gene names kept (Python kept them too — both drop now);
  empty CSV coordinate cells became `0` (dropped now); **2-D parquet produced
  zero molecules** (absent optional z came back as an empty column)
- Python: obs dictionary JSON used `", "` separators (now compact like
  `JSON.stringify`); unused `""` gene category survived into the manifest;
  `scaling_factor` written as `1.0` not `1`.

## Not covered yet / follow-ups

- Xenium `analysis/clustering/*.csv` (both loaders read it; column naming of
  the resulting obs column differs) and MERSCOPE `cell_categories.csv` (only
  the browser reads it) are not generated — add cases when those rules are
  settled.
- MERSCOPE UMAP (Python drops it), embedding name casing, `cluster_count`
  semantics: Tier-2 per the spec, normalised away or not generated.
- `tests/fixtures/obs-binary/*.bin.gz` goldens were written by the old Python
  (`", "` dictionary separators); the TS decoder accepts both, so they still
  pass — regenerate when convenient.
- CI: wire `npm run parity` into a workflow (needs the Python deps).
