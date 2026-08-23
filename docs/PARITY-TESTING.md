# Pipeline parity testing (to be built/run on the Ubuntu server)

Goal: a harness that runs **both** pipelines on the same inputs and diffs the
outputs, per [docs/PIPELINE-SPEC.md](PIPELINE-SPEC.md). Branch: `feat/pipeline-parity`
(stacked on `feat/worker-memory-opt`; PR to `develop` after memopt merges).

Start a Claude Code session in the repo and say:
> "Read docs/PARITY-TESTING.md and build the parity harness; tier 1 parity code is
> already on this branch."

## What is done on this branch
- SM: counts always written (both), post-filter `total_molecules`, all-unassigned
  gene no longer aborts the browser upload, server auto-detect passes cell-id,
  `has_unassigned` semantics, non-finite molecules dropped, z optional in parquet,
  viewer shows the user's title, `manifestJson` stored for server SM.
- SC: binary obs format (both writers + reader, JSON fallback), 1 gene/chunk
  (both), browser writes `de/*.bin.gz` for all categorical columns, Xenium/MERSCOPE
  export every metadata column in the browser, Python filters control probes,
  palettes sorted + `uns` colours honoured by stored order (both), `""` null label,
  Python sparse-block header fix.
- Tests: `tests/unit/obs-binary.test.ts` (codec round trips + Python-written
  goldens in `tests/fixtures/obs-binary/`), `tests/unit/gene-filter-parity.test.ts`.

## Harness to build (`scripts/parity/`)
1. `generate-inputs.py` — small deterministic inputs under `tests/fixtures/parity/`
   (or generated at run time; keep them < 5 MB each):
   - h5ad: dense f64 X with `Blank-*`/`NegControl*` genes, categorical obs with NaN,
     numerical obs, `uns[col_colors]`, `obsm` X_spatial + X_umap(2D) + X_pca(50D);
     sparse CSR variant; obs-coordinate-fallback variant.
   - Xenium folder: `cells.csv(.gz)` with centroid + metadata columns + `cell_id`,
     `cell_feature_matrix/` (mtx/h5 as the Python loader expects), `analysis/clustering`.
   - MERSCOPE folder: `cell_metadata.csv`, `cell_by_gene.csv`, optional UMAP csv.
   - SM: parquet (Xenium columns, with/without `cell_id`, 2D/3D), CSV (MERSCOPE
     columns with `cell_id` incl. -1, NaN coords, a gene whose molecules are all
     unassigned).
2. `run-python.sh` — `python3 scripts/process_spatial_data.py <in> <out>` /
   `process_single_molecule.py … --dataset-type … --cell-id-col …` (mirror what
   `worker/entrypoint.py` passes: `--chunk-size 1` implicit now).
3. `run-js.ts` (Node, `npx tsx`) — drive the browser pipeline headlessly:
   - SC: `H5adAdapter` needs h5wasm (use its Node build) and `File` objects
     (`new File([buffer], name)` from `buffer`/`undici` in Node ≥ 20); then build a
     `StandardizedDataset` the way `standardized-dataset.worker.ts` does, and call
     `GeneChunkProcessor.processGenes/processCoordinates/processObservations/
     processPalettes/processDeStats` + the modal's `createManifest` (extract it to a
     module if needed). Xenium/MERSCOPE: `XeniumAdapter`/`MerscopeAdapter` with File
     objects + `loadAllClusters()`.
   - SM: `SingleMoleculeDataset.fromParquet/fromCSV` (hyparquet works in Node via
     `hyparquet/src/node.js`; PapaParse needs a File/Blob) then
     `SingleMoleculeProcessor.createManifest` + gene file blobs.
   - `CompressionStream` exists in Node ≥ 18; `DecompressionStream` too.
4. `diff.py` — compare the two output trees: decompress `.gz`, normalise
   `manifest.name`, `created_at`, `dataset_id`, `processing.created_by`, gene-table
   `uncompressedSize`; compare JSON structurally and binaries byte-for-byte; report
   per-file PASS/FAIL with the first differing offset / JSON path.
5. `npm run parity` → all of the above; CI later.

## Acceptance (what must be identical)
coords/*.bin.gz · expr/index.json (genes order + chunk ids) · expr/chunk_*.bin.gz ·
obs/*.bin.gz · obs/metadata.json (type/unique_values/format) · palettes/*.json ·
de/*.bin.gz · manifest statistics + files.de_stats · SM manifest (statistics,
genes.unique_gene_names, molecule_counts, has_unassigned) · genes/*.bin.gz set + bytes.

## Known expected differences to normalise
`manifest.name` (title vs file stem), `created_at`, `dataset_id`, `version` ("2.0" vs
"1.0"), JS-only `files.expression{}`/`observations`/`palettes:[]`, JS `expr/index.json`
extra per-gene stats, gene-table `uncompressedSize`, embedding names (case) and
MERSCOPE UMAP (Python drops it) — see the Tier-2 list in the spec.
