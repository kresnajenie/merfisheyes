# Processing pipeline output spec (browser JS ≡ server Python)

Two pipelines produce the chunked S3 layout the viewer reads: the browser
("Upload & Save": `lib/utils/GeneChunkProcessor.ts`, `lib/utils/SingleMoleculeProcessor.ts`,
the adapters + `lib/workers/standardized-dataset.worker.ts`) and the server
(`scripts/process_spatial_data.py`, `scripts/process_single_molecule.py` run by
`worker/entrypoint.py`). **The same input must yield the same output from both.**
The Python scripts are the reference; where they were wrong the fix is applied to
both. Consumers (`lib/adapters/ChunkedDataAdapter.ts`, `lib/SingleMoleculeDataset.ts`)
stay tolerant of the legacy layouts already in S3.

## Single cell (`datasets/{id}/`)

| File | Format | Notes |
|---|---|---|
| `manifest.json` | JSON | `statistics.{total_cells,total_genes,spatial_dimensions,available_embeddings,cluster_count}`, `files.de_stats` (list), `files.expression_chunks`, `files.observation_columns`, `normalized: false`, `processing.spatial_scaling_factor: 1.0` (Python) |
| `coords/spatial.bin.gz`, `coords/{emb}.bin.gz` | `[u32 n][u32 dims][f32…]` gz | raw microns rounded to 2 dp; embeddings truncated to ≤3 dims |
| `expr/index.json` | JSON | `genes[{name, chunk_id, position_in_chunk}]`, `total_genes`, `num_chunks`, `chunk_size` |
| `expr/chunk_%05d.bin.gz` | binary gz | **one gene per chunk** (`chunk_size` = 1 for both pipelines); header `[u32 1][u32 numGenes][u32 chunkId][u32 totalCells]`, 24-byte gene table, per-gene sparse block `[u32 numCells][u32 nnz][u32 idx…][f32 val…]` |
| `obs/metadata.json` | JSON | per column `{type: "categorical"\|"numerical", unique_values, format: "bin-v1"}` (legacy files have no `format` → JSON arrays) |
| `obs/{col}.bin.gz` | **binary obs format v1** gz | see `lib/utils/obs-binary.ts` / `encode_obs_binary`: `[u32 1][u32 n][u8 kind][u8 width][u16 0][u32 dictBytes][dict JSON, 4-byte padded][codes \| f32]`; categorical dictionary in **sorted label order**, missing → `""`; numerical f32, missing → NaN |
| `palettes/{col}.json` | JSON | categorical only; default colours (`lib/utils/color-palette.ts` = Python list) assigned in **sorted label order**; h5ad `uns[col_colors]` overrides by **stored category order** when one colour per category and not all identical |
| `de/{col}.bin.gz` | binary gz | every categorical column with ≤ 20 000 values; `[u32 1][u32 G][u32 C][u32 0]`, C×(u32 len + utf-8), u32 counts[C], f32 means[G·C], f32 pct[G·C] (`encodeDeStatsBuffer` ≡ `write_de_stats_binary`) |

Rules shared by both pipelines:
- **Gene filtering**: control probes / blanks / codewords are dropped in every format
  (`lib/utils/gene-filters.ts` ≡ `should_filter_gene` in both scripts; `tests/unit/gene-filter-parity.test.ts`).
- **Obs columns exported**: h5ad — every obs column except coordinate columns and
  all-empty ones; Xenium / MERSCOPE — every metadata column except coordinate, id
  (`cell_id`, `id`, `barcode(s)`, `cell`, `EntityID`), gene-name columns and all-empty ones
  (`lib/adapters/table-obs-columns.ts` ≡ the Python loaders).
- **Type detection**: `lib/utils/column-type-detection.ts` ≡ `is_categorical` (name lists, ≥80 % unique → numerical, >50 % floats → numerical).
- **Missing categorical label**: `""` (browser used to write `"Unknown"`).
- **Memory**: the Python script streams genes chunk by chunk and reads only X/obs/obsm/var/uns-colours from the h5ad (see `feat/worker-memory-opt`).

## Single molecule (`datasets/{id}/`)

| File | Format | Notes |
|---|---|---|
| `manifest.json.gz` | JSON gz | `statistics.{total_molecules, unique_genes, spatial_dimensions}`, `genes.{unique_gene_names (sorted), molecule_counts}`, `has_unassigned`, `processing.{coordinate_range: "raw_rounded_2dp", scaling_factor: 1, …}` |
| `genes/{sanitized}.bin.gz` | f32 `[x,y,z,…]` gz | raw microns rounded 2 dp, 3 floats per molecule even in 2D; **every kept gene has an assigned file** (empty if all its molecules are unassigned) |
| `genes/{sanitized}_uuuuuuuuuu.bin.gz` | same | unassigned molecules; **only written when non-empty** |

Rules shared by both pipelines:
- `molecule_counts` is **always** present: `{assigned}` per gene, plus `unassigned`
  for every gene iff the input had a cell-id column.
- `total_molecules` = molecules of **kept** genes (post control-probe filter, post
  non-finite drop) — this is what becomes `Dataset.numCells` and the card count.
- `has_unassigned` = cell-id column present **and** ≥1 unassigned molecule among kept genes.
- Molecules with a missing gene or a non-finite coordinate are dropped (never
  written at 0,0 or as NaN).
- Cell-id default follows the preset: MERSCOPE → `cell_id` when present, Xenium/custom → none
  (the server auto-detect path passes `--cell-id-col` accordingly).
- The viewer shows the **user's title** (`Dataset.title`), not `manifest.name`.
- Server-processed SM rows get `manifestJson` copied onto the row at finalize.

## Rules pinned by the parity harness (`npm run parity`, docs/PARITY-TESTING.md)

- **Coordinate rounding** (spatial, embeddings, molecules): round the value in
  float64 as `np.round(x, 2)` does (round-half-even on `x*100`), then narrow to
  float32 — JS `round2` in `lib/utils/coordinates.ts`, Python `round_coordinates`
  / `_round_column`. Never narrow to float32 first. Parquet columns keep their
  own precision until rounding (float32 stays float32, float64 stays float64).
- **Expression values are float32** on both sides before anything is derived
  from them (DE means/pct accumulate float32 values in float64).
- **h5ad obs**: bool columns label as `"True"` / `"False"`; the fixed coordinate
  names (`center_x/y/z`, `centerX/Y/Z`, `x/y/z`, `X/Y/Z`, `x_centroid…`,
  `centroid_x…`) are never exported as obs; all-empty columns are skipped;
  sparse `X` (CSR/CSC) is supported by both.
- **Embeddings** are truncated to ≤ 3 dims and rounded like spatial.
- **Type detection** runs on per-cell values with missing (`null`, NaN, `""`,
  `"NaN"`) removed.
- **SM**: a blank / whitespace-only gene name is a missing gene (row dropped);
  an empty CSV coordinate cell is a missing coordinate (row dropped); an absent
  optional parquet z column means 2D; `processing.scaling_factor` is the
  integer `1`.
- **Obs dictionary JSON** is compact (`JSON.stringify` ≡ `json.dumps(…,
  separators=(",", ":"))`).

## Known remaining divergences (Tier 2, not in scope of `feat/pipeline-parity`)
- Two manifest `files.*` / `expr/index.json` schemas (JS writes both spellings of the
  keys consumers read; Python's `dataset_id` is `"local_dataset"`, overridden by the app).
- Embedding naming/dim-cap/rounding drift; Python drops MERSCOPE UMAP.
- Fingerprints: browser canonical content hash vs server raw-bytes hash (no cross-path dedupe).
- `statistics.cluster_count` semantics; `spatialDimensions` not persisted for server SM.
- Every server upload still presigns every chunk object on viewer open (perf).
