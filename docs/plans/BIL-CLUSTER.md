# Plan: the same pipeline, run on the cluster

[ONE-PIPELINE.md](ONE-PIPELINE.md) made one script, `scripts/prepare.py`, that
a person and the server both run to turn raw data into one SpatialData store.
This plan adds the **third runner: the PSC cluster**, and the things BIL needs
that a single upload never does — many sections in one store, a fleet of 87
datasets, resume, and molecule counts that reach the billions.

It does not restate ONE-PIPELINE's decisions; it depends on them. Judged by
[docs/PHILOSOPHY.md](../PHILOSOPHY.md).

This file is the brief for every session working on it: what was measured on
2026-10-23, the decisions, the steps with their gates. **Keep it current: a
step that changes a fact or a decision edits this file in the same commit.**

---

## 1. What is true today (2026-10-23, `feat/one-pipeline-server` at `dcfdc3c`)

### 1.1 The entry point exists and the server already runs it

```
python scripts/prepare.py INPUT OUT [--name] [--no-images] [--keep-z]
                          [--cluster] [--labels CSV] [--molecule-columns] [--workers]
```

Input recognised by content; output one store with `tables/table` (CSC `X`,
DE stats in `varm` since `feat/de-in-table`), `points/molecules`,
`shapes/`, `images/`, plus `mapping.json` and `stats.json`. Exit codes 0/2/3/4.
`scripts/requirements.txt` is **fully pinned (140 packages), locked against
`python:3.12-slim`**. ONE-PIPELINE steps 1–3 are done; 4 (docs) and 5 (retire
the chunked processors) are not.

**Measured there, and it matters here:** the new writer needs ~0.016–0.024 GB
per million molecules against the old processor's ~0.11. Beaumont set1 —
975 M molecules — went from 28:55 / **110.2 GB** (old, no tier fits) to
17:16 / **15.6 GB** (new).

### 1.2 What BIL needs that `prepare.py` does not do

| Need | Status |
|---|---|
| **Many sections in one store** | No `combine` step and no `--section` column. `prepare.py` takes one input. BIL items have up to **60** sections. |
| **Run MapMyCells** | `--labels CSV` takes MMC's *output*; running MMC (`scripts/map_my_cell.py`, 64 CPU / 512 G, mouse + human references) is outside the script. |
| **A cluster runner** | Nothing. The BIL pipeline is still `scripts/bil-scripts/`: 6 bash launchers, 10 `.sbatch`, 5 inline `--wrap`, no SLURM arrays, **no state or resume**, S3 sync logic in 4 shapes, hard-coded personal paths. |
| **Python ≥ 3.12** | BIL's venv is built on the `python/3.10.2` module (`bil-scripts/SETUP.md:5-29`). ⚠ **The pinned stack cannot run in `~/merfisheyes_env` as it stands** — verify in Step 0. |
| **Billion-scale molecule builds** | §1.4 below. |
| **Fleet-level publishing** | 87 BIL items currently expand to **446 catalog entries, 361 of them one-per-section**; `scripts/new-datasets/bil-completed*.json` is hand-maintained. |

### 1.3 The BIL pipeline's working state is scattered across unmerged branches

**Tier A — in the repo, not on the mainline. ✅ Done on `chore/bil-consolidate`.**
34 `bil-scripts/` files existed only on `pipeline/final-rerun` / `mmc-update` /
`s3-update` / `testing/full-pipeline` / `feature/spinal-rna-spots`
(`develop` carried 24; `testing/full-pipeline` 38). **29 restored**, each from
the branch holding its newest version: the `analysis-scripts/` directory,
`launch_mmc.sh`, `launch_xenium_pipeline.sh`, `launch_h5ad_sc_only.sh`,
`launch_sync_all.sh`, `s3_sync_sm_output.sbatch`, `collect_stats.sbatch`, the
three `count_*_per_year.sbatch` wrappers, the BIL scrapers
(`scrape_bil_merfish.py`, `scrape_bil_xenium.py`, `populate_bil_*.py`,
`merge_genes_into_bil_examples.py`), the `xenium-/csv-/ace-pie-/
samples_reprocess_border` lists, and **`scripts/combined_output_to_h5ad.py`** —
the converter without which MERSCOPE input cannot reach MapMyCells.

Deliberately **not** restored: six top-level copies of files that were later
moved into `analysis-scripts/` (`all_datasets.csv`,
`collect_correlation_sections.py`, `compare_mmc_methods.py`,
`count_{cells,genes,xenium_cells}_per_year.py`) — byte-identical to their newer
copies, and the `.sbatch` wrappers already point at `analysis-scripts/`; and
`analysis-scripts/h5ad-samples.csv`, identical to the top-level list and
referenced by nothing.

**The three shared scripts. ✅ Resolved (user, 2026-10-23) — but not by "take
the newest".** In two of three, no single version is a superset, so the newest
file would have deleted working code:

| File | Resolution | Why |
|---|---|---|
| `map_my_cell.py` | **`mmc-update` + the `center_x` fallback from `s3-update`** → 722 lines | The tip (2026-07-22, 292 lines) is 15 days *newer* but was stripped in `9f4f981` and has **zero** correlation support (`correlat` ×0 vs ×30). Dual-method needs the branch. Of the three dual-method branches `mmc-update` (07-07) alone carries the **case-insensitive gene-symbol fallback** — without it a mouse panel against the human reference "would otherwise drop nearly every gene", which is exactly the macaque→human path. `s3-update` (07-09) is 2 days newer but dropped that fix and instead added the **`obs['center_x']` fallback for BIL SISSeg h5ad exports with no `obsm['X_spatial']`** — the h5ad sample family. Both are needed, so `mmc-update` is the base and the 14-line fallback is ported onto it. |
| `process_spatial_data.py` | **tip + `92671a3`** | The tip (10-22) is genuinely newest and far richer. But it never received the 09-04 BIL fix: a `ProcessPoolExecutor` crash on `pd.Categorical` with int16+ codes, and — worse — `is_categorical()`'s `series.dtype == "object"` test, which **silently no-ops under pandas ≥ 3.0** (plain strings infer as `StringDtype`), misclassifying every high-cardinality string column as numerical and encoding it all-NaN. The tip still had the buggy line at `:258`. |
| `combine_slices_v3.py` | **tip + `836a73a`** | Here "newest" would have been actively destructive: `pipeline/final-rerun` (09-09) is 5 days newer than the tip (09-04) and has the panel-mismatch fill, but contains **zero occurrences of `sample_offsets`** — the per-slice translation table §2.3 depends on to put molecules in the combined cells' frame. Tip is the base; the fill is ported. |

Both ports applied with `git apply --3way` cleanly; every symbol other scripts
import from `process_spatial_data` (`OBS_BINARY_FORMAT`, `fmt_elapsed`,
`peak_rss_gb`, `write_coordinate_binary`, `build_obs_palette`,
`encode_obs_binary`, `write_de_stats_binary`) still resolves.

**The `process_spatial_data.py` port is verified empirically** (2026-10-23, conda
`base`, pandas 3.0.1 — the version where the bug bites). pandas 3 infers plain
strings as `str`, not legacy `object`, so on a 5,000-row barcode column:

| Version | `is_categorical(barcodes)` | Consequence |
|---|---|---|
| tip, unfixed | **False** | classified numerical → encoded all-NaN |
| tip + `92671a3` | **True** | correct |

⚠ Still unrun: `pytest tests/python` (needs `pytest`, `jsonschema`, `tifffile`,
`spatialdata`) and `npm run parity` (needs `node_modules`). The `map_my_cell.py`
and `combine_slices_v3.py` merges are verified only by compile, symbol
resolution and feature presence — **not** by processing a dataset. Step 0's
control run on one BIL dataset is the first real exercise.

`main`'s pipeline still cannot run: `launch_pipeline.sh:195` hands a directory
to an h5ad-only `map_my_cell.py`; `launch_sm_pipeline.sh:190` passes named flags
to a positional-only sbatch. Fixed on `pipeline/final-rerun` (`41cc34a`,
`317deb4`) — **and the restored `map_my_cell.py` accepts a directory again**, so
the first of those two is now moot on this branch.

**Tier B — in no branch of this clone.** `align2.py` and `build_store.py`,
which built the 11.24 B-molecule `12_16_2025_all_sets` store that
[BUDGETS.md](../specs/BUDGETS.md) measures against
([molecule-pyramid.md:254,281](../molecule-pyramid.md)). Searched every ref and
all reachable history: absent. They are on the author's machine or the cluster.
**`scripts/bil-scripts/find_missing_scripts.sh` searches for them by name, by
`PointsWriter` use, and by parallel-spill shape** — run it over SSH on BIL
before Step 4 reimplements `build_store.py`, which is exactly the parallel-spill
wrapper §1.4 is missing.

**Tier C — built, documented, discarded.** The `labels` element kind
(`lib/spatialdata/labels.ts`, `LabelsLayer.ts`,
`scripts/testdata/make-tiny-labels.py`) is described as working and
conformance-passing in
[SECOND-IMPLEMENTER-2026-10-17.md](../specs/SECOND-IMPLEMENTER-2026-10-17.md)
and exists on no branch. **Decision (2026-10-23): leave it gone.** That document
is a dated record of an experiment and its references stay as written; the
backlog it produced is what matters, and that is already in
[LAYER-CONTRACT.md](../specs/LAYER-CONTRACT.md). Rebuild the kind only when a
store needs labels drawn.

**Tier D — stale references. ✅ Fixed.** `docs/data-requirements.md:228` named
`scripts/process_h5ad.py`, removed in `9f69d63`; it now names
`process_spatial_data.py`. No documented `scripts/*` path is dangling except
`make-tiny-labels.py`, which is Tier C's historical record.

Also known broken for BIL: `ace-low-cop.h5ad` fails the Zarr writer on
`uns['cell_polygons']` (ONE-PIPELINE §1; decision g leaves such entries out
with a warning — **verify it actually covers the BIL files**).

### 1.4 The points writer is proven at 11.24 B, but only by hand

`scripts/points_zarr.py` works out of core and has written 11,239,601,084
molecules in 8 levels — ~2 h, peak ~42 GB, 248 GB / 38,030 files
([molecule-pyramid.md:294-321](../molecule-pyramid.md)). The algorithm is not
in question. The CLI around it is:

| Gap | Detail |
|---|---|
| **No parallel spill** | The 32-worker spill lived in the missing `build_store.py`. `points_zarr.py` / `transcripts_to_points.py` / `legacy_to_points.py` contain no `multiprocessing`. `PointsWriter.spill()` is picklable and designed for it. |
| **Temp-file explosion** | `transcripts_to_points.py:210` tags per 2 M-row batch → 11 B rows ≈ **240,000 temp files**; pass B then concatenates ~5,620 small arrays per band. The real build tagged per worker. |
| **`band_um` not on any CLI** | Constructor-only, 800 µm. `transcripts_to_points.convert` and `legacy_to_points.convert` never pass it; a dense dataset in a short y-extent OOMs (band working set ≈ 70–90 B × N / bands) with no remedy but a source edit. **`scripts/testdata/make-medium-store.py:152` does pass it**, driving `PointsWriter` directly — it is the working template for the wrapper, still serial and still tagging per part (`:174`). |
| **Temp dir fixed** | `store.parent/.{name}.tmp`, ~150–240 GB at 11 B. No `--tmp`. |
| **No resume** | `transcripts_to_points.py:181` exits if the output exists; `PointsWriter.resume()` must be called by hand and loads a whole level into RAM. |
| **Validation unskippable** | `prepare_folder` always ends in `validate_store.validate()` → `pz.validate(sample=50_000)`, reading a full 4 M-row shard per sampled span. Hours on a big store. |

**Hard ceiling: one gene must fit in RAM** — pass C loads a whole gene's
level-0 rows, ≈ 56 B × its count. Slc17a7 at 589 M = **~33 GB**. For one BIL
dataset (~1 B molecules, 500-gene panel, top gene ~5 %) that is **~3 GB** —
inside the 512 G the pipeline already requests. ⚠ Decision D1 (one item = one
store, every section together) pushes a large item's total molecule count up,
so the fleet's **largest gene after combining** is the number Step 0 must
produce.

### 1.5 Three input families, each with its own sample list

| Family | Sample list | Input | Cells | Molecules |
|---|---|---|---|---|
| **MERSCOPE multi-slice** | `samples.csv`, `csv-samples.csv` (`sample_name,species,percentile`) | a BIL dir of N slice subdirs | `cell_metadata.csv` + `cell_by_gene.csv` per slice → **combined** by `combine_slices_v3.py` | `detected_transcripts` per slice, **not combined** |
| **H5AD** | `h5ad-samples.csv` (`sample_name,input_path`) | one dir | `cell_by_gene.h5ad` | `segmented_spot_table.csv`, columns `gene_names,x,y,z,cell_ids` (`process_h5ad_sm.sbatch:58-67`) |
| **Xenium** | `xenium-samples.csv`, `ace-pie-samples.csv` | a raw Xenium run folder | spatialdata-io's reader | `transcripts.parquet` |

`prepare.py` already recognises Xenium and MERSCOPE column conventions natively
and takes `--molecule-columns GENE,X,Y[,Z]` for the h5ad family's spot table.
The per-sample **species** (mouse/human reference for MapMyCells, macaque → the
human reference) and **mask percentile** travel in `csv-samples.csv` columns 2
and 3 — on `pipeline/final-rerun`, not on `main`.

### 1.6 Test data on this machine (measured 2026-10-27)

Under `~/data/merfisheyes-test-data/elvis-bil-data/bil-download/`. All three are
real BIL MERSCOPE items with **cells and molecules**, and all name their files
`cellpose_metadata.csv` / `cellpose-cell-by-gene.csv` /
`cellpose-detected_transcripts.csv` — **not** the canonical names, which is why
`combine_slices_v3.py`'s fuzzy matching exists and what any new runner must keep.

| Dataset | Slices | Transcripts | Reference artefacts | Use it for |
|---|---|---|---|---|
| `ace-dry-dip` | 1 (`region_0`, canonical names) | 1.2 GB total | `combined_output/`, `mmc_output_corr/` | the **fast loop** — whole pipeline in minutes |
| `ace-dip-tug` | 3, all complete | 6.4 + 17 + 17 GB | `combined_output/`, **`mmc_output/` and `mmc_output_corr/`** | the **real gate** — combine, dual-method MMC, SC+SM in one store |
| `ace-dip-tub` | 3, but `1231122263` has **no `cellpose_metadata.csv`** | 3.0 + 13 + 2.7 GB | none | the **failure case** — a slice the combiner must skip or fail loudly on |

Schemas (from `ace-dry-dip/region_0`):

- transcripts — `barcode_id, global_x, global_y, global_z, x, y, fov, gene,
  transcript_id`: textbook MERSCOPE, recognised natively, and **no cell
  column**, so molecules here are written unsplit (D3).
- cell metadata — an unnamed index column plus `fov, volume, center_x, center_y,
  min_x, max_x, min_y, max_y`: no `obsm/spatial`, so the `center_x/y` path
  matters — the fallback ported into `map_my_cell.py` from `s3-update` (§1.3).

`ace-dip-tug/mmc_output/` holds `mapping_output.csv`, `enriched_metadata.csv`,
`mapping_output.json`, `query.h5ad` and `check_celltype_*.png`; the restored
(pre-strip) `map_my_cell.py` writes exactly that set, while the tip's stripped
version writes only the CSV — independent confirmation that the restored version
is the one that produced this data.

Also here: `~/data/merfisheyes-test-data/ace-low-bag-points/` holds the **same
sample in both molecule formats** — `sm_output/` (legacy per-gene, including
`*_uuuuuuuuuu.bin.gz` unassigned files) and `sm_points/sm_points.zarr` (the
tiled points element) — a ready-made before/after for §2.2.

---

### 1.7 The pipeline run end to end on this machine (2026-10-27)

`tests/bil/run_dry_run.sh` runs the cluster's stages without SLURM — the
chain `launch_combine_mmc.sh` and `launch_process_sync.sh` build with
`sbatch --dependency=afterok`, chained by `set -e`. Measured on a 64-core
box, python 3.12 with `scripts/requirements.txt` pinned:

**`ace-dip-tug` — 3 slices, the shape the fleet needs**

| Stage | Time | Result |
|---|---|---|
| combine | 1 m 41 s | 407,859 cells, `sample_offsets.json`, masks p5–p80 |
| mask | 13 s | |
| annotate | — | reused `mmc_output/mapping_output.csv` (references are 14 GB) |
| prepare (cells) | 17 s | `tables/table`, CSC, MMC labels merged |
| **points (molecules)** | **8 m 28 s** | **433,023,021 molecules, 8 levels** |
| validate | 1 m 15 s | PASS |

**11 m 58 s, peak RSS 4.9 GB. 40.6 GB of CSV in → a 6.0 GB store of 2,251
files.** Peak memory is the number that matters for the cluster: half a
billion molecules at under 5 GB, because nothing is ever concatenated —
each slice streams into `PointsWriter`, which spills and builds the levels
once.

Per slice, every molecule landed on its own cells (100 % inside the cells'
box + 50 µm): `1232698492` 69.7 M molecules / 140,958 cells shifted
(+40, +2); `1232698496` 177.7 M / 120,378 (+11241, −3); `1232698499`
185.7 M / 146,523 (+129, −9614). 434,118,398 source rows − 1,095,377
blanks and control probes = 433,023,021 written, exactly.

**The other two families, through `prepare.py` directly** — neither has a
combine step (`launch_h5ad_pipeline.sh` and `launch_xenium_pipeline.sh` go
straight to processing), so `slices_to_points.py` does not apply to them:

| Input | Time | Result |
|---|---|---|
| Xenium `Human_Pancreas` | 26 s | 103,901 cells, 377 genes, 13.9 M molecules, both boundary sets, a 6-level morphology image |
| h5ad family `ace-low-cop` | 1 m 45 s | 151,391 cells, 500 genes, 73.9 M molecules (57.4 M in a cell, 16.5 M not) |

**`ace-dip-use` is the migration case.** On this machine it is only a
previous run's `combined_output` — no raw slices, no molecules, no
`sample_offsets.json`, which is what all 87 existing datasets look like. Its
cells converted (448,084 cells, DE stats on 4 obs columns) and the points
stage refused with a message naming the combine. That is the evidence
behind D5.

**What this run changed in the code**

- `combine_slices_v3.py` finds each slice's `*detected_transcripts*` with a
  third fuzzy pattern and records its path beside that slice's `(dx, dy)` in
  `sample_offsets.json`. The molecules are never combined into a file — the
  combine is the only record of where they are.
- `combine_slices_v3.py` skips the directories this pipeline writes
  (`combined_output`, `mmc_output*`, `sm_output`, `meyes_output`, `*.zarr`).
  Without that, a re-run over a folder holding its own output reads it back
  as another slice and **doubles every cell**, exit 0, silently — 815,718
  rows instead of 407,859 on `ace-dip-tug`.
- `scripts/slices_to_points.py` is new: the stage that reads the offsets and
  writes every slice's molecules into one element, each tagged with its
  `sample`.
- `transcripts_to_points.convert_many()` does the writing; `convert()` is a
  one-line wrapper on it.
- `gene_names` joined `MOLECULE_COLUMN_NAMES` in all three copies
  (`meyes_common.py`, `worker/entrypoint.py`,
  `lib/ingest/molecule-columns.ts`) — without it the h5ad family's spot
  table is not recognised as molecules at all. `cell_ids` joined
  `cell_column()`.
- `tests/bil/` is the test space: 12 tests on generated fixtures (~30 s, no
  data) plus the runner above.

### 1.9 ace-dip-vat: four slices, 646 M molecules (2026-10-28)

Run with `tests/bil/run_dry_run.sh`, output in
`~/data/merfisheyes-test-data/bil-data/merfish_data_processed/ace-dip-vat/`.

| | |
|---|---|
| in | 60 GB, 4 slices |
| cells | 474,986 |
| molecules | **646,439,087** of 648,044,571 source rows (1,605,484 blanks/controls dropped) |
| in no cell | 218,565,924 (33.8%) |
| store | 8.8 GB, 3,034 files, 8 levels |
| prepare + points + validate | 15 m 52 s, peak RSS 6.4 GB |
| MapMyCells (correlation) | 474,986 cells x 5,322 clusters, ~2 min |
| `tests/bil/verify_store.py` | every check passed |

Each of the four slices placed 100% of its molecules on its own cells.

**It found a real bug, which is why it was worth running.** 54,553 of those
474,986 cell ids appear in more than one slice — the instrument numbers cells
per slice — and `prepare_merscope.py` was the one place still treating the raw
id as an identity. `counts.loc[ids]` on a duplicated index multiplies the
repeated rows together: 474,986 cells against 584,092 rows of expression, and
AnnData refused them, *after* the combine and MapMyCells had already run. The
table is now indexed by `<_sample_id>_<id>`, the same identity
`map_my_cell.py` and `process_spatial_data.py` use, with the slice's own id
kept in `slice_cell_id`.

A label file is now only joined on a key that is unique on both sides.
MapMyCells' compound ids identify a cell, so they **join by id**; the artifact
mask is keyed by the per-slice id that two slices share, so it falls back to
the order the combine wrote it in, and says so. Pinned by
`tests/bil/test_duplicate_cell_ids.py`.

### 1.8 The cluster scripts for this workflow (2026-10-28)

`scripts/bil-scripts/` gained **Pipeline D**: the seven stages above as Slurm
jobs, chained with `--dependency=afterok` exactly as the existing launchers
chain theirs.

| File | Does |
|---|---|
| `env.sh` | the one place paths, the bucket and the interpreter live. `meyes_run` dispatches to a venv **or** `apptainer exec`, so step 1's container is a one-variable change (`MEYES_IMAGE`). |
| `prepare_store.sbatch` | combined cells → store, with every label CSV merged |
| `slices_to_points.sbatch` | every slice's molecules → that store |
| `validate_store.sbatch` | the gate: non-zero exit, nothing syncs |
| `s3_sync_store.sbatch` | the store → `s3://bil-psc-xarr/<sample>/` |
| `launch_spatialdata_pipeline.sh` | combine → mask → annotate → prepare → points → validate → sync |

**The bucket is `bil-psc-xarr`** (`MEYES_S3_BUCKET` in `env.sh`), separate from
the legacy `merfisheyes-bil/bil-psc-data2`. It needs public read and CORS, or
the viewer cannot open a store by URL — a public URL cannot be listed, so the
consolidated metadata is the only way it learns what is in the store
(`docs/specs/store-layout.md`).

Every stage was run off-cluster by executing the sbatch bodies directly, with
their failure paths: a missing `combined_output` (exit 3), a missing
`sample_offsets.json` (exit 3, naming the combine), a missing or empty store
(exit 3). `--dry-run` prints the whole chain without submitting, which is the
only way to check the dependency graph without a scheduler.

**The artifact mask is kept, not applied.** The legacy chain dropped the cells
`artifact_mask_p25.csv` flags; the new chain merges it as an `is_artifact` obs
column (101,982 of 407,859 on `ace-dip-tug`) and keeps every cell, so the
filtering is the viewer's to do and is reversible. `--labels` is repeatable,
which is what carries both the mask and MapMyCells' labels. **Say so if cells
should be dropped instead** — that is a different decision, not an oversight.

**MapMyCells now measured here**: 407,859 cells against 5,322 reference
clusters in **1 m 07 s**, 5.6 GB peak, correlation method, using the 14 GB
reference on this machine. The `hierarchical` method additionally needs
`cell_type_mapper`, which the image must carry. The labels it produced merged
onto the right cells — verified by compound id, all matched, all identical.

**`stats.json` is refreshed by the points stage**, which otherwise left
`"molecules": 0` on a store holding hundreds of millions — it is the file a
catalog reads.

**Corrected from the 2026-10-26 reading.** The MapMyCells join was reported
as silently mislabelling 84.7 % of cells. That was a test against the wrong
table: MMC reads `combined_output`, not a raw region. Verified by compound
id on both `ace-dip-tug` and `ace-dip-use` — every cell matched, every label
identical. `cell_by_gene.csv`, `cell_metadata.csv` and `mapping_output.csv`
are all in one order, so the positional merge in `prepare_folder.merge_labels`
is correct. One narrow bug remains: `process_spatial_data.reorder_df_to_reference`
strips the compound id at the **first** underscore, so a `_sample_id` that
contains one (`region_0`) yields `0_101…` and `--reorder` raises. It is off
by default, so the cluster has never hit it.


## 2. The format change: how gene storage differs

The thing to understand before writing any code.

### 2.1 Cells — expression

| | Old (`expr/`) | New (table `X`) |
|---|---|---|
| Layout | `expr/index.json` + one gzipped object **per gene** (`chunk_size: 1`) | one AnnData table; `X` CSC, or `layers/X_csc` when `X` is CSR |
| One gene costs | a whole-object `GET`, ungzip, rebuild a dense `Float32Array(n_cells)`, return a boxed JS `number[]` (`ChunkedDataAdapter.ts:784-838`) | read `indptr` once, then a **slice** of `indices`+`data` → ranged reads of overlapping chunks; **17 reads/gene on the ten-section store** (`SpatialDataTableAdapter.ts:438-454`) |
| Scaling | one S3 key per gene, each presigned and fetched whole | one set of sharded arrays, any panel size |
| Failure mode | — | a CSR `X` with no CSC copy scans the whole matrix per gene — `validate_store.py:178` makes it a FAIL |

### 2.2 Molecules — the real change

| | Old (`genes/<gene>.bin.gz`) | New (`points/molecules`) |
|---|---|---|
| Layout | one gzipped `float32[x,y,z,…]` per gene, 12 B/molecule even in 2-D; unassigned in a `_uuuuuuuuuu` file | 8 levels (level 0 all, each level a uniform **quarter**), each in **two row orders** (`tile_major`, `feature_major`) + a row-number index |
| One gene costs | **the whole gene, always** — Slc17a7 at 589 M is a **7 GB** download for any view | index slices, then contiguous row ranges for the tiles in view, at the level whose count fits the budget |
| Bytes scale with | the dataset | **the screen** — Auto = 8 M molecules / 12 MB per view (`PointsLayer.ts:32-41`) |
| Zoom out | nothing; one resolution exists | a coarser uniform sample, density honest at every zoom |
| Measured | — | a **whole-atlas view of three genes on the 11.24 B store reads 29 MB** |
| Sections | one dataset per section, linked by `mapping.json` + right-click | one element, a `section` column, all sections drawn together |
| Build cost | ~0.11 GB per million molecules | ~0.016–0.024 GB per million |

**What it buys:** a dataset physically unopenable in the old format (a 1 B-molecule
gene is a 12 GB file) becomes interactive, with reads bounded by the view.

**What it costs:**

- No `_uuuuuuuuuu` concept — unassigned becomes a column value or a `cell_id`
  column (**D3**).
- Coordinates are **float32** in the spill record and the Zarr columns: fine at
  µm scale, lossy for raw nanometres over a large extent.
- The extent and gene set are needed up front → **two passes** over the transcripts.
- `spatialdata` 0.8 **cannot read** the points element (`KeyError: '0.3-zarr'`);
  `read_zarr(..., on_bad_files="warn")` reads everything else and skips the
  molecules (ONE-PIPELINE step 1).

### 2.3 The asymmetry this reprocessing exists to fix

**Cells are combined across slices; molecules are not.**

- `combine_slices_v3.py` merges every slice into one `cell_metadata.csv` +
  `cell_by_gene.csv`, grid-shifting each slice by a pure translation, tagging
  each cell with `_sample_id` (`:434`, the slice's top-level folder name), and
  writing the exact per-slice shift to **`sample_offsets.json`** (`:548`):
  `{sample_id, dx, dy, orig_bbox, grid_origin, num_rows}`, documented as
  "original = combined − (dx, dy); pure per-slice translation".
- `process_single_molecule.py` in directory mode writes **one output folder per
  slice** (`{out}/{sample_id}/manifest.json.gz` + `genes/`) plus a
  `mapping.json` with `linkColumn: "_sample_id"` pointing each slice at its own
  S3 URL. That is where the 361 per-section catalog entries come from, and why
  the only way to see a slice's molecules today is right-click → split panel.

**The join is already exact.** `process_single_molecule.py:202` states its
`sample_id` "is the top-level directory name relative to target (matches
`_sample_id` in `combine_slices_v3.py`)" — both derive it the same way, by
construction.

So putting molecules into the combined cells' frame is a **lookup, not a
search**: for each slice, add that slice's `(dx, dy)` from `sample_offsets.json`
and write every slice into one `points/molecules` with a `section` column equal
to `_sample_id`. No cross-correlation. `align2.py`'s search was needed only for
the 12_16_2025 retrofit, where both folders had been normalized to [-1,1]
**without their centre** ([molecule-pyramid.md:228](../molecule-pyramid.md)) —
provenance was already lost. It is not needed for a dataset the pipeline
combines itself.

⚠ **Verify once in Step 0** that a slice's transcripts and its
`cell_metadata.csv` share one frame before the shift (for MERSCOPE both are
microns in the mosaic frame, so they should), and that the h5ad family's
`segmented_spot_table.csv` does too.

---

## 3. Decisions

**Proposed (2026-10-23 — to confirm):**

1. **No new CLI.** `prepare.py` is the entry point on the cluster too. What
   PSC needs is *orchestration and packaging*, not another interface.
2. **One store per BIL item**, every section in one `points/molecules` with a
   `section` column — as `12_16_2025_all_sets` does with `set`.
3. **One artifact everywhere**: the worker image, run under Apptainer on PSC
   and on Fargate as today. The pinned `scripts/requirements.txt` is the
   contract; the cluster venv is retired.
4. **A per-dataset state file is the unit of work** — stage, input
   fingerprints, outputs, timings, exit status. Replaces `sync-remaining.csv`.
5. **A SLURM array over a sample manifest** replaces the six bash launchers.
6. **`validate_store.py` is the gate**: non-zero exit, nothing syncs.
7. **No backfill** (ONE-PIPELINE decision c): existing datasets stay readable.

**Open — needed before Step 3:**

| # | Question | Why it blocks |
|---|---|---|
| ~~D1~~ | ~~Is a BIL "dataset" one item or one section?~~ | **Answered (user, 2026-10-23): one item, one store** — the combined cells and every section's molecules together. That is the point of the reprocessing (§2.3). |
| ~~D2~~ | ~~Does MapMyCells run inside the pipeline?~~ | **Answered (user, 2026-10-27): yes, inside.** So `annotate` is a CLI stage, the image carries `cell_type_mapper`, and the references are **bind-mounted, not baked in** (14 GB: `/bil/data/meyes/mapmycells-reference` on the cluster, `~/data/merfisheyes-test-data/mapmycells-reference` locally). This diverges from the app's path, which has no server-side cell typing and takes a user CSV — `prepare.py --labels` is the seam: the cluster's `annotate` stage produces the CSV that `prepare.py` consumes, so neither side grows a second implementation. |
| ~~D3~~ | ~~How are unassigned molecules represented?~~ | **Answered by the code (`572ba19`, on develop):** a points element may carry an `unassigned` column copied into `feature_major`, with one shared rule (`meyes_common.unassigned_mask`: any number ≤ 0, NaN, blank, or the words UNASSIGNED / none / nan / na / null). A cell column is **never guessed**; MERSCOPE's `cell_id` is used when present, anything else names one through `--molecule-columns gene=G,x=X,y=Y[,z=Z],cell=C`. **Measured again 2026-10-27 (§1.7):** the multi-slice sets *do* carry `cell_id`, written `-1` for a molecule in no cell — 140.7 M of 433.0 M on `ace-dip-tug`. The single-region `ace-dry-dip` export does not. The BIL h5ad family's `segmented_spot_table.csv` names it `cell_ids`, now recognised. |
| ~~D4~~ | ~~Apptainer on PSC; S3 from compute nodes?~~ | **Decided (user, 2026-10-27): one workflow — the container does the processing *and* the S3 sync.** So `sync` is a CLI stage run inside Apptainer on the compute node, not a separate login-node job, and today's four copies of `aws s3 sync` collapse into one. **Partly answered by the scripts (2026-10-27):** `launch_sm_pipeline.sh` already runs `aws s3 sync` in a `--wrap` job on `--partition=compute` after `module load aws-cli`, so compute nodes **do** have the CLI and outbound S3, and sync is already an in-pipeline stage. ⚠ **Still to verify on the cluster** (Step 0): that Apptainer exists, and how credentials reach the container without the `~/.aws/config` race that was fixed three times (`2c15173`). If compute nodes turn out to have no egress, `sync` becomes a dependent login-node job and this decision narrows to "one implementation, two schedulers". |
| D5 | Re-run all 87, or only new datasets? | Decides whether §1.4's gaps must be closed before Step 5. **Evidence (2026-10-27): re-run.** No existing `combined_output/` on this machine holds `sample_offsets.json` — it did not exist when they were made — and without the per-slice shifts the molecules cannot be placed in the combined cells' frame. The cell half of a re-run is cheap: 1 m 41 s for `ace-dip-tug`'s three slices. |
| D6 | Does ONE-PIPELINE step 5 (retiring the chunked processors) wait for the BIL re-run? | `scripts/bil-scripts/` still calls them. |

---

## 4. Steps

Each step is one branch, one gate. **A step is not done until its gate passes
and this file records what changed.**

### Step 0 — recover the ground truth (`chore/bil-consolidate`, then `chore/bil-census`)

No product code.

- ✅ **Tier A additive files restored, Tier D fixed, Tier C decided** (§1.3).
- **Settle the three shared scripts** (§1.3): decide whether the stripped
  `map_my_cell.py` regains the dual-method MMC work, confirm the tip's
  `process_spatial_data.py` wins, and read `combine_slices_v3.py` hunk by hunk.
  This is the one piece of Tier A that needs judgement, not mechanics.
- Run `scripts/bil-scripts/find_missing_scripts.sh` on the BIL cluster for
  Tier B; record the outcome in §1.3 either way.
- Census the 87 items: sections, cells, genes, molecules, **largest gene per
  item**, image, species, and whether cells and transcripts share one frame (§2.3).
- On a PSC node verify: the Python floor (§1.2), Apptainer (D4), S3 from
  compute nodes, scratch for a 150–240 GB temp dir.
- Run `prepare.py` unchanged on one BIL h5ad + its transcripts, as a control.

**Gate:** `perf/bil/census.json` + `CENSUS.md` committed; the fleet's largest
gene is a number in this file; D4 answered; nothing in §1.3 still missing.

### Step 1 — the image runs on PSC (`chore/bil-container`)

The worker image (now carrying `prepare*.py` and pinned requirements) built for
`linux/amd64` and converted for Apptainer; one `sbatch` wrapper that runs it.

**Gate:** `prepare.py --help` and a tiny-fixture `prepare` + `validate` run to
completion **inside the container on a PSC compute node**; `pytest tests/python`
passes in the same image.

### Step 2 — sections into one store (`feat/bil-sections`)

**The core of the reprocessing** (§2.3). Combined cells + every section's
molecules in one store: the table keeps `_sample_id` as an obs column; each
section's transcripts get that section's `(dx, dy)` from `sample_offsets.json`
and all go into one `points/molecules` with a `section` column.

**Gate:** a real multi-section BIL item builds one store; `validate_store.py`
passes; cell counts match the per-section sum; **molecules land on their cells**
— for each section, the median per-cell correlation between a cell's expression
and the molecules within 10 µm of it, as `align2.py` checked (r ≈ 0.82 on
12_16_2025, falling to ≈ 0.5 when shifted 50 µm); the viewer draws every
section and colours by `section`.

### Step 3 — the fleet runner (`feat/bil-runner`)

SLURM array over a manifest; per-dataset state with `--resume` / `--force`;
per-stage resources from a table; **every stage exits non-zero on any per-item
failure** (fixes `process_single_molecule.py:1100`); one sync implementation.

**Gate:** a tiny fixture goes end to end as a PSC array job; killing it
mid-stage and re-running resumes without redoing finished work; an induced
failure in one element does not let `sync` run.

### Step 4 — billions (`feat/points-scale`)

Close §1.4 in order: per-worker spill tagging → parallel spill → expose
`band_um` / `tile_size` / `levels` / `--tmp` → resume → a `--sample`
passthrough so validation is affordable.

**Gate:** the census's largest item builds inside its allocation;
`points_zarr.py validate` passes; bytes-per-view at Auto within
[BUDGETS.md](../specs/BUDGETS.md). Numbers land in BUDGETS.md with their run.

### Step 5 — publish the fleet (`chore/bil-rerun`)

The pipeline **emits** the catalog JSON instead of it being hand-maintained
(446 → 87 entries); `scripts/ingest-bil-projects.ts` consumes it unchanged.
Run the fleet. Retire `scripts/bil-scripts/` and the two unused reversers
(`chunked_to_h5ad.py`, `meyes_to_h5ad.py`).

**Gate:** every item COMPLETE with a passing validation, or listed here with
its reason; `ingest-bil-projects.ts --apply` idempotent against the emitted
JSON; no orchestration code outside the runner.

---

## 5. What shows it is right

- `validate_store.py` exits 0 on every published store.
- `spatialdata.read_zarr` (a library we do not maintain) opens every store's
  table, images and shapes.
- The cluster and the worker write **identical stores** from the same input —
  ONE-PIPELINE step 2's parity gate, extended to Apptainer.
- A real multi-section BIL item opens from Explore with every section's
  molecules over its cells.
- BUDGETS.md gains a measured row per real BIL store.
- A second person runs the fleet from this file and `prepare.py --help` alone.

---

## 6. Building context as we go

The repo already has the mechanism; use it rather than inventing one.

- **The decision log is this file.** §3 is append-only with dates; a step that
  changes a decision edits it in the same commit — already [CLAUDE.md §2](../../CLAUDE.md).
- **Specs, not prose.** A format change edits `docs/specs/` or
  `points-zarr-spec.md`; a measurement edits BUDGETS.md. CLAUDE.md §5 ("two
  generations") gets rewritten when the chunked processors retire.
- **Hooks, sparingly:** a `Stop` hook noting that `scripts/` changed without
  this file changing. Everything else belongs in the code as a gate — a check
  in `prepare.py` beats a check in a hook.
- **Subagents** earn their keep on fan-out, not judgement: the Step 0 census
  over 87 items, sweeping the cluster for stranded scripts, mapping an
  unfamiliar script family (how §1 was written). Decisions stay in the main
  session.

---

## 7. How to brief a session on a step

> Read `docs/plans/BIL-CLUSTER.md` (this file), `docs/plans/ONE-PIPELINE.md`,
> `docs/PHILOSOPHY.md` and the specs they name. You are doing **Step N** on
> branch `<branch>`. Its gate is `<gate>`. Measure before you change anything;
> put the numbers in the file the plan names. When a fact or decision here
> turns out wrong, edit this file in the same commit. Do not commit or push
> unless asked.
