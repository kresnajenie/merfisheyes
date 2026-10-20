# Plan: one way to prepare and upload data

A tutorial that takes anyone from a platform export — or from parts they
processed themselves — to a folder MERFISHEYES opens as it is, by drag and
drop, by *Upload & Save* or from a public bucket. Judged by the rules in
[docs/PHILOSOPHY.md](../PHILOSOPHY.md): one layout, specified first,
produced by tested scripts, checked by the validator, read back by the
`spatialdata` library.

This file is the brief for every session working on it: what was measured
on 2026-10-19, the decisions, the chapters with their gates, and the prompt
for a session. Keep it current: a chapter that changes a fact or a decision
edits this file in the same commit.

## 1. What is true today (measured on develop `82150a9`)

- Three ways in, producing different cell formats. *Preview in browser* +
  *Upload & Save* takes already-processed parts as they are
  (`lib/local/classify-bundle.ts`: a SpatialData store or a lone AnnData
  zarr, a merfisheyes molecule folder, an OME-Zarr image; parts recognised
  by content). *Upload & process on server* takes raw H5AD, Xenium or
  MERSCOPE plus a transcripts file plus a `.zarr` image
  (`lib/ingest/classify-folder.ts`) and writes the legacy chunked cell
  folder, not a zarr table. The hand route is `aws s3 sync`,
  `scripts/spatialdata_points_to_sm.py`, `scripts/spatialdata_cluster.py`
  and a `mapping.json` (`docs/spatialdata/README.md`).
- Nothing shipped wrote the CSC copy the table profile requires. The
  `spatialdata` writer passes the table to anndata unchanged
  (`_io/io_table.py`); spatialdata-io 0.7.1's Xenium and Visium tables
  are CSR (`sc.read_10x_h5`), its MERSCOPE table dense from the CSV. Only
  the test fixtures added `layers/X_csc`. The validator fails a CSR table
  without a copy.
- Shapes: `spatialdata` writes one `shapes.parquet` per element
  (`_io/io_shapes.py`, geopandas); the profile requires exactly that. A
  directory of parts (dask-geopandas) is a validator failure and is what
  the skin reference store has.
- Molecules: spatialdata-io writes a `points.parquet` directory
  (`_io/io_points.py`, dask) the viewer does not read. Converters:
  `spatialdata_points_to_sm.py` → `bin.gz` folder (kept, decision 4);
  `legacy_to_points.py` → tiled points element (later).
- Visium: no reader anywhere in our pipelines; the only route is
  spatialdata-io `visium()` / `visium_hd()` into a store (steps 2–3 of
  `docs/plans/SDDB.md`, merged).
- Images: spatialdata 0.8 writes NGFF `0.5-dev-spatialdata` in Zarr v3
  (`_io/format.py` `RasterFormatV03`), which the reader accepts; chunk size
  is the caller's (`image_models_kwargs={"chunks": …}`), and the defaults
  gave the 4096² Xenium tiles the SDDB census flagged.
- anndata 0.12.8 writes Zarr **v2** by default
  (`anndata.settings.zarr_write_format == 2`). A v2 table opens from a drop
  (the older local path) but not inside a store or by URL (store-layout
  requirement 1). anndata 0.13 writes v3 only.
- `mapping.json` links were absolute URLs only; `fetchMappingConfig` passed
  them straight to the loaders.
- Environments: spatialdata-io 0.7.1 needs Python ≥ 3.12 and spatialdata
  ≥ 0.7.3; the ingest worker is Python 3.12 / numpy 2.2.6; the format tests
  install `tests/python/requirements.txt`.
- The skin reference store (`spatialdata/skin_adapt_v2/`, bucket
  `merfisheyes-bil`, audit 2026-10-17): molecules as a `bin.gz` folder
  (works); table CSR with root attributes naming a `layers/X_csc` that
  does not exist; shapes a directory of parts (403, not drawn); image
  4096² uint16 chunks of ~20 MB; a `points.parquet` directory nobody reads.
- Existing digests that overlap and disagree in places: five `/docs`
  pages, the ten-tab format guide, `README.md`, `scripts/README.md`,
  `scripts/UPLOAD_GUIDE.md`, `docs/data-requirements.md`,
  `docs/spatialdata/README.md`.

## 2. Decisions (user, 2026-10-19)

1. **The target layout per part**: image OME-Zarr (NGFF 0.5, `c, y, x`,
   1024² chunks); shapes one GeoParquet file per element; molecules the
   `bin.gz` folder (no tiling yet: the points element is still in
   progress); cells an AnnData table in a Zarr v3 SpatialData store with a
   CSC copy in `layers/X_csc` when `X` is not CSC. The viewer finds any CSC
   layer of `X`'s shape; the validator reports a missing one.
2. **Both audiences**: people who click *Upload* on a raw export keep the
   server route as it is; people who preprocess themselves run one script
   per platform and drop, save or sync the folder.
3. **Our own readers stay maintained** (`process_spatial_data.py`,
   `process_single_molecule.py`) beside spatialdata-io.
4. **Clusters and DEGs stay the separate `annotations/` sidecar**
   (`spatialdata_cluster.py`); SpatialData has no place for DEG stats. A
   simple Leiden is run only when the table has no categorical column.
5. **Every combination opens**: cells alone, molecules alone, image alone,
   any two, all three; a store with several elements and no table opens
   as molecules with the image under them. Each combination is named and
   tested, not assumed.
6. **Xenium and MERSCOPE first, Visium and Visium HD after**.
7. **The tutorial lives in `/docs` and in `README.md`**, generated from the
   scripts' help text so the two cannot drift.
8. **The skin reference store is re-prepared** by the tutorial's own
   scripts, so the reference data passes the profiles it is judged by.

## 3. Chapters

Branches off `develop`, one PR each, every PR with unit, Python and
`@quick` gates green, this file updated. Commit only on the user's
go-ahead; never push or open PRs unless asked.

### Chapter 1: parts you already have (`feat/upload-tutorial`)

The folder of three parts — an AnnData table, a `bin.gz` molecule folder,
an OME-Zarr image — made ready by one script, with the combination table.

- `scripts/prepare_folder.py`: finds the parts by content (the rules of
  `classify-bundle.ts`), wraps a lone table into a Zarr v3 store
  (`<name>.zarr/tables/table`; a v3 table is moved, a v2 table or an
  `.h5ad` is rewritten), writes `obsm/spatial` from `obsm/X_spatial`,
  adds `layers/X_csc` when `X` is not CSC and no CSC layer of its shape
  exists, optionally runs Leiden into `annotations/` when the table has no
  categorical column, writes `mapping.json` with names relative to the
  folder (`--url` for absolute), consolidates the root metadata last, and
  runs `validate_store.py`.
- `spatialdata_cluster.py` gains a `cluster()` function the script calls;
  its command line is unchanged.
- `lib/config/dataset-links.ts`: `resolveMappingPaths` resolves relative
  `links` and `image` against the folder the mapping came from, so a
  prepared folder opens from any bucket without editing the file.
  `docs/specs/mapping.md` says so.
- Tests: `tests/python/test_prepare_folder.py` (each table form, the CSC
  copy's counts equal to `X`, a CSC table untouched, a rerun a no-op, two
  candidates refused, no coordinates refused, molecules + image without
  cells, clustering only without a categorical column);
  `tests/unit/dataset-links.test.ts`.
- Docs: `/docs/prepare` (new page: the layout, the script, the
  combinations, the bucket) and a README section pointing at it;
  `docs/spatialdata/README.md` "Preparing a store" points at the script.

Gate: the Python tests green locally and in the `format-specs` job; the
unit suite green; the `@quick` local-bundle specs unchanged and green; a
folder prepared by the script dropped on the homepage opens with the
three parts overlaid (verified in the browser).

Result (2026-10-19): built as above. Python 32 passed (12 new), unit 300
passed, tsc clean, the `@quick` local-bundle and SpatialData specs 29
passed. A folder made of a Zarr v2 AnnData table, the tiny molecule
folder and the fixture image, prepared with the real command line, was
dropped on the homepage by a throwaway Playwright spec: 40 cells, the
image under them, three genes loaded whole. Found on the way: anndata
refuses to write into a group whose `zarr.json` carries consolidated
metadata (a `write_zarr` v3 table has it), so the script strips it from
the root and the table before editing and consolidates last; the Next
dev server refuses a symlinked `node_modules`, so a worktree needs its
own `npm install`.

Cold walkthrough (2026-10-19, a fresh `python -m venv`, following only
`/docs/prepare` as a bioinformatician with an h5ad, a transcripts parquet
and a TIFF). Gaps found and closed in the same chapter:

- No way to get the scripts: the page assumed a checkout. Now a first step
  (clone, venv, `scripts/requirements-prepare.txt`, one file for every
  script on the page). The fresh install resolved anndata 0.13.4 and zarr
  3.4.1, so the tests run against those too (anndata 0.13 writes v3 only).
- No command for the molecules or the image: the page named the scripts
  and left the TIFF to "any OME-Zarr writer". Now `process_single_molecule.py`
  with its column flags, and a new `scripts/image_to_ome_zarr.py` (a plane,
  a channel stack or RGB, uint8/uint16, 1024² chunks, pyramid, `scale` and
  `translation` from `--um-per-pixel` and `--origin`), checked by the
  validator's image rules through `check_image_folder`.
- Nothing said the three parts must share units and frame, and nothing
  checked it. Now the script compares each part's extent with the cells'
  (`molecules_box` from the first gene files, `image_box` from level 0):
  no overlap or a span over 5× or under 0.2× is a `WARN` naming both
  extents; the page explains the line and shows a sample report.
- Memory: a v2 or h5ad table is read whole; the page says so and how to
  avoid it (write v3, the table is moved).
- The docs-drift test now pins all four commands on the page and in the
  README.

Left open for a later chapter: a Seurat object (convert to h5ad first),
labels, a TIFF larger than memory (tile-wise writing), Windows paths.

### Chapter 2: Xenium (`feat/upload-tutorial-xenium`)

- `scripts/prepare_xenium.py`: `spatialdata_io.xenium(path,
  cells_as_circles=False, image_models_kwargs={"chunks": (1, 1024, 1024)})`
  → `write()` → `spatialdata_points_to_sm.py --drop-z` for the transcripts
  → `prepare_folder.py` on the result (CSC, mapping, validator). Labels
  are written but not read (profile backlog). Shapes come out as one
  file each from `spatialdata`.
- A pinned `scripts/requirements-prepare.txt` (spatialdata-io 0.7.1,
  spatialdata 0.8.0, Python ≥ 3.12), separate from the ingest worker's.
- Test on the `scverse/spatialdata-sandbox` Xenium sample in the Python
  job (downloaded, cached); `spatialdata.read_zarr` reads the result.
- `/docs/prepare` gains the Xenium chapter; README likewise.

Gate: validator green on the sandbox sample and on one real Xenium export
on this machine; the folder opens with image, polygons, molecules and
cells overlaid.

Result (2026-10-20): built as `scripts/prepare_xenium.py`, with two
departures from the sketch above. No labels by default (`--labels` to
keep them): the viewer does not draw them and they are the largest part
of the store; the reader then still ties the table to `cell_labels`, and
the viewer places the image relative to the region it names, so a region
that is not in the store would have put the image in pixels under cells
in microns (4.7× off). The script re-points the table at
`cell_boundaries`, which carries the same cell ids. And no sandbox
download: the test (`tests/python/test_prepare_xenium.py`) runs on an
export named by `MERFISHEYES_XENIUM` and is skipped otherwise; CI does
not download gigabytes. `spatialdata_points_to_sm.py` gained a
`convert()` function. The image gets `scale_factors` explicitly: without
them spatialdata-io writes a single level, and the viewer loads the
coarsest level whole as the backdrop.

Measured on `Xenium_V1_FF_Mouse_Brain_Coronal_outs` (Xenium 1.0.2,
130,870 cells, 248 genes, 7 M transcripts, 827 MB morphology focus):
69 s wall, 8.9 GB RSS at the peak, a 2.9 GB store (image 6 levels of
1024² chunks, one 13.7 MB GeoParquet per boundary element) and 293 MB of
molecule files; `--cluster` 47 s and 1.8 GB for 33 Leiden clusters.
Validator: every MUST OK; the expected WARNs are the unread
`points/transcripts`, the polygon files' single row group and
`nucleus_boundaries` not being a name the viewer draws. The gated test
(validator, store contents, `spatialdata.read_zarr` reading the table and
the polygons back, 93 s) passes; dropped on the homepage by a throwaway
spec the folder opens with 130,870 cells, the image, the polygon layer
and three genes loaded.

Cold walkthrough (fresh venv, `pip install spatialdata-io` on top of the
requirements file): the page's commands ran as printed. Gaps closed: the
reader's `region` column counted as a cell-type column, so `--cluster`
silently did nothing on every spatialdata-io store (`has_categorical_column`
now skips the join keys and single-category columns, with a test); the
page states the measured time, memory and sizes and which WARN lines are
expected. Left open: polygon files as one row group (the shapes profile's
later ask: row groups by tile), `nucleus_boundaries` never drawn, an
export without `cells.zarr.zip` (pre-1.3 CSV bundles) untested, memory
on a laptop with 8 GB (the reader holds the table and the transcripts
while writing).

### Chapter 3: MERSCOPE (`feat/upload-tutorial-merscope`)

- `scripts/prepare_merscope.py` on `spatialdata_io.merscope()` (the
  reader tested on 138 BIL directories on 2026-10-16: index fixes may be
  needed until spatialdata-io PR 428 lands), the same tail as Xenium; the
  dense table becomes CSC in `layers/X_csc`.
- The existing `process_spatial_data.py` MERSCOPE path is documented as
  the alternative for the legacy chunked folder (decision 3).

Gate: as chapter 2, on a sandbox MERFISH sample and one real MERSCOPE
export.

### Chapter 4: Visium and Visium HD (`feat/upload-tutorial-visium`)

- `scripts/prepare_visium.py` on `visium()` / `visium_hd()` with
  `fullres_image_file`; spots and bins draw as the cells' footprint (SDDB
  step 3); the default table gets its CSC copy; the hires and lowres
  images stay (the viewer picks the finest).

Gate: validator green, the H&E under the spots on a sandbox sample.

### Chapter 5: the skin reference store, re-prepared

- Run chapter 2's script on the skin Xenium export (or rewrite the shapes
  as one file, add the CSC copy, rechunk the image, with
  `prepare_folder.py`), sync to the bucket, record the validator's
  verdict in `docs/specs/AUDIT-2026-10-17.md`'s successor.

Gate: `validate_store.py` passes on the skin store through its URL.

### Later, out of this plan

- The tiled points element as the molecule format (when it is finished).
- Raw cells beside processed parts previewed together (backlog in
  `docs/spatialdata/README.md`).
- Retiring any of the overlapping guides.

## 4. What shows it is right

- `scripts/validate_store.py` passing on every folder the scripts prepare,
  in CI and on the skin store.
- The `spatialdata` library reading each prepared store back.
- The `@quick` upload and local-bundle specs green on prepared fixtures,
  through the real homepage cards.
- One person outside the team producing a working upload from the pages
  alone, on their own export, with no chat.
- Before that person: every chapter ends with a **cold walkthrough** — a
  fresh virtual environment, the page followed line by line as a
  bioinformatician who has only the platform's export (or an h5ad, a
  transcripts table and a TIFF), every command run as printed, every gap
  (a missing step, an unnamed dependency, an unstated assumption about
  units, memory or file layout) written into the chapter's result and
  closed in the same PR.

## 5. How to brief a session on a chapter

```
Read docs/plans/UPLOAD-TUTORIAL.md, docs/PHILOSOPHY.md and docs/specs/
store-layout.md, table-profile.md, shapes-profile.md, image-profile.md and
mapping.md. Do chapter <N> of the plan: <its title>. Start from the raw
material the plan names (the spatialdata-io readers at 0.7.1, the
spatialdata writers at 0.8, the sandbox samples, the exports on this
machine), not from the existing guides. Every command in the tutorial is a
tested script that ends in the validator. Then do the cold walkthrough of
§4: a fresh venv, the page followed as a bioinformatician with only the
export, every gap closed and recorded in the chapter's result. It is done
when the chapter's gate holds, the Python, unit and @quick suites are
green, and this plan, /docs and README are updated. Commit only when
told; never push or open a PR unless told. List every place you had to
touch outside the files the chapter names.
```
