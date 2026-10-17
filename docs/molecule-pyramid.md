# Molecule pyramid (proposal) — a multiscale, spatially indexed molecule store

Status: **prototype working, branch `feat/points-zarr`.** Written 2026-10-16
after the storage research for very large molecule datasets. The format is
specified in [points-zarr-spec.md](points-zarr-spec.md); the proof of concept
is the 12_16_2025 mouse brain dataset, all ten sets aligned into one
SpatialData store ([below](#proof-of-concept-the-12_16_2025-dataset)).

## Goal

One way of storing molecules that:

- **stays interactive in the viewer at any size**: a gene can be drawn across
  the whole tissue or zoomed into a few cells, and each view downloads an
  amount bounded by the screen, not by the dataset;
- **works for analysis too**: per gene, per region, per cell, with room for
  more per-molecule fields later;
- **builds on Zarr conventions** others are converging on, so SpatialData,
  spatiato, Celldega and we can read the same thing.

## Today

`process_single_molecule.py` writes one file per gene,
`genes/{gene}.bin.gz`: a gzipped flat float32 `[x, y, z, …]`, in bit-reversal
order so the first bytes are an even sample of the whole tissue (see
[SM-STREAMING-PLAN.md](SM-STREAMING-PLAN.md)). Genes with ≥ 2 M molecules are
streamed.

What that can't do:

- **No spatial index.** Zoomed into 100 µm, the viewer still downloads the
  gene's molecules from the whole tissue.
- **Nothing per molecule but x, y, z.** No cell id, no QC.
- **All genes in a small area** means downloading every gene file in full.
- A gene of 1 B molecules is a 12 GB file.

## What we measured

### Skin (`skin_adapt_v2.zarr`, SpatialData / Celldega)

| | |
|---|---|
| Molecules | 74,011,892 |
| Gene names | 10,017 |
| Extent | about 9.1 mm across (x, y in µm) |
| Molecules per gene | median 5,439 · p90 33,477 · p99 136,847 · max 1,766,363 (EEF1G) |
| Source | `points/transcripts/points.parquet`, 36 parts, 1.6 GB, 13 columns, one row group per 53 µm tile (171 × 82 tiles) |

### Each order alone fails one kind of query

**Space-first** (rows grouped by tile, as the skin parquet is): one gene
across the tissue decodes nearly every row group. Spatiato measured this on a
Xenium 5K section (136.6 M molecules, 5,122 genes): the median gene is in
2,968 of 7,294 tiles, so drawing it decodes ~85 % of the section to show
~10 k molecules.

**Gene-first** (our first idea: a per-gene pyramid where each molecule sits in
exactly one level, 4,000 per box). Simulated on the skin molecules
(`scratchpad/pyramid/sim.py`):

| Query, all genes | Molecules inside | Molecules the gene-first pyramid reads |
|---|---|---|
| 50 µm square | 20,651 | 30.1 M (41 %) |
| 200 µm square | 354,377 | 38.0 M (51 %) |
| 1 mm square | 6.8 M | 46.6 M (63 %) |

It was excellent for one gene (the biggest, 1.77 M molecules, in 932 boxes;
a 19 k-molecule gene in 7), and grouping genes into shards brought the file
count from 44,731 to 274 per field. But it can't answer "all genes here", and
capping each box at 4,000 points flattens a gene's hotspots until you zoom in.

**Conclusion (shared with spatiato's proposal for SpatialData): store both
orders.**

| View | Rows it needs | Order that keeps them together |
|---|---|---|
| All genes, small region | Every row in a few tiles | Space-first (tile-major) |
| Few genes, whole tissue | One gene's rows in thousands of tiles | Gene-first (gene-major) |
| Few genes, region | One gene's rows in some tiles | Gene-first, as short runs |
| All genes, whole tissue | Too many to draw | A coarser sampled level, either order |

## Design

### Conventions it builds on

| Convention | Status | Used for |
|---|---|---|
| [multiscales](https://github.com/zarr-conventions/multiscales) (uuid `d35379db-88df-4056-af3a-620245f8e347`) | v0.1, pilot; breaking changes expected before v1 (planned before end of 2026) | Declaring the levels |
| [dataframe](https://github.com/ilan-gold/zarr-convention-dataframe) `df:` (uuid `1f1b36d2-3185-6da8-91f6-9a0823af2d33`) | maturity 0, gathering feedback | Each copy of the molecules is a dataframe of 1-D columns |
| [enum](https://github.com/ilan-gold/zarr-convention-enum) `enum:` (uuid `c906c423-56ed-413c-9943-7b2ff52d18f2`) | maturity 0 | Gene and cell id columns: integer codes → names |
| Ours, `points:` ([spec](points-zarr-spec.md)) | — | Column names, tile grid, index arrays |

No Arrow: x, y, z are fixed-size floats, so a decoded chunk is a
`Float32Array` as it is. The only copy is interleaving x/y/z for the GPU.

### Layout

The format is specified in **[points-zarr-spec.md](points-zarr-spec.md)**
(normative). In short:

```
points/transcripts/                  zarr.json: ngff:points + multiscales + points: attributes
├── 0/                               every molecule
│   ├── tile_major/                  df: x, y, z, gene (enum), set (enum) — sorted by (tile, gene, r)
│   ├── feature_major/               df: x, y, z, gene, id — sorted by (gene, tile, r)
│   └── index/                       tile_rows, feature_rows, feature_tile_ptr, feature_tile, feature_tile_count
├── 1/                               a uniform quarter of 0, same layout
└── …                                until a level holds about a million molecules
```

- `r` is one random 32-bit number per molecule; level *k* is the molecules
  with `r < 2^(32 − 2k)`, so levels are nested uniform samples.
- `id` is a molecule's row in level 0 `tile_major`, at every level.
- Tiles: square, 100 µm at level 0, doubling per level; numbered row-major.
- Point columns: 65,536-row chunks, 64 chunks per shard; coordinates
  `blosc` (zstd + shuffle), ids delta-coded.

## How the viewer reads it

| Action | Reads |
|---|---|
| Open | Root metadata (consolidated), plus each level's `feature_rows` and `gene_tile_ptr` — about 80 KB per level for 10 k genes |
| Select a gene | Its tile list and counts at each level (small). Pick the **finest level whose count in view fits the point budget** (shared across selected genes), then read `feature_major` x, y, z for the tiles in view as a few row ranges. Whole tissue = one range |
| Zoom / pan | The same, for newly visible tiles; loaded ranges stay cached |
| All genes, zoomed in (later) | `tile_major` rows of the visible tiles at a fitting level |
| Hover a molecule | `id` → level 0 `tile_major` row → cell id and other fields |

A gene with few molecules is always drawn from level 0, completely; only
abundant genes fall back to a sample, and every sample keeps relative density.

## How analysis reads it

- **The full table** is level 0 `tile_major`. The `df:` layout mirrors
  AnnData's dataframe layout, so a pandas reader is a few lines.
- **A region**: `tile_rows` → the rows of the tiles overlapping it, then an
  exact filter.
- **A gene**: `feature_rows` → one row range of `feature_major`.
- **Molecules per cell / the cell × gene matrix**: one scan of `gene` and
  `cell_id`.
- **A new per-molecule field** (another segmentation's cell id, QC flags): one
  more `column_n` in level 0 `tile_major` and a longer `columns`. Nothing else
  is rewritten. Copying it into `feature_major` is optional and listed in its
  metadata.

## Writer

1. Read gene, x, y, z, cell id from the source (SpatialData points parquet,
   Xenium `transcripts.parquet`, MERSCOPE `detected_transcripts.csv`). Gene and
   cell id become enum codes.
2. Draw `r` per molecule with a fixed seed (recorded in metadata), so a rerun
   on the same input gives the same rows.
3. For each level *k* (molecules with `r < 4^-k`): compute tiles; sort by
   (tile, gene, r) → `tile_major` (level 0 defines `id`); sort by
   (gene, tile, r) → `feature_major`; write the index arrays from the counts.
4. Check: level 0 holds every molecule once in each order, ids match, index
   counts sum to the row counts, every level is contained in the one below.

The skin data (74 M) fits in memory. Datasets of billions need the sorts done
out of core (partition by tile or gene first) — part of the large-dataset
test.

## Decisions and alternatives

| Choice | Alternative considered | Why this one |
|---|---|---|
| Both orders in one store | Gene-first only, parquet beside it for regions | One store, shared ids, analysis without a second format |
| Nested levels (each a sample of the one below) | Disjoint levels, each molecule once (Neuroglancer annotations) | `multiscales` readers expect level 0 to be complete; uniform samples keep density honest. Cost: levels duplicate (~1.33× total) and zooming in re-downloads instead of adding |
| Row numbers in index tables | Box = chunk, capped at 4,000 points | Index independent of chunk layout; matches the `df:` convention's 1-D columns |
| float32 µm | uint16 relative to the tile (half the size) | Readable by anything without decoding; revisit as a codec once measured |
| Point budget per view | Cap per box / per level | The counts give it before reading; one setting for all genes |
| Zarr | Parquet | Columns appended without rewriting; read granularity independent of the index; one Zarr store with the image and cells |

## Costs, risks, open questions

- **Storage**: two orders × about 1.33 for the levels. Compressed sizes are
  the first thing to measure (skin parquet: 1.6 GB for 13 columns).
- **Draft conventions**: names and layout may change.
- **`multiscales` for points**: the transform is the identity and the tile
  size is our own property. Worth raising with the convention's owners.
- **Requests**: zarrita fetches each chunk separately; an abundant gene's
  whole-tissue read could be dozens of requests until neighbouring ranges are
  merged in `lib/storage/fetchRange.ts`.
- **Where it lives**: beside the SpatialData store, linked from
  `mapping.json`, for the prototype. Spatiato puts its cache at
  `points/<name>/transcripts_vis_zarr` inside the store; to align later.
- **Not yet decided**: 3D (atlas space) pyramids; more fields than x, y, z,
  gene, cell id; one store for several sections.

## Prototype

| Part | Where |
|---|---|
| Writer, reader, validator (Python) | `scripts/points_zarr.py` |
| Tiny fixture | `scripts/testdata/make-tiny-points.py` → `tests/data/tiny/zarr/points_v1.zarr` |
| Viewer reader | `lib/spatialdata/points.ts` (element, index, row ranges; a byte cache so a chunk is fetched once) |
| Viewer layer | `lib/spatialdata/PointsLayer.ts` (per gene: finest level that fits a 4 M point budget shared by the selected genes; one read per tile row in view; a view stays up until the next has arrived) |
| Wiring | `lib/spatialdata/pointsMolecules.ts`; `components/spatialdata-viewer.tsx` and `lib/local/open-bundle.ts` open a store's points; `components/three-scene.tsx` draws them with the gene picker and legends |
| Tests | `tests/unit/points-zarr.test.ts` |

Success: any gene, any zoom, a bounded download per view; every molecule
recoverable from level 0.

## Test datasets

1. **12_16_2025 mouse brain, sets 1–10 — the proof of concept.** See below.
2. Skin — `s3://merfisheyes-bil/spatialdata/skin_adapt_v2/` (74 M molecules,
   10,017 genes). Not run yet.

## Proof of concept: the 12_16_2025 dataset

**Goal:** one SpatialData store holding the combined cells of all sets and
the molecules of every set, aligned to the cells, with the molecules in this
new format under `points/` — opened in merfisheyes with seamless molecule
viewing at any zoom. Processed locally first and previewed in the browser.

**Decisions (2026-10-16):** all ten sets; the cells rebuilt from the
merfisheyes cell folder (the original AnnData is gone and was processed
further); molecules carry x, y, z, gene and set only (custom MERFISH+, no cell
ids); the molecules live in the store under `points/`.

### The data

| Part | Where | What |
|---|---|---|
| Cells, all sets combined | `s3://merfisheyes-bil/12_16_2025_all_sets_meyes/` | merfisheyes cell folder: 1,587,309 cells, 500 genes, 2D, `batch` `set1`…`set10`, `class_name`, `subclass_name`, `supertype_name`, `cluster_name`, `leiden`; expression log-normalized |
| Molecules, one folder per set | `s3://merfisheyes-bil/12_16_2025_sm_mereyes_set{1…10}/` | merfisheyes molecule folders: 500 genes, 3D, `genes/{gene}.bin.gz` |

Both folders keep coordinates **normalized to [-1, 1]** with a
`scaling_factor` but **not the centre**: multiplying back gives µm about an
unknown centre. The cells' combined frame is used as is; each set's
molecules get a translation into it (step 1), which absorbs the unknown
centre.

Two things in the molecule folders: the manifests count about 2 % more
molecules than the files hold (the files are taken as the truth), and two
gene files are named with `_` for `-` (`Nkx2_1`, `Nkx2_4`).

| Set | Cells | Molecules (manifest) |
|---|---:|---:|
| set1 | 21,153 | 170,430,511 |
| set2 | 60,400 | 492,861,658 |
| set3 | 109,356 | 762,810,910 |
| set4 | 129,745 | 844,837,332 |
| set5 | 272,340 | 1,796,173,910 |
| set6 | 231,618 | 1,676,279,217 |
| set7 | 194,703 | 1,454,802,111 |
| set8 | 259,852 | 1,978,845,274 |
| set9 | 172,412 | 1,268,049,582 |
| set10 | 135,730 | 1,018,264,205 |
| **all** | **1,587,309** | **11,463,354,710** (the files hold 11,239,601,084) |

### Step 1 — align each set's molecules to the cells (done)

Translation only (x, y; z kept), per set, by `align2.py`:

1. Molecules × the set's `scaling_factor`; cells in µm (× 22,173.96).
2. For the 60 genes with most molecules: the set's cells weighted by their
   expression of the gene, against the gene's molecule density. Summed
   cross-correlations, coarse at 40 µm over the whole set, then 4 µm within
   ±200 µm, then a 1 µm search for the most molecules within 5 µm of a cell.
3. Check: each molecule to its nearest cell (≤ 10 µm); per cell, the
   correlation of those counts with the cell's expression.

| Set | Translation (µm) | Fine-peak sharpness | Molecules ≤ 10 µm of a cell | Median per-cell r | r if shifted 50 µm |
|---|---|---:|---:|---:|---:|
| set1 | (−12,139.0, 2,728.7) | 6.3 | 0.73 | 0.84 | 0.51 |
| set2 | (−3,495.2, 2,420.8) | 12.1 | 0.76 | 0.84 | 0.52 |
| set3 | (6,279.8, 1,734.6) | 3.2 | 0.81 | 0.83 | 0.52 |
| set4 | (17,323.1, 1,101.6) | 3.1 | 0.81 | 0.83 | 0.53 |
| set5 | (18,025.6, −11,452.9) | 6.3 | 0.81 | 0.81 | 0.50 |
| set6 | (9,261.5, −10,016.3) | 6.0 | 0.81 | 0.81 | 0.49 |
| set7 | (2,301.7, −9,092.1) | 5.6 | 0.79 | 0.82 | 0.49 |
| set8 | (−4,597.2, −18,801.0) | 4.6 | 0.78 | 0.82 | 0.48 |
| set9 | (−5,020.9, −7,293.1) | 3.5 | 0.78 | 0.82 | 0.51 |
| set10 | (−11,600.8, −7,218.4) | 3.4 | 0.77 | 0.82 | 0.51 |

Sharpness is the best match over the next best outside ±100 µm (1 = no
clear match). Every set lands on its cells: the molecules reproduce each
cell's expression (r ≈ 0.82) and lose it 50 µm away (r ≈ 0.5).

### Step 2 — one SpatialData store (`build_store.py`)

```
12_16_2025_all_sets.zarr/
├── zarr.json            spatialdata_attrs, consolidated metadata
├── tables/table/        the 1.59 M cells (written by spatialdata 0.8)
└── points/transcripts/  11.24 B molecules in the points-zarr format,
                         columns x, y, z, gene, set; 8 levels
```

Out of core: 32 workers spill each (set, gene) file into 800 µm bands on
disk; then per level, band by band, sort and append (tile-major), spill per
gene, and append per gene (feature-major).

**Full build (2026-10-16):** 11,239,601,084 molecules, 8 levels (11.24 B →
2.81 B → 702 M → 176 M → 43.9 M → 11.0 M → 2.74 M → 686 k), tiles 100 µm to
12.8 mm on a 43.2 × 33.9 mm grid. About 2 h on this machine (32 workers to
spill, one process to sort and write; peak memory about 42 GB):

| Step | Time |
|---|---|
| Spill all files into bands | 7 min |
| Level 0, tile-major (35 bands) | 72 min |
| Level 0, feature-major (500 genes) | 20 min |
| Levels 1–7 | 25 min |

The run stopped on level 6's last band, which sent no molecule on to level 7
(fixed in the writer); levels 6–7 were rebuilt from level 5 with
`PointsWriter.resume`, a fresh uniform quarter as the spec allows.
`points_zarr.py validate` passes on every level (200,000 sampled rows per
level matched against level 0; levels 4–7 checked to sit inside the level
below).

| Part | Size |
|---|---:|
| Table (1.59 M cells) | 1.1 GB |
| Level 0 `tile_major` | 82 GB |
| Level 0 `feature_major` | 84 GB (x 25, y 24, z 36, gene 0.02, id 0.27) |
| Level 0 index | 31 MB |
| Level 1 / 2 / rest | 59 / 17 / 6 GB |
| **Store** | **248 GB, 38,030 files** |

Level 0 holds about 14.8 bytes per molecule for both orders; the source
molecule folders are 64 GB (gzipped x, y, z, one order, no levels). z is the
largest column — its values span only about 30 µm, so quantizing it is the
first saving to try.

Test build on set 1 alone (167 M molecules): 130 s, 3.7 GB of points.

### Step 3 — view it

- By URL: serve the folder (with byte ranges) and open
  `/spatialdata-viewer/from-s3?url=<store URL>`.
- From disk: drop the store folder on the homepage (browser preview).

Measured on the set-1 store, three abundant genes (17.6 M, 6.7 M, 6.6 M
molecules), headless Chromium:

| View | Levels (per gene) | Molecules drawn | Read for this view |
|---|---|---:|---:|
| Whole section (5 mm) | 2, 2, 2 | 1.9 M | 18 MB |
| 1.5 mm | 2, 1, 1 | 2.7 M | 20 MB |
| 400 µm | 0, 0, 0 (every molecule) | 2.0 M | 15–20 MB |
| 100 µm | 0, 0, 0 | 0.3 M | nothing new |
| Back to the whole section | 2, 2, 2 | 1.9 M | nothing new |

Opening the folder from disk: 17 s from drop to molecules on screen (most of
it the 1.59 M cells).

Measured on the **full store** (11.24 B molecules), the three most abundant
genes — Slc17a7 589 M, Slc1a3 286 M, Gfap 284 M molecules:

| View | Levels (per gene) | Molecules drawn | Read for this view |
|---|---|---:|---:|
| All ten sections (40 mm) | 5, 4, 4 | 2.8 M | 29 MB |
| One section (9 mm) | 3, 2, 2 | 3.6 M | 58 MB |
| 1.5 mm | 0, 0, 0 (every molecule) | 4.8 M | 55 MB |
| 400 µm, 200 µm, 100 µm | 0, 0, 0 | 0.4 M, 0.14 M, 0.09 M | nothing new |
| Back to all ten sections | 5, 4, 4 | 2.8 M | nothing new |

Zoomed out, molecules are kept at least 1 px on screen; without that, the
viewer's sub-pixel culling hides them (the same happens with molecule
overlays today, which never had a whole 40 mm dataset on screen).

### Compatibility

spatialdata 0.8 reads the table (`sd.read_zarr(path, selection=("tables",))`)
and lists `transcripts` as a points element it can't load; reading everything
fails on the unknown points version (`0.3-zarr`), as the spec says.

## References

- [zarr-conventions/multiscales](https://github.com/zarr-conventions/multiscales)
- [ilan-gold/zarr-convention-dataframe](https://github.com/ilan-gold/zarr-convention-dataframe)
- [ilan-gold/zarr-convention-enum](https://github.com/ilan-gold/zarr-convention-enum)
- [spatiato multiscale points cache format](https://github.com/vibspatial/spatiato/blob/main/src/spatiato/core/multi_scale_cache_points_zarr/CACHE_FORMAT.md)
- [Neuroglancer precomputed annotations and sharding](https://github.com/google/neuroglancer/tree/master/src/datasource/precomputed)
- Report: [Molecule Pyramids](https://claude.ai/artifact/TRzerL6E9tBjEZTNLhGxYQ) (earlier research, gene-first design)
- [SM-STREAMING-PLAN.md](SM-STREAMING-PLAN.md) — today's progressive per-gene files
