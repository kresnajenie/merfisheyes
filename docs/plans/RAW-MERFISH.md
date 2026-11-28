# Plan: raw MERFISH runs, checked and converted on the researcher's computer

A researcher with a raw MERFISH run — the lab pipeline's `.npz` output —
asks the agent to sanity-check it, convert it into one finalized
SpatialData store, and view it. Today they convert the `.npz` files into
something MERFISHEYES accepts (h5ad, a molecule table…) themselves, upload
it, and wait for it to process and load. Locally, the MERFISHEYES extension
does all of it in its own Python, once, and the viewer streams the store
from their disk: nothing uploaded, nothing reprocessed.

This is a pillar of the local extension ([MCP.md](MCP.md) pillar D item 7).

Keep this file current: a step that changes a fact or a decision edits it
in the same commit.

> **Status (2026-11-28): steps 0–2 done** (the check, the conversion; the
> converted store viewed and confirmed by the user). Next: step 3. §4 is
> answered: Jess's processed runs are the baseline — the shape a run must
> have before conversion.

## 1. Facts

### 1a. The baseline: Jess's runs (2026-11-27)

Jess's runs under `/media/guacamole/Jess/NeuronTracing/` (read only):
`20260326_WT_BRBB_Mer` (18 sets, 115 GB of spots), `MouseA56` (20 sets,
with GFP), `MouseA54/20260522_Sample3` (with GFP). Her notebooks
(`BRBB_MERFISH_analysis/NeuronTracing/Step4_RNAMosaic`, `Step5_cellbyGene`,
`Step6_Scanpy`) are how they were made.

**`final_spots/spots_allth3_final_set<K>.npz`** — one file per set, every
member stored uncompressed (so read in place, in pieces):

| member | content |
|---|---|
| `XF` | float32 N × 5: `z` (plane, × 0.4 µm), `x`, `y` (stitched pixels, × 0.108333 µm), `scoreA`, field index into `fovs` |
| `genes` | `<U13`, one name per molecule (4 GB for 78 M molecules — larger than `XF`); 500 genes + 125 `blankNNNN` |
| `pix_size` | `[0.4, 0.108333, 0.108333]` µm for z, x, y |
| `fov_size` | `[50, 3000, 3000]` planes / pixels |
| `pos_fovs_pix`, `fovs` | each field's stitched position (pixels) and name (`Conv_zscan9_00`) |

- **De-duplicated already**: a molecule is kept only in the field whose
  centre is nearest (Chebyshev distance, `KDTree(...).query(p=inf)`,
  Step4). Checked in the data: density in overlaps 12.2 vs 13.8 per bin
  elsewhere; 97 % of overlap molecules sit in their nearest field.
- **Score-filtered already**: `scoreA > -0.5`. Blanks are still there.
- A set is one tissue section: the 18 sets of `20260326_WT_BRBB_Mer` are
  serial sections imaged at separate places on one slide, ~2 mm apart, none
  overlapping (`POS/pos_K.txt`; confirmed with Jess, 2026-11-28). Nothing is
  stitched between sets.
- One stitched frame per set, origin at the set's lowest field
  (`Pmin`); µm = `XF[:, :3] × pix_size`.
- Rows are ordered by field; within a field they are not spatially ordered,
  so evenly spaced blocks are a fair sample. Reads: 370 MB/s contiguous, ~20
  ms per scattered read on the share; 64 blocks of 4,000 rows ≈ 3.5 s a set.
- Molecules carry **no cell id**.

**`final_mats/<field>--mat.npz`** — per field: `Mcells` (`Conv_zscan9_00--1`,
unique across the run), `Mgenes` (the 625), `Mcounts` uint32 cells × 625,
`Xum` (z, x, y µm, the set's frame). **Not de-duplicated**: 402,576 cells
in the 1,220 files of `20260326_WT_BRBB_Mer` vs 313,076 after Step6's rule
(set 9: 18,552 → 14,126; every dropped cell is in an overlap). Loading all
1,220 takes 2.4 s. Molecules were assigned to cells by distance to the
segmentation's pixels (≤ 2.5 µm, Step5).

**`final_matsWithGFP/`** (GFP runs) — the same files plus `ctsGFP`, GFP spots
per cell. **`GFP_xfit/XGFPF_set<K>.npz`** — GFP spots, float32 N × 8
(columns not yet known).

**`Segmentation2d/<field>--<run>_set<K>--segm.npz`** — `segm` uint16
(50, 1000, 1000) per field, downsampled 3×.

**Her finished h5ad** (`combined_adata_control.h5ad`) — 313,076 × 500,
`obs_names` equal to `Mcells`, MapMyCells labels, leiden, UMAP, blanks in
`obsm['X_blank']`, `X_spatial_z` per-set µm, `X_spatial` sets side by side.

### 1b. Other layouts (2026-11-24, the runs under `/data`) — not the first target

**Spot files** — `saveDecodedFinalV1--MERFISH_Analysis[_groupN]--Conv_zscan__NNN--_setK.npz`,
one per field of view (FOV) per set, in `FinalDecodedSpots/`,
`final_spots*/` or `final_spots_corrected/`:

| key | content |
|---|---|
| `XF` | float32 (float64 in Snonia's runs), N × 17–18, one row per molecule |
| `header` | `XF`'s column names: `z, x, y` (in the field), `abs_brightness, cor, brightness, color, mean_bightness_variation, mean_distance_variation, index_from_XH_pruned, score, cell_id, ifov, iset, [in_nuc], z_um_abs, x_um_abs, y_um_abs` |
| `genes` | one name per molecule; about half the codewords are blanks (`blank0062`) |
| `cellinfo` | one row per cell |
| `header_cells` | `cell_id, volm, zc, xc, yc, z_um_abs, x_um_abs, y_um_abs` |

**Variants** — the 2025-04-01 3×3 run uses `XFall` / `genesall` and has
no `cellinfo`; `in_nuc` is in some runs only; dtypes differ. Detection must
go by the keys and columns a file holds, not its name.

**Beside the spots** — `Segmentation/*dapi_segm.npz` (label masks, int32
(30, 1000, 1000) for a (30, 3000, 3000) field: downsampled 3×),
`mosaics/*.tiff` (stitched images; 51 GB in Roy's run), `DataFrames/*.pkl`
(the lab's downstream tables; pickles, not opened: loading one runs code).

**Scale**

| run | size | FOV files | molecules |
|---|---|---|---|
| Roy, `RoyEdu_C_ASO_4week__1mopInj_11_21_2023` | 29 GB | 1,220 | 661 M (sets 206 M / 230 M / 226 M) |
| Stephen, `ALSspinalcord_DCSSBB_11_30_23/final_spots_corrected` | 100 GB | 10,210 | not counted |
| Snonia, `HttHuman_Snonia__08_03_2023/final_spots4` | 1.7 TB (folder) | 75,642 (folder) | not counted |

**Geometry and identity**
- Fields are ~326 µm wide on a ~312 µm step: neighbours overlap by
  ~14 µm, so molecules and cells in the overlaps appear in two files.
- `cell_id` repeats across fields (the same range in every FOV; 0 =
  unassigned). A cell is (set, FOV, cell_id).
- 62–78 % of molecules are assigned to a cell (Roy, first FOVs of set 1).
- Molecule counts come from the `.npy` headers inside each `.npz` without
  loading the arrays: all 1,220 of Roy's files in ~35 s on a cold disk.

**What exists to build on** — `transcripts_to_points.py` / `points_zarr.py`
write the streamed molecule element out of core (70 M rows ≈ 1 min, 1 GB;
Roy's run ≈ 10 min, to be measured); `prepare.py` dispatches by content;
`de_stats.py`, `validate_store.py`, `image_to_ome_zarr.py`; the extension's
analysis jobs (`scripts/meyes_analysis.py`, `src/local/analysis.ts`).

## 2. Shape

```
run folder (final_spots/ + final_mats[WithGFP]/, the baseline in §1a)
  │  meyes_raw_check          a minute: headers + blocks + the mats → summary, plots, flags
  │  (the agent asks the person: sets, which h5ad's labels, where the output goes)
  │  prepare.py → prepare_merfish_run.py   a job, resumable, progress per set
  ▼
<run>.zarr   points (molecules) · table (cells × genes) · one section per set
  │          [later] GFP spots · image (mosaics) · shapes (segmentation)
  ▼
meyes_open_local → viewer on the person's computer → meyes_view_* / analysis tools
```

## 3. Steps

One module, `scripts/raw_merfish.py`, reads the baseline (detect, members
in place, decode gene names, the overlap rule); the check and the
conversion both use it.

**0. The format written down.** *Done 2026-11-27.* `docs/specs/raw-merfish-profile.md`: what a
run folder must hold (§1a), required and optional members, units, sets →
sections, the overlap rule, cell identity, blanks, GFP. Detection goes by
members and columns, not names; a folder that is not the baseline is
refused with what it lacks (the §1b layouts among them). Gate: Jess's three
runs recognised; a §1b run refused with a reason. **Met**: `raw_merfish.py
detect` reads `20260326_WT_BRBB_Mer` (18 sets, 1.71 B molecules, 1,220
mats), `MouseA56` (20 sets, 908 M, 1,158 GFP mats) and `MouseA54/…Sample3`
(21 sets, 502 M); Roy's run is refused ("per-field saveDecodedFinalV1
files: that layout is not read yet").

**1. The sanity check — `meyes_raw_check` (about a minute).** *Built
2026-11-27* (`scripts/raw_merfish.py check`, the tool in
`src/tools/analysis.ts`, a job like the analysis steps;
`tests/python/test_raw_merfish.py`). Python in the extension's environment
reads headers, 8 blocks of 50,000 molecules per set and the mats:
- summary per set: fields, molecules, cells before / after the overlap rule,
  **blank fraction** (misidentification estimate), score range, extent in µm;
  for the run: the panel, genes with too few molecules;
- plots, returned as images: (a) every set's cells at their µm position,
  coloured by one gene's molecules per cell (from the mats: seconds); the
  `XF` + `genes` plot — every molecule of each chosen gene in one set, one
  full pass through that set's `genes` and `XF`; (b) the field grid per set,
  coloured by molecules per cell or blank fraction (a failed or dim field
  stands out); (c) per-gene counts against the blanks;
- flags: missing fields (in `fovs` but no mat, or the reverse), fields with
  few molecules or a high blank rate, spots that look not de-duplicated,
  sets of very different depth.

Gate: `20260326_WT_BRBB_Mer` checked in under a minute, plot (a) included.
**Met** (2026-11-27, read from the lab share): 11 s for all 18 sets (cells,
fields, codewords), 31 s with every Slc17a7 and Gad2 molecule of set 1
(35 M molecules scanned); 2.5 GB peak. MouseA56 10 s, MouseA54 11 s.

What it taught (kept here, the reasons for the code):
- Molecules can't be sampled spatially in seconds on the share: random
  reads cost ~15 ms each and are served one at a time (np.memmap page
  faults worse still), while contiguous reads run at ~370 MB/s; within a
  field, rows run along x, so a block is a stripe. Hence cells (mats) for
  the all-set picture, contiguous blocks for the statistics, a full pass
  for the molecule plot. Members are read with `os.pread`.
- Gene names are mapped by a hash of their characters, not `np.unique`
  (a set's `genes` is 4–7 GB of `<U13`).
- The overlap rule on the mats reproduces her h5ad to within a tie per run
  (the profile's "Fields imaged twice").
- Thresholds that held on three runs: a field is flagged below half its
  set's median molecules per cell or above 3× its blank fraction; a set
  outside ½–2× the run's median molecules per cell or per field; genes
  below the median blank; blanks above 10× it. MouseA54's flags were real:
  set 13 at 0.3× the run's molecules per cell, set 21 at 9.2× (another
  library?), 5 fields with no mat, 29 with no GFP mat.

**2. Full conversion → one SpatialData store — `prepare_merfish_run.py`, a job.**
*Built 2026-11-27* (`scripts/prepare_merfish_run.py`, `prepare.py`'s
`merfish_run` kind with `--sets`, `--labels-from`, `--resume`;
`meyes_open_local(path, sets, labels_from, output)` converts as a job
(`trackJob`) that resumes from `<output>.partial`;
`tests/python/test_prepare_merfish_run.py`). Decisions taken with the user
(2026-11-27): the output goes where the person says (default the cache);
a few sets can be converted first; molecules get no cell id until the
outlines exist (step 3). Blanks follow the shared gene filter
(`should_filter_gene`): out of the drawn molecules and of `X`, counted per
cell as `n_blank`. `score` is not written on the molecules (4 bytes × every
molecule for a column the viewer does not use; the run is already
filtered) — a change from §4.2's first wording.
`prepare.py` recognises the run by content, so `meyes_open_local(<run>)`
works. Per set (resumable — a laptop that sleeps continues at the next set):
- molecules: the points element — gene, µm x / y (z with `--keep-z`),
  `section` (copied to the feature-major order), `field`; blanks out;
- cells: the table from `final_matsWithGFP` when present, else
  `final_mats`, with Step6's overlap rule applied to `Xum`; counts × the 500
  genes; blank counts per cell (`n_blank`), `gfp_count` from `ctsGFP`;
  DE stats as for every store. An h5ad of the same run, when the person
  names one, adds its labels and embeddings by cell name;
- sections: one per set, on a grid read left to right, top to bottom (as
  `combine_slices_v3.py` does: rows along −y; 250 µm apart), y flipped so
  each section reads as in the lab's plots (§4.6): x = x_set + dx,
  y = dy − y_set, the shifts in the table's `uns['merfisheyes_sections']`.
  One spatial frame only (`obsm['spatial']`); from an h5ad, only
  embeddings (`X_umap`, …) are taken, not its layouts (`X_shifted`);
- then `stats.json` and `validate_store.py`.

Gates: `20260326_WT_BRBB_Mer` converted here, time and peak memory in
[BUDGETS.md](../specs/BUDGETS.md); molecule counts equal the files',
cell counts equal her h5ad's (313,076); molecules land on their cells;
the same conversion on a 16 GB laptop.

Measured (2026-11-27, lab share → local NVMe, 4 workers): sets 1 + 9 of
`20260326_WT_BRBB_Mer` — 21,022 cells (all matched to her h5ad's labels),
111.8 M molecules (1.03 M blanks left out) — in 252 s: 44 s reading, 167 s
building the first level (single-threaded, `points_zarr.finish`), 18 s DE
stats; 1.7 GB store; 5.2 GB peak; `validate_store.py` passes. Molecules
land on their cells: per gene, a cell's count correlates with the molecules
within 5 µm of its centre at r = 0.70–0.78 (0.05–0.26 with the cells
moved 40 µm).

The whole run (2026-11-27, same setup): 313,075 cells (313,075 matched to
her h5ad's labels; her 313,076th is the tie in the profile), 1,695,800,334
molecules (13.5 M blanks left out) — equal to the spot files' count less
the blanks — in 47 min (9 min reading, 30 min the first level); a 22 GB
store (13 bytes a molecule); validates. **Peak 26 GB**: pass B of the
points writer sorts one 800 µm band at a time, and 18 sections side by side
make a band 21 mm wide — 120 M molecules at ~140 bytes each. Fixed by sizing
the bands to the data (`points_zarr.band_for`: ~25 M points a band, a whole
number of tiles): sets 1 + 9 then peak at 3.7 GB (was 5.2). The whole run
again: **6.0 GB peak, 23 min** (7 min reading, 9 min the first level —
thinner bands sort faster too; the share was partly cached the second
time), the same counts and store, validates. Under the 16 GB a laptop has.
The store, for checking by hand: `/data/merfisheyes_claude_extension/raw
conversion tests/brbb_full` (22 GB).

Checked by hand (2026-11-28, the user's Claude session, an audit of that
store against the raw files): cells' positions and counts equal the mats'
exactly; molecules on their cells (log-correlation of gene fractions
0.997); the fields' orientation right as written. Found and fixed: her
h5ad's `X_shifted` (another layout of the sets) was carried as an
embedding — now only embeddings are taken; the grid's rows went up and each
section showed upside down (the viewer draws +y upwards) — rows now go
down and y is flipped (§4.6); 80 fields imaged but never decoded in 7 sets
(`POS/pos_K.txt` lists more than the spot files hold: a decode stopped at
field 90) — upstream, now flagged by the check and logged by the
conversion. Left as is: the repeated fields at each set's end (her
de-duplication split those tiles' molecules between the two copies, so
dropping one leaves holes); the viewer's default molecule genes (Igfbpl1,
Drd1, Th, hard-coded for every dataset in `lib/utils/auto-select-genes.ts`)
— the viewer's, not this plan's. Reconverted with the fixes (2026-11-27):
32 min, 6.0 GB peak, same counts, validates; per section the molecules'
extent equals the cells' within 2 µm; set 1 top-left, sets 16–18 the
bottom row; the copy in `raw conversion tests/brbb_full` replaced.

Still open: the 16 GB laptop itself; the first-level build is
single-threaded (`points_zarr.finish`, shared by every molecule
converter) and is most of the time.

**3. GFP spots, images and outlines.** GFP spots (`GFP_xfit`) as their
own entry in the points once the columns are known; mosaics → OME-Zarr image
under the cells; `Segmentation2d` masks → cell polygons in µm.

**4. The agent's workflow.** A section of the guide (`src/guide.ts`):
check → ask (sets, an h5ad's labels, where the output goes and the disk it
needs) → convert as a job → open → `meyes_view_*` for QC pictures (hide a
bad field) → the analysis tools on the cells. Gate: a researcher in Claude
Desktop goes from "here is my run folder" to a viewed store with no
terminal.

**5. Later.** §1b layouts turned into the baseline by an adapter; the same
`prepare_merfish_run.py` in the server worker
([ONE-PIPELINE.md](ONE-PIPELINE.md)) and on the cluster
([BIL-CLUSTER.md](BIL-CLUSTER.md)).

## 4. Questions — answered 2026-11-27

These were required before any step started.

1. **Overlaps.** Does the lab pipeline already drop duplicates in FOV
   overlaps, or do we? If we do: keep each molecule and cell in the field
   whose centre is nearest?
2. **"Final" spots.** Are `saveDecodedFinalV1` / `final_spots*` already
   filtered by score, or should a threshold be offered? Which folder is
   the one to read when a run has several (`final_spots2` vs
   `final_spots4`, `FinalDecodedSpotsFineSegmCorrect`, `…_corrected`)?
3. **Sets.** Is each set a separate section, and are the `*_um_abs`
   coordinates one stitched frame per set?
4. **The 3×3 run** (`XFall` / `genesall`, no cell table): molecules only,
   or is its cell table elsewhere?
5. **Who produces this format.** Only this lab, or collaborators running
   the same pipeline? Decides how forgiving detection must be.

**Answers (user, 2026-11-27):** Jess's runs are the baseline — "some
people process differently but hers is a good baseline for what the data
needs to look like before conversion."

1. **Overlaps** — done upstream for molecules (`final_spots`, nearest field
   centre). The per-field mats are not; the converter applies the same rule
   to cells, as her Step6 does. The check flags spots that look
   un-de-duplicated.
2. **"Final" spots** — `final_spots/`, already `scoreA > -0.5`. No threshold
   offered. (`score` is left off the converted molecules: step 2.)
3. **Sets** — one section each; each set is its own stitched frame
   (pixels × `pix_size`); sections go side by side in one store.
4. **The 3×3 run and the other §1b layouts** — not targeted first: they are
   refused with a reason until an adapter brings them to the baseline.
5. **Who produces it** — the lab, processed in different ways; the baseline
   is what a run must look like. Detection is by members and columns, and a
   refusal says what is missing.
6. **Orientation** — first "fixed downstream, kept as written"; changed by
   the user on 2026-11-28 after viewing the converted store (the sections
   showed upside down, the grid bottom to top): y is flipped in the store
   and the grid's rows go down, so sections read as in the lab's plots.
7. **Cell table** (Claude's suggestion, taken 2026-11-27) — built from the
   mats, not from a finished h5ad: the check and conversion come before
   analysis, and not every run has one. An h5ad of the same run (`obs_names`
   = `Mcells`) can add its labels and embeddings by name.
8. **GFP** (Claude's suggestion, taken 2026-11-27) — `final_matsWithGFP`
   when present, its `ctsGFP` as a `gfp_count` cell column; the GFP spots
   (`GFP_xfit`) wait for their columns to be confirmed (step 3).
