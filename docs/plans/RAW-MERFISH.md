# Plan: raw MERFISH runs, checked and converted on the researcher's computer

A researcher with a raw MERFISH run — the lab pipeline's per-field `.npz`
files — asks the agent to sanity-check it, convert it into one finalized
SpatialData store, and view it. Today they convert the `.npz` files into
something MERFISHEYES accepts (h5ad, a molecule table…) themselves, upload
it, and wait for it to process and load. Locally, the MERFISHEYES extension
does all of it in its own Python, once, and the viewer streams the store
from their disk: nothing uploaded, nothing reprocessed.

This is a pillar of the local extension ([MCP.md](MCP.md) pillar D item 7).

Keep this file current: a step that changes a fact or a decision edits it
in the same commit.

> **Status (2026-11-25): blocked on §4.** No step starts until every
> question in §4 is answered by the user and the answers are written into
> this file. The plan below is a draft to be revised by those answers.

## 1. Facts (2026-11-24, from the runs under `/data`)

Nothing in the repository reads this format yet.

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
run folder (.npz per FOV per set)
  │  meyes_raw_check          seconds: headers + a sample of FOVs → summary, plots, flags
  │  (the agent asks §4-style questions of the person)
  │  prepare.py → prepare_merfish_npz.py   a job, resumable, progress per FOV
  ▼
<run>.zarr   points (molecules) · table (cells × genes) · one section per set
  │          [later] image (mosaics) · shapes (segmentation)
  ▼
meyes_open_local → viewer on the person's computer → meyes_view_* / analysis tools
```

## 3. Steps (draft — revised by §4 before any starts)

**0. The format written down (spec + census).**
`docs/specs/raw-merfish-profile.md`: what a run folder must hold, required
and optional columns, units, sets → sections, the overlap rule, cell
identity, blanks. A census script over every run in `/data` lists the
variants (keys, columns, dtypes, folder names). Gate: every run in `/data`
recognised, or refused with a reason.

**1. The sanity check before converting — `meyes_raw_check` (seconds).**
A step in the extension's Python reading headers plus a sample of FOVs:
- summary: sets, FOVs per set, molecules per set, the gene panel, **blank
  fraction** (the misidentification estimate), fraction assigned to cells,
  cells, score / brightness ranges, extent in µm;
- plots, returned as images: (a) each set's spots at `x_um_abs, y_um_abs`
  sampled and coloured by chosen genes (the `XF` + gene plot); (b) an FOV
  grid map coloured by molecules or blank fraction (failed or dim fields
  stand out); (c) one FOV at full detail; (d) per-gene counts against the
  blanks;
- flags: empty or missing FOVs, blank-rate outliers, stitching jumps, sets
  of very different depth.

Gate: Roy's run checked in under a minute, plot (a) included.

**2. Full conversion → one SpatialData store — `prepare_merfish_npz.py`, as a job.**
`prepare.py` recognises the run by content, so `meyes_open_local(<run
folder>)` works. Pass 1 decodes FOV by FOV into parquet parts (bounded
memory, parallel, **resumable** — a laptop that sleeps continues). Pass 2
writes:
- molecules: the points element, every decoded molecule — gene, µm
  x / y / z, set, FOV, unique cell id, score, in_nuc; blanks out of the drawn
  element, kept for QC;
- cells: a table from `cellinfo` (centroid, volume, set, FOV) and counts
  per gene from assigned molecules; blank counts per cell as a QC column;
  DE stats as for every store;
- one section per set (a `section` column);
- then `stats.json` and `validate_store.py`.

The conversion of `meyes_open_local` becomes a job with progress ("FOV
412 / 1,220"): a large run takes minutes to hours and Claude Desktop cuts
calls off after a few minutes.

Gates: Roy's run converted here with time and peak memory recorded in
[BUDGETS.md](../specs/BUDGETS.md); molecule and cell counts equal the source
files' after the overlap rule; pixel checks in the viewer; the same
conversion on a 16 GB laptop.

**3. Images and outlines (after step 2).** Mosaics → OME-Zarr image under
the cells (`image_to_ome_zarr.py`); segmentation masks → cell polygons in
µm (needs the mask → µm placement rule from step 0).

**4. The agent's workflow.** A section of the guide (`src/guide.ts`):
check → ask (sets, score threshold, blanks, where the output goes and the
disk it needs) → convert as a job → open → `meyes_view_*` for QC pictures
(hide a bad FOV) → the analysis tools on the cells. Gate: a researcher in
Claude Desktop goes from "here is my run folder" to a viewed store with no
terminal.

**5. Later.** The same `prepare_merfish_npz.py` in the server worker
([ONE-PIPELINE.md](ONE-PIPELINE.md)), so the identical folder can be
uploaded; the cluster ([BIL-CLUSTER.md](BIL-CLUSTER.md)) can use it.

## 4. Questions — required before any step starts

The user is finalizing how to approach this pipeline. Each answer is
written here (with its date) before work begins.

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

Answers: *(none yet)*
