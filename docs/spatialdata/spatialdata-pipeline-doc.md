# The SpatialData pipeline

How BIL datasets become one SpatialData store per dataset, run by you on the
PSC cluster or on any machine with the pinned environment.

This is the reference for the whole process: what it needs, what it does,
what it writes, and how to check it worked. For *why* it is being built this
way and what is still open, see
[docs/plans/BIL-CLUSTER.md](../plans/BIL-CLUSTER.md).

---

## 1. What changed, in one page

**The old format.** One gzipped file per gene, for the whole dataset:
`genes/Slc17a7.bin.gz` held every molecule of that gene everywhere. To draw
one gene the viewer downloaded all of it — at BIL scale `Slc17a7` is ~589 M
molecules, roughly 7 GB, before anything appeared. Cells were the same idea:
one `expr/chunk_*.bin.gz` per gene.

**The new format.** One **SpatialData store** (a Zarr v3 folder) per dataset,
holding the cells and the molecules together in one coordinate system:

```
<sample>.zarr/
  tables/table        the cells: counts (CSC), obsm/spatial, obs columns
  points/molecules    the molecules: a tiled pyramid, 8 levels
  shapes/…            polygons, when the input has them
  images/…            a microscopy image, when the input has one
```

The molecules are stored **twice at every level, on a regular grid**:

| | sorted by | answers |
|---|---|---|
| `tile_major` | tile | "everything in this region" |
| `feature_major` | gene | "this gene, everywhere" |
| `index/` | — | turns either question into byte ranges |

Level 0 is 100 µm tiles holding every molecule; each level up is a uniform
quarter of the one below, to 8 levels. On `ace-dip-tug`:

| level | tile | grid | molecules |
|---|---|---|---|
| 0 | 100 µm | 193 × 178 | 433,023,021 |
| 3 | 800 µm | 25 × 23 | 6,764,744 |
| 7 | 12,800 µm | 2 × 2 | 26,494 |

So zoomed out the viewer reads 26 thousand points, not 433 million; zoomed in
it reads only the tiles on screen. **The unit of reading became the view
instead of the gene.** That is the whole point of the change.

> **"Store" means the `.zarr` folder** — SpatialData's word for it, used
> throughout the scripts and this doc.

**The other half of the change:** cells from many slices were already
combined into one frame; molecules were not, so the two could never be drawn
together. Now the combine records each slice's exact shift and a later stage
moves that slice's molecules by the same amount.

---

## 2. Components

Everything lives in `scripts/`; the cluster wrappers in `scripts/bil-scripts/`.

### The work

| Script | Does | Notes |
|---|---|---|
| `combine_slices_v3.py` | many slices → one cell frame | grid-shifts each slice so they do not overlap; tags every cell `_sample_id`; writes **`sample_offsets.json`** (each slice's shift *and* where its molecules are) and `artifact_mask_p*.csv` |
| `map_my_cell.py` | cells → cell types | `--method correlation` (fast, needs only the reference `.h5`) or `hierarchical` (needs `cell_type_mapper`) |
| `prepare.py` | a table → a store | the one entry point, also used by the web app's server. Recognises h5ad, Xenium, MERSCOPE, Visium, a molecule table, or a folder of parts, **by content** |
| `slices_to_points.py` | every slice's molecules → `points/molecules` | reads `sample_offsets.json`; shifts each slice onto the cells' frame; refreshes `stats.json` |
| `validate_store.py` | the gate | applies the same rules the viewer applies; non-zero exit means do not publish |

### The cluster wrappers

| File | |
|---|---|
| `env.sh` | **the only place paths, the bucket and the interpreter live.** `meyes_run` dispatches to a venv or `apptainer exec` |
| `combine_slices.sbatch` | stages 1 and 2 (`--mask-only` for the second) |
| `map_my_cell.sbatch` | stage 3 |
| `prepare_store.sbatch` | stage 4 |
| `slices_to_points.sbatch` | stage 5 |
| `validate_store.sbatch` | stage 6 |
| `s3_sync_store.sbatch` | stage 7 |
| **`launch_spatialdata_pipeline.sh`** | **chains all seven** with `--dependency=afterok` |

### Testing

| | |
|---|---|
| `tests/bil/` | 14 tests on generated fixtures, ~30 s, no data needed |
| `tests/bil/run_dry_run.sh` | the same seven stages without SLURM, for a real dataset on one machine |

---

## 3. What it needs

### Input

Three families, and they do **not** all take the same route:

| Family | Looks like | Combine? | Route |
|---|---|---|---|
| **MERSCOPE multi-slice** | `<slice>/…/cellpose_metadata.csv`, `cellpose-cell-by-gene.csv`, `cellpose-detected_transcripts.csv` | **yes** | all seven stages |
| **h5ad** | `cell_by_gene.h5ad` + `segmented_spot_table.csv` | no — one slice per sample | `prepare.py` directly |
| **Xenium** | the folder holding `experiment.xenium` | no | `prepare.py` directly |

File names are matched **fuzzily**, so `cellpose-detected_transcripts.csv`
and `detected_transcripts.csv` are both found. Slices may sit several
directories down; the combine does a breadth-first search for them and skips
anything this pipeline wrote itself (`combined_output`, `mmc_output*`,
`sm_output`, `*.zarr`).

### Environment

Python **3.12** with `scripts/requirements.txt` pinned — the same list the
web app's worker installs, so a store made here is the store the server would
make. This is **not** the old `~/merfisheyes_env`, which the legacy
`process_*.py` still use. Keep both.

```bash
python3.12 -m venv ~/merfisheyes_xarr_env
source ~/merfisheyes_xarr_env/bin/activate
pip install -r ~/merfisheyes/scripts/requirements.txt
```

`spatialdata` is required — `prepare.py` cannot write a store without it.

MapMyCells needs the **14 GB reference** (`/bil/data/meyes/mapmycells-reference`),
read in place, never copied into an image.

### Resources and time

Measured on `ace-dip-tug` — 3 slices, 40.6 GB of CSV, 407,859 cells,
433,023,021 molecules — on a 64-core machine:

| Stage | Request | Measured |
|---|---|---|
| combine | 2 CPU, 64 G | 1 m 41 s |
| mask | 2 CPU, 64 G | 13 s |
| annotate | 64 CPU, 512 G | 1 m 07 s, 5.6 G peak (correlation) |
| prepare | 8 CPU, 128 G | 17 s |
| points | 8 CPU, 128 G | 8 m 28 s |
| validate | 2 CPU, 32 G | 1 m 15 s |
| **total** | | **11 m 58 s, 4.9 GB peak RSS, 6.0 GB store** |

**Memory does not grow with the molecule count** — nothing is ever
concatenated; each slice streams in batches. **Time** is what grows.

**Disk:** the points stage spills to `<store_dir>/.molecules.tmp` before
building the levels, ~16 bytes per molecule (7 GB for 433 M, ~180 GB at the
11.24 B scale). That is on the **store's** filesystem, not `TMPDIR`.

---

## 4. How to run it

### On the cluster

```bash
cd ~/merfisheyes/scripts/bil-scripts

./launch_spatialdata_pipeline.sh samples.csv --dry-run   # print the chain, submit nothing
./launch_spatialdata_pipeline.sh samples.csv             # run it, no upload
./launch_spatialdata_pipeline.sh samples.csv --sync      # run it and upload
./launch_spatialdata_pipeline.sh ace-dip-tug /bil/data/18/aa/.../input   # one sample
```

**Always `--dry-run` first.** It prints every `sbatch` line and the dependency
graph without submitting, which is the only way to check the chain before it
costs anything.

`samples.csv` is `sample_name,input_path`, the same file the other launchers
take. A **third field sets that dataset's mask percentile**:

```csv
ace-dip-tug,/bil/data/18/aa/.../input,75
ace-low-cop,/bil/data/.../h5ad_sample,0
ace-dip-vat,/bil/data/.../input
```

| Option | |
|---|---|
| `--sync` | upload to `s3://bil-psc-xarr/<sample>/` (off by default) |
| `--species human` | MapMyCells taxonomy (default `mouse`) |
| `--method correlation` | MapMyCells method (default `hierarchical`) |
| `--mask-percentile N` | default for every row; `0` attaches no mask |
| `--skip-combine` | `combined_output` is already there |
| `--skip-mmc` | reuse `mmc_output/mapping_output.csv` |
| `--dry-run` | print, submit nothing |

Each stage waits on `afterok` of the one before, so **a failure stops that
sample's chain and leaves the others running**. Nothing half-written is
published, because sync waits on validate.

### On one machine (no SLURM)

```bash
MEYES_PYTHON=~/merfisheyes_xarr_env/bin/python \
  tests/bil/run_dry_run.sh INPUT OUT NAME \
    --reference ~/data/merfisheyes-test-data/mapmycells-reference \
    --method correlation --mask-percentile 25
```

`--labels mapping_output.csv` reuses an MMC result instead of computing one.

### Settings

All in `env.sh`, overridable from the environment:

| | default |
|---|---|
| `MEYES_BASE` | `/bil/data/meyes` — output root |
| `MEYES_S3_BUCKET` | `bil-psc-xarr` |
| `MEYES_REFERENCE` | `/bil/data/meyes/mapmycells-reference` |
| `MEYES_PYTHON` | `~/merfisheyes_xarr_env/bin/python` |
| `MEYES_IMAGE` | empty; set it and every stage runs under Apptainer instead |

---

## 5. What it outputs

```
/bil/data/meyes/<sample>/
  combined_output/
    cell_metadata.csv        every slice's cells, shifted into one frame
    cell_by_gene.csv         counts, same row order
    sample_offsets.json      each slice's (dx, dy) + where its molecules are
    artifact_mask_p*.csv     p5 … p80, one row per cell
    check_*.png              sanity plots
  mmc_output/
    mapping_output.csv       cell_id = <sample_id>_<cell_id>, one row per cell
    check_celltype_*.png
  store/
    <sample>.zarr            THE OUTPUT
    mapping.json             names the store
    stats.json               cells, genes, molecules, versions
```

**Only `store/` is synced.** The bucket holds one folder per dataset.

> The bucket needs **public read and CORS**, or the viewer cannot open a store
> by URL: a public URL cannot be listed, so the store's consolidated metadata
> is the only way the viewer learns what is in it. See
> [docs/specs/store-layout.md](../specs/store-layout.md).

---

## 6. How to check it worked

**The gate, first:**

```bash
python scripts/validate_store.py <store>.zarr --mapping <store_dir>/mapping.json
```

Every line is `OK`, `WARN` or `FAIL`. Any `FAIL` means do not publish.
`WARN tables/table: no instance_key` is expected for MERSCOPE combined output
(there are no polygons to match).

**Then the numbers:**

```bash
cat <store_dir>/stats.json     # cells, genes, molecules — molecules must not be 0
```

**Then by hand, if you want to be sure:**

```python
import sys; sys.path.insert(0, "scripts")
import points_zarr as pz, anndata as ad

store = "…/ace_dip_vat.zarr"
g, attrs = pz.open_element(f"{store}/points/molecules")
print(attrs["points:n_points"])                 # total molecules
df = pz.read_dataframe(f"{store}/points/molecules", level=0,
                       order="tile_major", columns=["x", "y", "sample"])
print(df["sample"].value_counts())               # molecules per slice

obs = ad.read_zarr(f"{store}/tables/table").obs
print(obs["_sample_id"].value_counts())          # cells per slice
```

The three checks that matter:

1. **Molecules per slice ≈ the source row counts**, less blanks and control
   probes (which are dropped by one shared rule).
2. **Each slice's molecules sit on that slice's cells** — compare the x/y
   range of `df[df["sample"] == s]` with `obs[obs._sample_id == s]`'s
   `center_x`/`center_y`. They should overlap almost exactly.
3. **The cell-type labels are on the right cells** — the MMC `cell_id` is
   `<sample_id>_<cell_id>`, so join on that, never on row order.

---

## 7. The artifact mask

MERSCOPE CSV datasets carry a lot of cells with little or no expression, and
cells outside the section boundary. The combine computes a mask at every
percentile from p5 to p80 (`artifact_mask_p*.csv`, one row per cell with
`is_artifact`) and plots them in `check_artifact_thresholds.png`.

**How much to filter is a judgement about that tissue** — some datasets need
p75, h5ad datasets usually need none — so it is not automated. Today:

- the percentile is **per dataset**, the third field of `samples.csv`;
- the mask is attached as an **`is_artifact` obs column**, and **no cell is
  dropped**. The filtering is the viewer's to do, and it is reversible.

This differs from the legacy pipeline, which deleted the masked cells.

---

## 8. Known gaps

| | |
|---|---|
| **Apptainer image** | not built. `env.sh` is already shaped for it (`MEYES_IMAGE`) |
| **The sync stage has never run** | written, preflight checks pass, never talked to S3 |
| **`cell_type_mapper`** | needed for MapMyCells `hierarchical`; `correlation` does not need it |
| **Parallel spill** | the points writer takes a per-worker `tag` but the driver is single-process: ~8 m 30 s per 433 M molecules, so ~3.6 h at 11.24 B |
| **Old `combined_output` folders** | predate `sample_offsets.json`, so their molecules cannot be placed — the combine must be re-run (it is cheap) |
| **h5ad `cell_polygons.geojson`** | not read, so those stores have no `shapes/` |
| **Xenium `cell_id`** | deliberately not used to split molecules, to match the browser preview |

---

## 9. If something breaks

| Message | Means |
|---|---|
| `no sample_offsets.json` | the combine did not run, or ran before this feature existed. Re-run it. |
| `no cell_metadata.csv in …` | the combine did not produce output; check stage 1's log |
| `exit 2` from `prepare.py` | the input was not recognised — check the file names against §3 |
| `exit 3` from `prepare.py` | recognised but incomplete (e.g. molecules with no table) |
| `exit 4` from `prepare.py` | written but the validator failed it |
| `ModuleNotFoundError: cell_type_mapper` | use `--method correlation`, or install it |
| cells doubled | the combine read its own output — only possible on a build older than this one |
