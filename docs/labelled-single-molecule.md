# Labelled single-molecule datasets (`/lm-viewer`)

A third dataset type alongside single-cell and single-molecule. Points are
molecules, each carrying **N categorical labels and no expression matrix**.

Built for the spiralia embryo data, where a molecule's identity is three labels
— gene, RNA domain, cell — rather than one gene plus a value. Neither existing
viewer fits: the single-cell viewer expects an expression matrix, and the
single-molecule viewer knows only about genes.

First dataset: `MER6-2_E3_1` — 3,140,608 molecules, 213 genes, 11 domains,
26 cells.

---

## The one idea

All three identifiers are categorical, so all three are **obs columns**, and the
dataset has no `expr/` directory at all. That is why the existing
`ChunkedDataAdapter` could load it with a single change: making
`expr/index.json` optional, gated on `manifest.has_expression`.

Filtering is **OR within a menu, AND across menus**:

```
shown = gene ∈ {…}  AND  domain ∈ {…}  AND  cell ∈ {…}
```

The OR is not a design choice — a molecule has exactly one value per column, so
checking two genes can only mean "either". An empty menu imposes no constraint.

---

## Producing a dataset

Two steps. The first is dataset-specific, the second is generic.

### 1. Source → parquet (spiralia-specific)

`scripts/spiralia/export_parquet.py` reads the dill `*_RNA_domain_annotated`
object and writes a flat parquet. It resolves several things that are not
obvious from the source data — see `scripts/spiralia/QUESTIONS.md` for how each
was established:

- `Xh` is in **pixels**; the voxel size is **0.5 µm in z, 0.177 µm in x/y**,
  derived by solving against the µm columns in `domain_meta.csv` (exact to
  1e-10 across every domain).
- The per-molecule gene index is **`icodesN`**, not `node_ids` — the latter is a
  306-entry lookup table.
- The 306-entry panel is **213 real genes + 80 blank controls + 13 duplicate
  probe entries**. The duplicates carry 657,877 raw molecules but **zero** after
  `molecules_keep`, so that flag alone is sufficient.
- `analysis_mask` is exactly `molecules_keep` minus the blanks, which also means
  the `unassigned` RNA domain *is* the blank control probes, not unassigned
  tissue.

```bash
python scripts/spiralia/export_parquet.py MER6-2_E3_1
```

### 2. Parquet → chunked folder (generic)

`scripts/process_labelled_molecules.py` is generic in its columns and has no
spiralia knowledge. It imports `encode_obs_binary`, `build_obs_palette` and the
gene-filter rules from the existing scripts rather than reimplementing them, so
the binary layout stays single-sourced with the browser's `decodeObsColumn`.

```bash
python scripts/process_labelled_molecules.py molecules.parquet output/ \
    --category gene=feature_name \
    --category domain=rna_domain_anno \
    --category cell=cell_id --label cell=cell_name \
    --drop-unassigned cell_id \
    --drop-control-genes feature_name \
    --name "MER6-2_E3_1"
```

`--label` handles columns keyed on an id but displayed by name: cell names
repeat (two polar bodies are both `pb`), so the menu keys on `cell_id` and shows
`pb (1)` / `pb (3)`.

### Output layout

```
manifest.json          coord_format: "q16", has_expression: false
coords/spatial.q16.bin.gz
obs/metadata.json
obs/{gene,domain,cell}.bin.gz    dictionary-encoded, uint8 codes
palettes/{gene,domain,cell}.json
```

`obs/*.bin.gz` is **byte-identical to the single-cell obs format** (`bin-v1`):
header, JSON dictionary, then one uint8/16/32 code per molecule, gzipped.

**Coordinates are quantised to uint16 per axis.** Float coordinates barely
compress — 37.7 MB gzipped to only 32.1 MB, since the mantissa is high-entropy —
so this took the dataset from 34.7 MB to **19.5 MB**. Max error 0.003 µm against
a source already rounded to 0.01 µm. This is LM-only; the shared f32 path is
chosen whenever the manifest does not declare `q16`.

### Publishing

```bash
aws s3 sync output/ s3://<bucket>/<prefix>/<name>/
```

Two things that silently break it:

- **Never set `--content-encoding gzip`.** The app gunzips in JS via
  `DecompressionStream`; if S3 advertises the encoding, the browser decompresses
  first and every binary read throws.
- The bucket needs **CORS** allowing GET from the app origin.

When *replacing* a dataset in place, set `Cache-Control: no-cache` on
`manifest.json` and `obs/metadata.json`. The adapter now revalidates both, but
that only protects caches populated after the header exists.

---

## Viewing

```
/lm-viewer/from-s3?url=https://<bucket>.s3.<region>.amazonaws.com/<prefix>/<name>
```

`expr/index.json` returning 403/404 on every load is expected and absorbed.

| | |
|---|---|
| rail | Gene / Domain / Cell menus, molecule-size slider |
| hotkeys | `G` `D` `C` switch the colouring column, `H` hides the UI |
| legends | active selections; eye to hide (⌘-click to solo), name to recolour, X to remove |
| scene | hover for all three labels; ⌘-click eases the camera to a molecule; double-click toggles that molecule's value for the colouring column |
| camera panel | selected/unselected size, reset view, owner "save current view as default" |
| share | full viewer state, including camera pose, encoded into `v=` |

Genes are coloured on selection from a 10-slot palette so the few you pick stay
distinct; domain and cell use fixed palettes generated at ingest.

---

## Rendering

The material is **opaque with depth writes** (`transparent: false,
depthWrite: true`). The first version used `transparent: true,
depthWrite: false`, which produced two symptoms at once: nothing occluded
anything, and every one of 3.1M points shaded even when buried. That is the same
tradeoff `createSmPointCloud` documents in `lib/webgl/point-cloud.ts`, and why
the single-molecule viewer never had the problem.

Consequence: **partial alpha cannot be expressed.** Unselected molecules are a
solid grey backdrop rather than a faint haze, and there is no opacity control.

Selection lives in three **lookup textures**, one texel per category, so a
checkbox toggle re-uploads under a kilobyte. The per-point index attributes are
static and uploaded once. `dotSize` must be back-calculated from the data extent
(`targetPx * distance / (baseSize * proj11)`) — getting it wrong by orders of
magnitude puts every point under the sub-pixel cull and silently blanks the
scene with no error.

---

## Verification

Ground truth computed in Python from the source parquet, before any viewer code
existed:

> gene ∈ {NANO2, ACTC} ∧ domain ∈ {3q, 3Q} ∧ cell ∈ {3B, 3D} = **32,665**
> of 3,140,608

The live UI reports the same number, and it survived a column rename and the
coordinate format change. Every obs column and the coordinate buffer round-trip
row-exact against the parquet. `tests/unit/labelled-molecule-lut.test.ts` covers
the LUT and intersection logic.

Re-check that number after any change to the filter, the LUT or the ingest.

---

## What still needs doing

### 1. Segmentation masks and 3D cell meshes — **blocked**

The largest remaining piece, and the reason the viewer currently shows molecules
with no cell geometry.

The data already exists. Under `Segmentation/Bogdan/` each embryo has:

| artefact | shape (E3_1) | voxel size |
|---|---|---|
| `Old/*_segm.npz['segm']` | 370 × 725 × 725 | 1.0 / 0.708 / 0.708 µm |
| `AutoRefined_Masks/*_cell_mask_DS4.npz` | 93 × 182 × 182 | 4.0 / 2.832 / 2.832 µm |
| `ManualCurate_Masks/*_curated_cell_mask_DS4.npz` | 93 × 182 × 182 | same |

and **prebuilt triangle meshes** in HDF5, 26 cells for E3_1, in two variants:

- `annotated_domain_viz_cell_meshes/` — marching cubes on the curated mask.
- `annotated_domain_viz_pointcloud_nonoverlap_cell_meshes/` — built from the
  molecule cloud, guaranteed non-overlapping, and **the only one carrying a
  `cell_identity` attribute** per mesh.

Each mesh group holds `faces`, `vertices_xyz` and `vertices_reoriented_xyz`.
Vertices are in **µm**, and `vertices_xyz` sits inside the molecule bounding box
— which is why the ingest keeps coordinates raw rather than normalising.

**The blocker (question A3):** nobody has confirmed what
`vertices_reoriented_xyz` is. If it is a per-embryo alignment, the meshes and the
molecules will not co-register without that transform, and it is **not stored in
the h5**. Until that is answered, building the mesh layer risks being wrong in a
way that looks plausible.

Also unresolved: which of the three masks is canonical, and which mesh variant
to render.

### 2. Upload and registration

The viewer is reachable only via `?url=`. There is no DB row and no
`/lm-viewer/[id]`, so datasets are not discoverable or linkable by id. Ownership
(claim banner, owner-saved camera) works against a by-URL lookup and needs the
dataset claimed once.

### 3. Per-value colour and size

`colorOverrides` and `sizeOverrides` are per menu **and per value** in the store,
and the shader already reads both out of the palette LUT. Colour has a picker in
the legend; **per-value size has no UI**. Worth having: gene abundance spans
21 → 256,402 molecules, so a rare gene is easily swamped by a common one.

### 4. The other 44 embryos

The ingest is generic and the other embryos are ready, but **question E1** —
whether the 306-probe panel is identical across MER1 / MER2 / MER5-1 / MER6 — is
unanswered. Note the mesh coverage split: **31 of 45 embryos have meshes, and the
14 without them are exactly the 14 MER2 embryos.**

### 5. Smaller items

- Explore listing and the plot panel, both deferred.
- Opacity was removed with the blended renderer; restoring it needs either
  blending back or dithered transparency.
- If the grey backdrop ever buries the selection at some angles, the fix is a
  two-pass render: selection opaque, backdrop depth-tested but not
  depth-writing. The usual ordering artefact does not apply because the backdrop
  is one flat colour.

### 6. Open questions for the data author

`scripts/spiralia/QUESTIONS.md` — 20 questions across segmentation, coordinates,
filtering, naming and scope. **Still unsent.** A3 and E1 above are the two that
block work.

---

## Files

| | |
|---|---|
| `scripts/process_labelled_molecules.py` | generic parquet → chunked ingest |
| `scripts/spiralia/export_parquet.py` | spiralia dill object → parquet |
| `scripts/spiralia/view_napari.py` | load one embryo in napari, per the source notebook |
| `scripts/spiralia/{PLAN,QUESTIONS}.md` | design decisions; open questions |
| `lib/webgl/labelled-molecule-{shaders,lut}.ts` | shader and lookup-table builders |
| `lib/stores/createLabelledMoleculeVisualizationStore.ts` | viewer state |
| `components/labelled-molecule-*.tsx` | scene, rail, legends, top controls |
| `lib/hooks/useLmVizUrlSync.ts` | URL state |
