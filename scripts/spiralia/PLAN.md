# Labelled single-molecule viewer — implementation plan

A third dataset type alongside `single_cell` and `single_molecule`. Points are
molecules; every molecule carries N categorical labels and no expression matrix.
First dataset: `MER6-2_E3_1`, 3.15 M molecules × (gene, domain, cell).

## Decisions already locked

| | |
|---|---|
| format | molecules × categorical columns, **no expression matrix** |
| columns (E3_1) | `gene` (213), `domain_anno` (12), `domain_id` (18), `cell` (26) |
| menus | **three buttons** — Gene / Domain / Cell. The domain menu toggles anno ↔ id; switching resets that menu's selection, since the value sets differ |
| colour | switchable per button; **default = cell** |
| non-matching | grey at alpha 0.2, with a toggle for alpha 0 (hidden). Literal grey, not a desaturated palette colour |
| empty menu | no constraint from that menu |
| combine | OR within a menu, AND across menus |
| excluded at ingest | `cell_id == 0`, blank control probes |
| route | new, so nothing existing changes |
| coordinates | **raw µm**, `normalized: false` — keeps molecules co-registered with the µm mesh vertices, no transform needed when segmentation lands |
| gene colours | assigned on selection from the lowest free slot, SM-style, so selected genes are always distinct. Domain (12–18) and cell (26) fit a fixed palette and get one at ingest |
| cell menu | keyed on `cell_id`, displaying `cell_name`, suffixed when names collide (`pb (1)`, `pb (3)`) |
| camera | reuse SC `viewMode` 2D/3D |
| sizing | global size + alpha, plus per-value size override on the colour-by column |
| meshes | later |
| /explore listing | later — direct link only for now |
| plot panel | later |

## Phase 0 — naming and registration

- `datasetType = "labelled_single_molecule"`, route `/lm-viewer/[id]`.
- Prisma enum + migration.

**Verify:** migration applies; a row with the new type round-trips.

## Phase 1 — ingestion (Python)

New `scripts/process_labelled_molecules.py`: parquet → chunked folder.
Generic in the columns — takes `--x/--y/--z` and `--category` flags, so it is not
spiralia-specific. The existing `scripts/spiralia/export_parquet.py` stays as the
dataset-specific step that produces its input.

```
datasets/{id}/
├── manifest.json.gz          # n_molecules, columns, dimensions, has_expression: false
├── coords/spatial.bin.gz     # Float32 flat [x,y,z,…], raw µm
├── obs/
│   ├── metadata.json         # column list + types + unique counts
│   ├── gene.bin.gz           # dictionary-encoded codes
│   ├── domain_anno.bin.gz
│   ├── domain_id.bin.gz
│   └── cell.bin.gz
└── palettes/{col}.json
```

Reuse `encode_obs_binary()` and `build_obs_palette()` from
`process_spatial_data.py` — the binary layout and the adapter's reader already
match. No `expr/` directory is written.

**Verify:** decode every `obs/*.bin.gz` back to labels and assert equality with
the parquet column, row for row — the same check that passed on the 345 SM gene
files.

## Phase 2 — loading (TypeScript)

`ChunkedDataAdapter.initialize()` currently `Promise.all`s `expr/index.json`
(line 117), so it hard-fails without an expression matrix. Gate that fetch on
`manifest.has_expression`; leave every other path untouched.

`StandardizedDataset` then materialises with `genes: []` and 4 cluster columns
carrying `valueIndices` + `uniqueValues` + `palette`.

**Verify:** dataset loads from S3; `clusters` has 4 columns with unique counts
213 / 12 / 18 / 26 and `getPointCount() === 3_140_608`.

## Phase 3 — store

New `lib/stores/createLabelledMoleculeVisualizationStore.ts`:

- `colorBy: string` — the active column key, default `cell`
- `selections: Record<string, Set<string>>` — one set per column
- `paletteOverrides` / `sizeOverrides`: `Record<string, Record<string, …>>`
- global `sizeScale`, `alphaScale`, `unselectedAlpha` (0.2), `viewMode`

The visible mask is **not** stored — it is computed in the shader.

**Verify:** unit test that toggling values produces the expected selection sets;
URL round-trip.

## Phase 4 — scene and shader

New `components/labelled-molecule-three-scene.tsx` and a new shader pair.

Static vertex attributes, uploaded once:

- `position` — Float32, 3 × N
- `aCol0…aColN` — the category index per column. `Uint16Array` with
  `normalized = false` converts to an exact float in the shader, so 4 columns
  cost 25 MB rather than 50 MB as Float32.
- `aSize`, `aAlpha` — per-point, static for now, reserved for continuous
  per-molecule values (`correlation`, `brightness`, `scoreA` from `Xh`).

Small textures, re-uploaded on interaction:

- selection LUT per column (R8, ≤293 px)
- palette LUT for the colour-by column (RGB)
- per-value size LUT for the colour-by column (R32F)

Uniforms: `uColorBy`, `uGlobalSize`, `uGlobalAlpha`, `uUnselectedAlpha`.

Vertex shader: `visible = selCol0[aCol0] * selCol1[aCol1] * …`, colour from the
palette LUT at the colour-by index, then the size/alpha combination above. The
existing distance-fade, sub-pixel cull and shape logic in `lib/webgl/shaders.ts`
carries over unchanged.

**Verify:** toggling a checkbox uploads < 1 KB and shows no frame hitch;
rendered visible count matches the Python ground truth for a fixed set of
selections (e.g. gene ∈ {NANO2, ACTC} ∧ domain ∈ {3q, 3Q} ∧ cell ∈ {3B, 3D}
→ **32,665**).

## Phase 5 — UI

- `labelled-molecule-controls.tsx` — one button per column, opens the panel.
- Panel — search, checkbox list with per-value counts, palette swatch with
  colour picker, per-value size slider on the colour-by column.
- Legends — active selections per column, clear-all.
- Tooltip — reuse `pinnedTooltipColumns`, pinned to all columns by default for
  this type. No new tooltip code.

**Verify:** the UI's displayed intersection count equals the Python number for
the same selections.

## Phase 6 — upload and registration

Output is SC-shaped, so `/api/datasets/initiate`'s pre-chunked path should carry
it with a `datasetType` flag rather than a new upload route.

**Verify:** `MER6-2_E3_1` uploads and opens end-to-end at `/lm-viewer/[id]`.

## Risks

**Overdraw.** 3.1 M grey points at alpha 0.2, blended every frame, is a lot of
fragments, and a small selection (32,665 visible against 3.1 M grey is a 1%
signal) may read as mud. The sub-pixel cull helps zoomed out, not zoomed in. The
alpha-0 toggle is the escape hatch; which of the two is the better default is a
judgement to make once the scene is on screen.

**Palette LUT size.** The colour-by texture is sized to the largest column
(213 genes). Fine as a 1D R8/RGB texture, but the per-value size LUT wants
R32F — worth confirming float texture support rather than assuming it.
