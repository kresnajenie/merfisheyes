# SC + SM in One Space — Non-`__all__` Mode (Future Work)

**Status:** `__all__` mode ships and works. The **non-`__all__` (per-sample) mode
still needs to be developed.** This doc is a pointer map, not a design spec.

## What this feature is

Render single-molecule (SM) data in the **same Three.js scene** as a single-cell
(SC) dataset, so molecules overlay the cells they belong to.

Which SM dataset a cell maps to is driven by a `mapping.json` fetched from the
SC dataset's custom S3 base URL. Its shape:

```ts
// lib/config/dataset-links.ts:1
interface DatasetLinkConfig {
  linkColumn: string;              // "__all__" OR an obs/cluster column name
  links: Record<string, string>;  // column value → SM S3 base URL
}
```

There are two modes, keyed by `linkColumn`:

- **`__all__`** — every cell links to **one shared** SM dataset. The SM overlay
  loads into the SC scene automatically. **Built and working today.**
- **non-`__all__`** — `linkColumn` is a real obs/cluster column (e.g. `"batch"`),
  so each column value (`set1`, `set2`, …) maps to a **different** SM dataset.
  See the registry example at `lib/config/dataset-links.ts:8` (`linkColumn: "batch"`,
  10 sets). **This is the per-sample case that is NOT yet built as an in-one-space
  overlay.**

## The vision for non-`__all__`

Because non-`__all__` is **per sample** (many SM datasets keyed off a column
value), it can't just auto-overlay one dataset like `__all__` does. It needs a
**new interactive feature** to let the user drive which sample(s) render in the
shared space — a purpose-built interaction, not the current behavior.

**Today's stopgap for non-`__all__`:** right-click a cell → opens the mapped SM
dataset in a **split panel** (separate scene), NOT an in-space overlay. That is
the interaction the new feature would replace/augment. See "Right-click split"
pointers below.

## Where everything is

### Mapping config
- `lib/config/dataset-links.ts` — `DatasetLinkConfig` type (`:1`), hardcoded
  `DATASET_LINK_REGISTRY` with the `batch` example (`:8`), `getDatasetLinkConfig()`
  (`:27`), and `fetchMappingConfig()` that pulls `mapping.json` from S3 (`:40`).

### The scene — `components/three-scene.tsx`
- SM overlay refs/state (`:119`–`:139`): `smDatasetRef`, `smPointCloudsRef`,
  `smDataset`, `smLoading`, plus SM gene selection read from the SM viz store.
- Link-config lifecycle effect (`:290`–`:368`): resolves registry → `mapping.json`,
  clears any stale overlay (`clearSmOverlay`, `:302`), and **the `__all__` load
  path** (`:339` `if (mappingConfig.linkColumn === "__all__")` → `fromCustomS3`,
  push into global SM store).
- **Right-click handler (the non-`__all__` path):** `:795`–`:835`.
  - `:807` `if (linkConfig.linkColumn === "__all__") return;` — `__all__` bails out
    (already overlaid).
  - `:809`–`:824` reads `linkColumn` from `dataset.clusters`, resolves the clicked
    cell's value → `linkConfig.links[setValue]`.
  - `:832`–`:834` `enableSplit()` + `setRightPanelS3(smUrl, "sm")` — **this is the
    current split-panel stopgap; the in-space per-sample overlay would hook in
    here** (or replace this branch).

### SC/SM overlay UI (only meaningful when an overlay exists)
- `components/sc-sm-mode-toggle.tsx` — `ScSmMode = "sc" | "sm"` toggle (`:8`).
- `components/single-molecule-gene-picker.tsx` — SM gene search/checkbox picker
  used on the SC overlay page.
- `components/visualization-controls.tsx:80`–`:89` — SC/SM layer visibility
  toggles; `hasSmOverlay` is gated on the **global SM store** having a dataset
  (`:89`), which today only happens for `__all__`.

### Stores
- `lib/stores/createSingleMoleculeVisualizationStore.ts` — overlay viz state:
  `selectedGenesLegend` (`:27`), `smLayerVisible` (`:35`). The picker + legends
  render off `selectedGenesLegend`.
- `lib/stores/singleMoleculeStore` (global) — holds the loaded overlay SM dataset;
  `three-scene.tsx` pushes into it (`:361`) and `visualization-controls.tsx` reads
  `currentDatasetId` from it to decide `hasSmOverlay`.

### URL persistence
- `app/viewer/from-s3/page.tsx:54`–`:64` — SM overlay URL state via a dedicated
  `ov=` slot; reads the loaded SM dataset from the global SM store (populated by
  the `__all__` path in `three-scene.tsx`).

## TODO for non-`__all__`

- Design + build the per-sample interactive feature (the "vision" above) —
  currently unbuilt; right-click split (`three-scene.tsx:832`) is the only
  non-`__all__` behavior.
- Decide how a per-sample overlay populates the global SM store so the existing
  overlay UI (`hasSmOverlay` at `visualization-controls.tsx:89`, gene picker,
  layer toggles) lights up for non-`__all__` too.
- Extend `ov=` URL persistence (`from-s3/page.tsx:54`) to encode which sample(s)
  are overlaid.

## History / context

- Original (now-removed) proof-of-concept plan: `PLAN-sc-sm-overlay.md` at repo
  root, added in commit `baafd41`, removed in `04fcffd`. `git show baafd41:PLAN-sc-sm-overlay.md`
  to read it.
- The `__all__` right-click split was intentionally dropped in `0924b62`
  ("Camera recenter on cell + reset; drop `__all__` right-click split").
