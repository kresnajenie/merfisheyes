# Loading and querying

How every reader gets its bytes, and the rules that bound what a view
costs. The decision record for step 3 of the loading plan; the budgets it
serves are in [BUDGETS.md](BUDGETS.md), the findings it answers in
[AUDIT-2026-10-17.md](AUDIT-2026-10-17.md) §1.

## Decision

Every reader — a SpatialData table, a points element, polygons, an image,
a legacy molecule folder — reads through one `ElementStore`
([lib/storage/ElementStore.ts](../../lib/storage/ElementStore.ts)). The
store is a zarrita `AsyncReadable` (`get`, `getRange`) plus `open` for a
streamed body, over a `Backend` that knows only how to fetch: a public URL,
a dataset in the app's bucket through presigned URLs, or files dropped in
the browser ([lib/storage/backends.ts](../../lib/storage/backends.ts)).
Everything that governs loading lives in the store, once:

| Concern | Rule |
|---|---|
| Cache | One byte cache per store, least recently used out first, 256 MB by default (`cacheBytes`); a read in flight is shared by everyone who asks for it. Missing objects (403/404) are remembered too. |
| Cancellation | Every read takes an `AbortSignal` that reaches `fetch`. A shared read is aborted when its last waiter has gone; a read aborted while queued never starts. |
| Ranges | Ranged reads of one object issued within a tick are merged into one request when they lie within `mergeGap` (64 KB) of each other and the merged span is under `mergeMax` (8 MB); each asker gets its slice. Suffix reads (shard indexes) are their own cached reads. |
| Concurrency | One priority queue per store, `concurrency` (32) requests in flight. Higher priority first; reads without a signal (metadata) before everything. Browsers cap connections per host below this, so the queue mostly decides the order reads reach the browser; a lower bound costs: the medium fixture's table opens in 1,092 small ranged reads, 11 s at 8 in flight, 5 s at 32. |
| Metadata | A Zarr v3 store's root `zarr.json` is read once; its consolidated metadata is seeded into the cache, so opening any node costs no request. `resolveZarrUrl` hands the root it read to the opener. |
| Accounting | `stats()`: requests, bytes, in flight, queued, cache hits, merged, aborted, cache bytes, and the same per reader (`byTag`). |

zarrita forwards only the signal of a read to the store, so a reader's
hints ride on the signal: `hint(signal, { priority, tag })`, or
`hinted({ … })` for a reader with nothing to cancel.

## Who reads what

| Reader | Reads | Signal |
|---|---|---|
| `SpatialDataTableAdapter` | coordinates, obs columns, one gene's CSC slice; a gene asked twice at once is read once | — |
| `points.ts` | index slices (`priority: Infinity`, tag `points`); rows of a piece (tag `points:<level>`, so bytes per row are known per level) | the layer's, per piece |
| `PointsLayer` | one `AbortController` per piece, `priority = level` (coarse first) | aborted when the plan no longer wants the piece, on `abort()`, on dispose |
| `ImageLayer` | one controller per tile and channel, `priority = level` | aborted when the tile leaves the view, on `abort()` |
| `shapes.ts` | the element's GeoParquet, whole (tag `shapes`) | — |
| legacy molecules (`lib/molecules/storeLoading.ts`) | the manifest (`get`), a gene whole or streamed (`open`), its count from the gzip trailer (`getRange`, 4 bytes) | the layer's, per gene |

## Budgets in the readers

- **Points per view and bytes per view.** `PointsLayer` picks, per visible
  gene, the finest level whose molecules in view fit the gene's share of
  the point budget *and* whose bytes fit its share of the byte budget
  (`DETAIL_BUDGET`, `DETAIL_BYTES`: Auto 8 M and 12 MB, More 16 M and
  36 MB, All 64 M and unbounded). Bytes per row come from what the store
  has read at that level so far, with a prior of 3 B before any read.
  A scene's `setBudget({ bytesPerView })` caps Auto further.
- **Coarse first.** A gene with nothing on screen yet queues its coarsest
  level's pieces for the view first, as stand-ins: a few small reads, so
  something is drawn in the first moments while the chosen level streams in.
- **Resident memory.** `setBudget({ residentBytes })` bounds what a layer
  keeps for views that scrolled away (points at 17 B each, image textures);
  `status().residentBytes` reports it. The store's cache is bounded on its own.

## What is not done here

- The h5ad-zarr cell path (`StandardizedDataset.fromS3Zarr`, zarrita 0.5.1)
  still reads through `PresignedFetchStore`; it is not a SpatialData reader.
- Image textures keep their CPU copy: the scene rebuilds its WebGL context on
  a view-mode switch and layers now outlive that, so textures must be able
  to upload again.
- A store opened by dataset id lists one level of the bucket to find the
  store's root (`/list?delimiter=/`) and only falls back to the whole key
  list when the root is deeper than one folder.
- Hidden genes keep what they loaded and get no share of the budgets.

## Verification

- `tests/unit/element-store.test.ts`: cache and sharing, cancellation to the
  backend, merged ranges, priority and bound, per-tag counts, consolidated
  metadata, streaming.
- `tests/unit/molecule-store-loading.test.ts`: a molecule folder through the
  store — manifest, whole load with progress, streamed batches, trailer
  count, missing unassigned file, cancellation.
- `tests/unit/registry.test.ts`, `points-zarr.test.ts`,
  `spatialdata-table-adapter.test.ts`, `local-bundle.test.ts`: the readers
  on the tiny stores.
- `e2e/perf/elements.spec.ts` (`@perf`): the budgets on the tiny and the
  generated medium store at 12 MB/s; the numbers in BUDGETS.md.
