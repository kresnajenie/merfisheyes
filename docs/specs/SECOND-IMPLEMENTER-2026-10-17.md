# Second-implementer test: a `labels` element kind from the written contract

**Date:** 2026-10-17. **Base:** `feat/format-specs` at `2ca5c8d`.
**Exercise:** add an OME-NGFF label image (`labels/<name>`, axes `y, x`,
`uint16`/`uint32`, 0 = background) as a fifth element kind, using only
`docs/PHILOSOPHY.md`, `docs/specs/*.md` and the named CLAUDE.md sections as
the design, and report where they were not enough. Nothing was asked of
anyone; every open point below was decided and is marked.

Outcome: the kind works end to end in the unit suites (opener on a real
Zarr v3 fixture, layer through the conformance suite and its own tests,
registry entry, controls), `tsc` is clean, the whole unit suite passes. Not
done: the URL slot (step 4), the `@quick` browser case (step 5, e2e was out
of scope), and a labels format profile / validator check (not in the five
steps at all). The one design rule I could not follow is the one the
exercise flagged: the tile machinery could not be reused by composition.

## (a) What was built, file by file

| File | What |
|---|---|
| `lib/spatialdata/labels.ts` (new, 104 lines) | `openSpatialDataLabels(store, name?)` → `SpatialDataLabels { name, bytesPerPixel, levels: LabelsLevel[] }`. Reads `labels/<name>`'s `attributes.ome.multiscales[0]`, requires axes exactly `y, x` and arrays `uint16` or `uint32`, composes each level's `coordinateTransformations` with the multiscale's, and undoes the cells' element transform the same way the image does. `LabelsLevel` is `ImageLevel` without the channel. Metadata only: one group open plus one array open per level. |
| `lib/spatialdata/image.ts` (modified, +27/−20) | Exported `Transform` and `composeTransforms`; extracted the "scene units → shared system from the shapes element, else points" block of `openSpatialDataImage` into an exported `cellsToShared({ store, keys })`. `openSpatialDataImage` calls it; behaviour unchanged (the three image unit tests and the registry tests pass as before). Done so the labels opener does not copy the transform rule. |
| `lib/spatialdata/LabelsLayer.ts` (new, 494 lines incl. comments) | `class LabelsLayer extends BaseLayer`, `kind = "labels"`. Typed API: `visible`, `opacity` (0–1, a uniform). `tick(frame)` reads `frame.view`: picks the finest level whose pixel is ≥ 0.75 screen pixel within 12 tiles, wants the coarsest level whole plus the chosen level's chunks in view, queues reads (4 in flight), one `AbortController` per tile with `hint(signal, { tag: "labels", priority: level })`, reads `array.getChunk([y, x], { signal })`, uploads a `DataTexture` of `RedIntegerFormat` / `R16UI` or `R32UI`, cancels tiles that left the view, evicts idle tiles over `min(256 MB, budget.residentBytes)`. One mesh per chunk; the fragment shader (GLSL ES 3) hashes the id (Knuth multiplicative, folded) into one of 24 `DEFAULT_COLOR_PALETTE` colours, `discard`s id 0, and writes `vec4(colour, uOpacity)`; the material is `transparent: false` with `CustomBlending` (src alpha / one-minus-src-alpha) so it blends in the opaque pass, `renderOrder = -999.9 - 0.01·level` (over the image at −1000−level, under the polygon fill at −999.5). `status()`, `abort()`, `dispose()`, the E2E hooks (`markLayerChanged` on a new wanted set, `markLayer`, `markLayerContent`, `markLayerSettled`). |
| `lib/scene/contract.ts` (modified) | `ElementKind` gains `"labels"`. |
| `lib/scene/registry.tsx` (modified, +16) | `OpenedElements.labels?: SpatialDataLabels`; `const labels: ElementDefinition<"labels">` with `detect` = `elementNames(keys, "labels")`, `open` = `openSpatialDataLabels(store, name)`, `createLayer` = `new LabelsLayer(opened)`, `Controls` = `SpatialDataLabelsControls`, `slot: "panel"`; inserted into `ELEMENTS` after `image` (draw order: image, labels, shapes, points, molecules). |
| `lib/scene/useLayers.ts` (modified) | `elements.labels` added to the effect's dependency list (see gap 1). |
| `components/spatialdata-labels-controls.tsx` (new, 76 lines) | "Labels" glass button beside the image button (`top-[4.25rem] left-[5rem]`) opening a glass panel: "Show labels" checkbox, opacity slider (5–100 %), a line with tiles drawn / wanted and a loading note. Re-renders through `useLayerStatus(layer)`; mutates the layer through `visible` / `opacity`. `data-testid="labels-button"` / `"labels-panel"`. |
| `tests/unit/layer-contract.test.ts` (modified, +42) | `fakeLabels()` (two levels, 256² and 64², ids in quadrants, 0 every eighth pixel, mirroring `fakeImage`'s geometry so the same single coarse tile loads) and a `labels` case in `CASES` whose change is `opacity = 0.3`. |
| `tests/unit/registry.test.ts` (modified, +58) | `fixtureStore` takes `labels?: boolean` and routes `labels/…` keys to `labels_v1/`; the `ELEMENTS` order assertion now includes `labels`; a new case opens `{ image, shapes, labels }` and asserts the name, `bytesPerPixel`, two levels, that every level's size, chunking, pixel size and origin equal the image's (both are placed through `cellsToShared` with the shapes' identity transform), and that chunk `[0, 0]` of level 0 holds exactly the ids `labels_expected.json` says. |
| `tests/unit/labels-layer.test.ts` (new, 166 lines) | The layer on a fake one-level uint32 image: four tiles read with an `AbortSignal`, drawn as `R32UI` textures in opaque-pass blending at a render order between the image and the fill; opacity and visibility change without a re-read; `abort()` aborts every signal and leaves the phase valid; `dispose()` frees tiles and bytes. |
| `scripts/testdata/make-tiny-labels.py` (new, 105 lines) | Writes `tests/data/tiny/zarr/labels_v1/cell_labels/` with numpy + zarr-python: Zarr v3, OME-NGFF 0.5, axes `y, x`, `uint16`, 64×64 and 32×32 (strided, never averaged), chunks 32×32, the exact `scale` / `translation` of `images/dapi`'s levels, an `image-label` block naming `../../images/dapi`. Ids 1–30 are the tiny table's cells in order (same RNG seed and draw order as `make-tiny-zarr.py`), each owning the pixels within 4 units of its centre, nearest first; 608 labelled pixels. Also `labels_v1/labels_expected.json` (ids present, pixel count, ids in the first chunk). 72 KB. |
| `docs/specs/SECOND-IMPLEMENTER-2026-10-17.md` | This report. |

Not touched: `ImageLayer`, `PointsLayer`, `ShapesLayer`, `MoleculeLayer`,
any scene or viewer, `url-state-writer.ts`, the e2e helpers and specs, the
validator, the Python spec tests, the other docs.

## (b) Where the written docs were insufficient or wrong

Numbered in the order I hit them. "Guessed" = decided from the docs' spirit;
"code only" = the answer exists only in source; "contradicts" = the docs say
one thing and the code or the task another.

### Step 1, open it

1. **The transform rule has no exported implementation** (code only).
   `store-layout.md` §8 and `image-profile.md` state "undo the shapes
   element's transform, else the points element's", but in code that lived
   as a private block inside `openSpatialDataImage`, with `composeTransforms`
   and `Transform` module-private. A second kind either copies ~35 lines or
   edits `image.ts`. I edited: exported both and extracted `cellsToShared()`.
   The spec should name the helper a kind placed in the shared system must
   call.
2. **No labels profile** (guessed). `image-profile.md` is written for
   `c, y, x` only and nothing says what a label image must be. I decided:
   axes exactly `y, x`; `uint16` or `uint32`; the multiscale and dataset
   transforms composed as for an image; placed by `cellsToShared`; the OME
   `image-label` block (colors, properties, source) read by nothing. A
   `labels-profile.md` with those rows is the missing document.
3. **`store-layout.md` is now wrong**: `labels/<image>/ (read by nothing
   yet)`. The CLAUDE.md "Element registry" section also lists the
   `OpenedElements` keys literally (image, shapes, points, molecules).
4. **The format side is absent from the five steps.** PHILOSOPHY §3 and the
   "Formats as specs" section say "change a format by changing its spec,
   schema and validator together", but "Adding an element kind" never
   mentions a profile, `scripts/validate_store.py` or
   `tests/python/test_specs.py`. I ran the validator on a store composed
   with the labels fixture: it passes and never mentions `labels/` — a
   labels element with the wrong axes would pass too. The Python spec
   tests compose points and shapes into the tiny store by hand; labels
   would need the same hand edit.
5. **`open(store, name, ctx)`: is `name` to be trusted?** (code only). The
   definition passes a `name`, but `openSpatialDataImage(store)` ignores it
   and re-detects. I made `openSpatialDataLabels(store, name = first)`.

### Step 2, draw it

6. **Composition of `ImageLayer`'s tile machinery is not possible** (code
   only; the finding the exercise asked about). Everything that would be
   reused — the level choice, `chunksInView`, the `Tile` record, the queue,
   `pump`, `loadChunk`, `makeTexture`, `evict`, `report` — is `private`, is
   written for N windowed channels (`textures[]`, `loading[]`,
   `controllers[]`, the 99.9th-percentile auto-window on the coarsest tile),
   indexes chunks as `[c, y, x]`, and bakes one fragment shader (additive,
   windowed channels). The public surface is `visible`, `setChannel`,
   `hasWindow`, `channels`, `maxValue`. The lifecycle half *can* be composed
   (an inner `ImageLayer` attached under a sub-`SceneHandle` whose `group`
   is the outer layer's root; `BaseLayer` needs only `group` and
   `invalidate`), and the reads could be routed through a fake 3-D array
   whose `getChunk([c, y, x])` returns the R, G or B byte of a hashed
   colour as three `uint8` "channels" — but that cannot give opacity, a
   transparent background, or a render order above the image, all three of
   which the kind needs. So `LabelsLayer` re-implements the same rules for
   one channel: about 200 of its lines are `ImageLayer`'s structure with
   the channel arrays removed. The fix is a `TileLayer` base (level choice,
   tiles in view, queue, eviction, hooks) with `loadTile` / `makeMaterial`
   as the per-kind hooks, extracted from `ImageLayer` — which this
   exercise's constraints put out of reach. Until then "one implementation
   per concern" is violated by any second tiled kind.
7. **"Blended in the opaque pass" has no recipe in `RENDERING.md`** (code
   only). The table says the fill is "blended in the opaque pass at −999.5"
   and why; the how (`transparent: false`, `blending: CustomBlending`,
   `SrcAlphaFactor` / `OneMinusSrcAlphaFactor`, `depthTest/Write: false`)
   is a comment in `ShapesLayer.ts`. A layer with an opacity must know it.
8. **No render-order rule for a new kind** (guessed). `ELEMENTS` is "in
   draw order" and `RENDERING.md` gives the existing numbers; nothing says
   where a new kind slots or that three.js `renderOrder` on a group does
   not propagate to its children (both `ImageLayer` and I set it on every
   mesh). I chose `-999.9 - 0.01·level`: above every image level, below
   the fill.
9. **The E2E hook semantics are in `ImageLayer.report()` only** (code
   only). Step 2 says "write the E2E hooks (`markLayerChanged/Content/
   Settled`, `markLayer`)". When "changed" fires (the set of wanted tiles at
   the chosen level differs from last frame's), what "settled" means (every
   wanted tile at that level visible) and that the counters must be in the
   layer's own unit are not written down; `BUDGETS.md` "How the numbers
   are read" describes the fields, not the obligations. The kinds listed
   in `test-hooks.ts`'s comment (`"molecules-legacy"`) do not match the
   registry's (`"molecules"`).
10. **`status().bytes` / `requests`: "when it can"** (guessed). The store
    counts wire bytes per tag; `ImageLayer` reports 0 / 0. I count
    requests (one per chunk read) and report 0 bytes, since decoded bytes
    are not network bytes. The spec should say either "the store's `byTag`
    is the source; report 0" or "report what you issued".
11. **`setBudget({ bytesPerView })` for a tile layer** (guessed).
    `LOADING.md` describes it for points. `ImageLayer` ignores it; so do I.
    `residentBytes` I honour exactly as the image does.
12. **A disposed layer must still answer `status()`** (code only). The
    contract test calls `status()` after `dispose()` in my own test; the
    spec's rule list says dispose "empties the group, drops subscribers",
    not that `status()` stays callable. It does, trivially, but it is a
    rule the suite relies on.

### Step 3, register it

13. **"No viewer or scene changes" is not quite true** (contradicts).
    `lib/scene/useLayers.ts` lists the kinds in its effect's dependency
    array (`[elements.image, elements.shapes, elements.points,
    elements.molecules]`). Without adding `elements.labels`, a labels
    element arriving on its own never creates its layer, and nothing
    fails loudly. Found only because I read `useLayers.ts` as one of the
    allowed entry points. The spec should either name this file in step 3
    or the hook should derive the list from `ELEMENTS`.
14. **A new kind must edit an existing assertion** (contradicts the
    task's "do not change existing tests except to add your cases").
    `registry.test.ts` asserts `ELEMENTS.map(kind)` equals a literal
    list. I changed it. Either the spec says "and update the order
    assertion" or the test checks uniqueness and the presence of an
    opener per kind instead of the exact list.
15. **A labels-only store shows nothing** (code only). `spatialdata-
    viewer.tsx` chooses the scene by cells / molecules / image; a store
    holding only a label image has no scene. Out of scope here, but the
    registry's promise that a viewer "only decides which scene to show"
    hides that the scene choice itself enumerates kinds.

### Step 4, control it

16. **`slot: "panel"` says nothing about placement** (code only / guessed).
    Panel controls position themselves absolutely: image at
    `top-[4.25rem] left-4`, the points status at `right-6`. A new panel
    has to read the others and find a gap. Labels went to `left-[5rem]`,
    beside the image button; the image panel covers it while open. A
    stacking container for the slot would remove the guesswork.
17. **The URL slot cannot be added without editing the writer**
    (contradicts "adding a kind is adding an entry"). Step 4 says keep
    settings in a URL slot through `lib/utils/url-state-writer.ts`; that
    file has six hard-coded pending slots (`v`, `rv`, `ov`, `rov`, `img`,
    `rimg`) and no generic one, and `imageUrlState.ts` is image-specific.
    I did not add a `lbl=` slot; labels settings are not in the link.

### Step 5, prove it

18. **Two hand-written composers of the tiny store** (code only).
    `tests/unit/registry.test.ts`'s `fixtureStore` and
    `e2e/helpers/spatialdata.ts` each merge points and shapes into the
    tiny store by key prefix; a new fixture folder is added to both by
    hand. I added it to the unit one only (letters `c m M p s i` in the
    server have no `l`), and `e2e/tests/viewer/elements.spec.ts` iterates a
    literal `["image", "points", "shapes"]`, so "asserts every kind reaches
    its layer" will not know about labels until edited. The spec should
    name both places, or they should share one composer.
19. **The conformance suite's hidden geometry** (code only). The "pushes a
    change once" case shows the layer, ticks, flushes, then expects the
    change to notify exactly once. Any read the change triggers would
    notify again. The camera (z = 100, 75°, 800 × 600) and
    `MAX_TILES_PER_VIEW` make the image fake fall back to its one coarse
    tile; I mirrored `fakeImage`'s sizes to get the same. Nothing says a
    fake must settle within one flush, or that "change" must be a pure
    appearance change.
20. **What a "fake" must provide** is only visible from `fakeImage`: an
    `array` with `shape` and `getChunk` cast `as never`. The spec could
    say "the opener's handle with its reads stubbed".

### Process

21. `git checkout feat/format-specs` fails when the branch is checked out
    in another worktree (it was: `…/worktrees/loading-standardization`).
    I used `git reset --hard feat/format-specs` on the agent's own branch
    from a clean tree, which gives the same content.
22. `LOADING.md`'s "Who reads what" table and `BUDGETS.md` have no row for
    the new kind; the spec does not say to add one. (Labels: one
    controller per tile, `priority = level`, tag `labels`.)
23. Not verified on a GPU: the shader (`usampler2D`, `uint` multiply and
    shift, `discard`, `R32UI`) and the blending are checked only by unit
    tests on the material's settings. A browser run was out of scope.

## (c) The five steps: effort and difficulty

| Step | Time | Why |
|---|---|---|
| 1 Open it | short–medium | The image opener is a good model; the time went to the private transform helpers (gap 1) and deciding the profile (gap 2). |
| 2 Draw it | **long, the hardest** | Composition could not meet the kind's needs (gap 6); the layer had to be re-derived from `ImageLayer` rule by rule, and the opaque-pass blending and render-order rules found in other layers' comments (gaps 7–8). |
| 3 Register it | short | One entry — plus two hidden enumerations (gaps 13–14) that took longer to find than the entry took to write. |
| 4 Control it | short | `useLayerStatus` and the glass primitives made the panel quick; the URL slot was not attempted (gap 17); placement was guessed (gap 16). |
| 5 Prove it | medium | The fixture script was straightforward (the image writer is a clear model); the registry fixture composer and the contract fake needed code reading (gaps 18–20). |

## (d) Commands run and their results

All from `/home/kjenie/merfisheyes/.claude/worktrees/agent-a7bcee6674caff70a`.

```
git checkout feat/format-specs
  fatal: 'feat/format-specs' is already checked out at '…/worktrees/loading-standardization'
git reset --hard feat/format-specs
  HEAD is now at 2ca5c8d feat(specs): the formats as specifications, with a validator and tests

python3 scripts/testdata/make-tiny-labels.py
  wrote 30 cell ids over 608 pixels to tests/data/tiny/zarr/labels_v1/cell_labels
  (9 files, 72 KB: zarr.json, 0/zarr.json, 0/c/{0,1}/{0,1}, 1/zarr.json, 1/c/0/0, labels_expected.json)

npx vitest run tests/unit/layer-contract.test.ts tests/unit/registry.test.ts
  ✓ tests/unit/layer-contract.test.ts (20 tests)   — 5 kinds × 4 rules
  ✓ tests/unit/registry.test.ts (8 tests)          — incl. "opens the label image of a store, placed where its image is"
  Test Files 2 passed, Tests 28 passed

npx tsc --noEmit -p tsconfig.json
  (no output; 28 s; no .next-e2e/dev/types folder existed)

npx vitest run tests/unit/labels-layer.test.ts
  ✓ 4 tests passed

npx vitest run
  Test Files 32 passed (32); Tests 282 passed | 1 skipped (283)

npx eslint lib/spatialdata/labels.ts lib/spatialdata/LabelsLayer.ts lib/spatialdata/image.ts \
  lib/scene/registry.tsx lib/scene/useLayers.ts lib/scene/contract.ts \
  components/spatialdata-labels-controls.tsx tests/unit/*.test.ts   (no --fix)
  0 errors; 3 warnings (no-console, the same pattern as ImageLayer / image.ts / registry.tsx);
  test files ignored by the project's ignore pattern

python3 <scratchpad>/compose_labels_store.py     # tiny store + labels_v1 + consolidate
  zarr-python reads labels: (64, 64) uint16 [y, x]
  consolidated names labels nodes: ['labels', 'labels/cell_labels', 'labels/cell_labels/0', 'labels/cell_labels/1']
python3 scripts/validate_store.py <scratchpad>/composed/store.zarr
  OK root / tables/table / images/dapi …   exit 0   — labels/ not mentioned (gap 4)
```

Not run, by instruction: the e2e suites (`@quick`, `@upload`, `@perf`),
`npm run lint` (would auto-fix), `tests/python/test_specs.py`.
