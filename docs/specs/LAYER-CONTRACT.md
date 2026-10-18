# The layer contract

**Status:** decided and implemented 2026-10-17 (`lib/scene/contract.ts`).
Every element drawn in a scene goes through it; the conformance suite
`tests/unit/layer-contract.test.ts` is what "implements the contract" means.

## Decision

A **layer** draws one element kind (an image, a points element, a shapes
element, legacy molecule files) inside the group a scene gives it. All
layers share one lifecycle, one way of reporting, one way of joining a
scene, one frame. What a layer is told to draw stays typed per layer.

```ts
interface Layer {
  readonly kind: ElementKind;                 // "image" | "points" | "shapes" | "molecules"
  attach(scene: SceneHandle): void;           // join; again after a scene rebuild
  detach(): void;                             // leave, keeping what was loaded
  tick(frame: FrameInfo): void;               // once per frame, from the scene, before render
  status(): LayerStatus;                      // phase, progress, drawn, wanted, resident, bytes, requests
  subscribe(listener: () => void): () => void; // push; nothing polls a layer
  setBudget(budget: LayerBudget): void;       // bytes per view, bytes resident
  abort(): void;                              // cancel every download; what is drawn stays
  dispose(): void;                            // free everything; attach afterwards throws
}
```

The scene's side is `LayerHost` (`lib/scene/LayerHost.ts`): it builds the
`SceneHandle` (group, camera, renderer, `view()`, `invalidate()`,
`attach()`, `cellColors?`), computes the view rectangle **once per frame**
for every layer (`lib/scene/view-rect.ts`), ticks the layers from the
scene's single `requestAnimationFrame` loop (`initializeScene`'s `onFrame`),
and is where a layer joins or leaves. There is no other way in.

`BaseLayer` (`lib/scene/BaseLayer.ts`) implements the shared parts: the root
group that joins the scene's group, `subscribe` with notifications coalesced
to one per task, the budget, `changed()` — which asks the scene for a frame
and tells the subscribers — and `dispose`.

## Reasons

- Four layers had four constructor shapes, four update styles, status by
  callback / a 200 ms poll / toasts, their own `requestAnimationFrame`
  loops each computing the same view rectangle, and three ways of being
  attached (`docs/specs/AUDIT-2026-10-17.md` §3). A fifth kind would have
  added a fifth of each.
- One loop means one frame: layers see the same view, tick before the
  render, and the scene can render on demand once every change goes
  through `changed()`.
- Push status means a controls panel costs nothing while nothing happens.
- `abort()` and `setBudget()` are in the contract so that cancellation and
  the byte budgets (steps 3 and 4 of the plan) arrive for every layer at
  once, and the conformance suite can check them for every layer.

## What is typed per layer

How a layer is told what to draw is its own API, because the state differs
in kind:

| Layer | Drives it | From |
|---|---|---|
| `ImageLayer` | `visible`, `setChannel(i, patch)`, `hasWindow(i)` | the image controls, the `img=` URL slot |
| `PointsLayer` | `visible`, `setGenes(Map<gene, GeneStyle>)`, `setDetail()`; `pointsStatus()` for the panel | the SM viz store (`three-scene.tsx` effect 5a), the status panel |
| `ShapesLayer` | `visible`, `mode`, `opacity`; `shapesStatus`, `zoomedIn` | the shapes controls |
| `MoleculeLayer` | `sync(selection)` | the SM viz store (`useMoleculeLayerSync`) |

The registry entry's `Bind` component maps the panel's stores to these
calls, so a viewer does not know them.

## The registry

`lib/scene/registry.tsx` holds one `ElementDefinition` per kind, in draw
order (`ELEMENTS`): how the kind is found in a store (`detect`), opened
(`open`, metadata only), drawn (`createLayer`), what a scene must offer
(`canAttach`, e.g. the cells' colours for the polygons), the component that
keeps the layer in step with the panel's stores (`Bind`) and its controls
(`Controls`, in the viewer's button rail or as a floating panel: `slot`).

Everything above the entries is generic:

- `openElements(store, ctx)` iterates the registry and returns an
  `OpenedElements` record. The one dataset loader
  (`lib/viewer/load-dataset.ts`) calls it for a SpatialData store;
  `lib/local/open-bundle.ts` for a folder dropped for a preview.
- `useLayers(elements, ctx)` (`lib/scene/useLayers.ts`), called by the
  scene component, creates one layer per opened element, publishes them to
  the panel's layers store (`lib/stores/createLayersStore.ts`) and returns
  the `onSceneReady` that attaches them to each (re)built scene. Layers
  outlive scene rebuilds: a view-mode switch detaches and re-attaches them,
  so what they loaded and how they were set stay.
- `LayerControls slot` (`components/layer-controls.tsx`) renders the
  controls of every layer the panel's scene is drawing, in the viewer's rail
  (`VisualizationControls`) or inside the scene; `LayerBindings` renders
  every `Bind`.

A viewer therefore only decides which scene to show and passes
`elements={…}` to it; `three-scene.tsx` adds the overlays it still resolves
itself (a `mapping.json`, a dropped folder) to the same record.

## Adding an element kind, in five steps

1. **Open it.** A function that turns a store (or files) into the element's
   handle, reading only metadata (`lib/spatialdata/image.ts`,
   `points.ts`, `shapes.ts` are the models). Budget: open in well under a
   second; everything else lazily, by view.
2. **Draw it.** `class XLayer extends BaseLayer` with `kind`, a typed API for
   what to draw, `tick(frame)` reading `frame.view`, `status()`, `abort()`,
   `dispose()` (calling `super`). Call `this.changed()` whenever what is drawn
   or reported changes, and never run a loop of your own. Write the E2E
   hooks (`markLayerChanged/Content/Settled`, `markLayer`) so the budgets
   can be measured on it.
3. **Register it.** Add the kind to `ElementKind` and `OpenedElements`
   (`lib/scene/contract.ts`, `registry.tsx`) and one `ElementDefinition` to
   `ELEMENTS`: `detect`, `open`, `createLayer`, and a `Bind` component when
   the layer follows a store. No viewer or scene changes.
4. **Control it.** A `Controls` component on the entry that uses
   `useLayerStatus(layer)` to re-render, mutates the layer through its typed
   API, and keeps its settings in a URL slot through
   `lib/utils/url-state-writer.ts`; `slot: "rail"` for a button in the
   viewer's rail, `"panel"` for a floating panel.
5. **Prove it.** Add a fake and a "change" to `CASES` in
   `tests/unit/layer-contract.test.ts`; a case to `tests/unit/registry.test.ts`
   (what `openElements` finds); the kind's own unit tests; and a `@quick`
   browser case on a fixture (`e2e/helpers/spatialdata.ts` serves stores
   with any combination of parts; `e2e/tests/viewer/elements.spec.ts` asserts
   every kind reaches its layer).

## Rules the suite enforces

- `kind` is set; `status()` is well-formed before and after attach.
- Attaching adds exactly one child to the scene group; detaching removes
  it; attaching again (a rebuild) works; a tick before attach is harmless.
- A change notifies subscribers once per task and asks for a frame;
  unsubscribing stops it.
- `setBudget` is accepted; `abort()` leaves the layer in a valid phase;
  `dispose()` empties the group, drops subscribers, and a later attach
  throws.
