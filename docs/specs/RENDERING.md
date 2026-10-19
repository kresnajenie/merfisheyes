# Rendering

How what a layer loaded reaches the screen, and the rules that keep a frame
within its budget ([BUDGETS.md](BUDGETS.md): 16 ms while the camera moves,
a slider within 100 ms, hover never blocking a frame). The decision record
for step 4 of the loading plan; the findings it answers are in
[AUDIT-2026-10-17.md](AUDIT-2026-10-17.md) §2.

## Rules

| Rule | Where |
|---|---|
| **A still scene costs nothing.** One `requestAnimationFrame` loop per scene, drawing only when something asked for a frame (`invalidate()`, a layer's `changed()`, the controls) or at the idle rate. No layer has a loop. | `lib/webgl/scene-manager.ts`, `lib/scene/LayerHost.ts` (step 2) |
| **Constants are uniforms.** A molecule cloud carries positions only (12 B a molecule); colour, size, alpha and shape are uniforms of its material (`createMoleculeMaterial`), so a slider tick changes a uniform and uploads nothing. The cell cloud keeps colour, size and alpha per vertex because every cell differs. | `lib/webgl/shaders.ts` (`moleculeVertexShader`, `moleculeFragmentShader`), `lib/webgl/point-cloud.ts`, `PointsLayer`, `MoleculeLayer` |
| **A cell tick allocates no GPU buffer.** New colours, sizes and alphas are copied into the attribute's own array and re-uploaded in place; the gene's percentile is computed once per gene and cached. | `updatePointCloudAttributes`, `calculateGenePercentile` |
| **Hover is constant time.** A uniform grid over the cells' positions as drawn (rebuilt when they change) answers the cell under the pointer from the cells around it; a 2D dataset is picked where the pointer's ray meets the data plane, a 3D one by screen distance among the cells near that spot. Hidden cells are not picked. | `lib/webgl/picking-grid.ts`, `three-scene.tsx` |
| **Sizes are CSS pixels.** The shaders clamp a point to [1, 200] CSS pixels and scale by `uPixelRatio`; the renderer's pixel ratio is capped at 2. | `shaders.ts`, `scene-manager.ts` |
| **Spots and bins at their real size.** When the cells have a footprint (Visium spots, Visium HD bins: `lib/spatialdata/footprint.ts`), the cell cloud's `dotSize` is that size in world units (`uWorldSized`), turned into pixels by the canvas half height (`uHalfHeight`, pushed with the per-frame uniforms), drawn as circles or squares (`uShape`) and capped at 8192 CSS pixels instead of 200 (the GPU's own point-size limit, typically 1024 px or more, still applies when one spot fills the view). A cell at the "selected" size multiplier is drawn at its real size. | `shaders.ts`, `setCellFootprint`, `three-scene.tsx` |
| **The second matrix multiply is skipped** unless the distance filter is on. | `shaders.ts` |
| **Draw order is render order.** Image tiles under everything (−1000 − level), the polygon fill blended in the opaque pass at −999.5, molecules by legend position (−999 upwards), cells at 0, the polygon outline over them. A transparent fill would be drawn after every opaque object, over the molecules. | `ShapesLayer`, `PointsLayer`, `MoleculeLayer` |
| **Stand-ins are modest.** A coarser level standing in for a view still loading is drawn at most 2× its size. | `PointsLayer` |
| **High-performance GPU asked for.** `powerPreference: "high-performance"`. | `scene-manager.ts` |

## Not done here

- **Cells through a lookup texture** (expression and cluster code as scalar
  attributes, colour from a colormap texture, every slider uniform-only).
  The polygons take their colour from the cells' colour attribute
  (`SceneHandle.cellColors`) and the tooltips read it too, so the LUT needs
  both moved first; with attributes reused in place and the percentile
  cached, a tick costs one copy of 20 B per cell and no allocation.
- **Pooled draws** for the points element (one geometry per gene and level,
  pieces as draw ranges). At Auto the medium store draws 79 calls; the
  ten-section store 887 at the old budget. Pieces leave holes when evicted,
  so a pool needs compaction; deferred until the draw count is the measured
  bottleneck.
- **Image levels stacked.** Coarser tiles stay under the current level so a
  view refines instead of blanking; the extra quads are few.
- **`preserveDrawingBuffer`** stays on: the export box and the thumbnail
  capture read the canvas after compositing. Rendering once into an
  offscreen target for them is the step after.
- **Resize** is already coalesced by `ResizeObserver` and the on-demand
  loop (a resize asks for one frame).

## Verification

- `tests/unit/picking-grid.test.ts`: the grid agrees with a pass over every
  point; a screen-space distance; empty and degenerate sets.
- `tests/unit/streaming-point-cloud.test.ts`, `molecule-layer.test.ts`
  ("never refetches for an appearance change, but applies it"): a molecule
  cloud carries positions only and its colour is a uniform.
- `e2e/perf/elements.spec.ts`: `sliderToFrame`, `sliderToFrameHeavy`,
  `avgRenderMs`, `drawCalls`, `hoverFrameMs`; numbers in BUDGETS.md.
