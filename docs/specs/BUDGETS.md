# Budgets: human time, and what it implies for data and GPU

The acceptance criteria every loading, query and rendering change is judged
by. They come from how people perceive delay, not from what the code
happens to do today. The perf suite measures them (`e2e/perf/elements.spec.ts`
on fixtures at a capped 12 MB/s; `scripts/perf/real-store.ts` on the real
stores), and the "today" columns below are filled from those runs.

## Human limits

| Limit | Meaning | Source |
|---|---|---|
| 16 ms | one frame at 60 Hz while panning or zooming; 33 ms is the floor on integrated GPUs | RAIL animation budget (Google, 2015) |
| 100 ms | an interaction must show its result or a visible start | Nielsen 0.1 s; RAIL response |
| 500 ms | added latency at this level measurably reduces exploration: fewer observations, less coverage, and the effect carries over | Liu & Heer, *The Effects of Interactive Latency on Exploratory Visual Analysis*, InfoVis 2014 |
| 1 s | the flow of thought stays uninterrupted; first content of a new view | Nielsen 1.0 s (Miller 1968; Card, Robertson & Mackinlay 1991) |
| 10 s | the attention limit; anything longer needs a percent-done indicator and a way to cancel | Nielsen 10 s; Myers 1985 for percent-done |

Point-cloud viewers settle on a per-view point budget for the same reason:
Potree's default is 1 M points per frame, scaled to the GPU.

## Targets

Measured at the 95th percentile on a 100 Mbit/s connection (12 MB/s) and a
laptop GPU, on the 11.24 B-molecule, ten-section store
(`spatialdata/12_16_2025_all_sets/`) with three genes selected.

| Event | Target | What it needs |
|---|---|---|
| Open store → cells on screen | ≤ 2 s | one metadata fetch; index arrays read lazily |
| Pan/zoom → every frame | ≤ 16 ms | no per-frame CPU work, few draws, render on demand |
| Pan/zoom → first new molecules | ≤ 300 ms | stale reads cancelled, coarse level first, small pieces |
| Pan/zoom → view settled at "Auto" | ≤ 1 s | ≤ 12 MB per view: a byte budget beside the point budget, 2–3 B per molecule on the wire |
| Gene selected → first molecules | ≤ 300 ms | per-gene index slice (≤ 50 KB), then the coarse level |
| Slider, colour, toggle → next frame | ≤ 100 ms, no refetch, no re-upload of positions | constants as uniforms, colour through a lookup texture |
| Hover → frame | ≤ 16 ms, never blocks a frame | a spatial index, only when the camera is still |
| Any download | cancellable within 100 ms; > 10 s shows percent-done | an AbortSignal through every store |
| Memory | ≤ 1.5 GB JS heap and ≤ 1.5 GB GPU, both reported | one budget, CPU copies released after upload |

## Limitations, stated plainly

1. **Bandwidth bounds everything.** 12 MB/s is about 4 M molecules per
   second at 3 bytes each. "All" (64 M molecules per view) cannot be
   seamless over the internet and is labelled so.
2. **A tab has about 4 GB; a GPU has 1–2 GB to spare.** 64 M molecules at
   17 B are already 1.1 GB; "All" also exceeds the memory budget.
3. **The public bucket allows only the production and standard dev ports in
   CORS, and lists nothing.** Real-store measurements run from an allowed
   origin or through `scripts/perf/cors-proxy.cjs`.
4. **`spatialdata` 0.8 cannot read the points element.** Until the spec
   lands upstream, compatibility means the table and the image round-trip.
5. **3D stores** cost a third more bytes per molecule, and the tile grid is
   2D; the ten-section store is 2D by decision.

## Today (2026-10-17, before any change; step 3's numbers follow below)

Fixtures at 12 MB/s, headless Chromium with SwiftShader
(`npm run test:perf`, `perf/results/latest.json`). Tiny: 30 cells, 3,000
molecules. Medium: 50,000 cells, 30 M molecules, the heaviest gene 18 M.

| Metric | spatialdata-tiny | spatialdata-medium |
|---|---|---|
| openToCells | 1.6 s | 2.7 s |
| firstMolecules (default genes) | 1.9 s | 3.6 s |
| geneToSettled (18 M-molecule gene, whole section → coarser level) | — | 17.8 s |
| geneBytesKB / geneRequests | — | 24.0 MB / 148 |
| viewSettled (zoom to a fifth of the width, one wheel event) | 0 (all resident) | **132 s** |
| bytesPerViewKB / requestsPerView | 0 | 59.4 MB / 364 |
| hoverFrameMs (light scene) | 48 ms | 56 ms |
| sliderToFrame (light) | 82 ms | 808 ms |
| sliderToFrameHeavy (heavy gene drawn) | — | **19–32 s** |
| frameOver16Pct / worstFrameMs (SwiftShader, 1 s pan) | 99 % / 26 ms | 100 % / 1.4–5 s |
| avgRenderMs (CPU inside render) | 0.1 ms | 1.6–2.5 ms |
| drawCalls | 7 | 282 |

Real stores through `scripts/perf/cors-proxy.cjs` on this machine's link
(about 1 MB/s to the bucket, so bytes and request counts are the durable
numbers; times are upper bounds), headless Chromium with SwiftShader
(`perf/real/*.json`):

| Metric | all_sets (11.24 B) | set1 (167 M) | atera (738 M) |
|---|---|---|---|
| openToCells | 3.4 s | 2.5 s | 3.0 s |
| openBytesKB / openRequests | 13.7 MB / 76 | 361 KB / 37 | 8.6 MB / 58 |
| gene selected (Slc17a7 / Slc17a7 / EPCAM) → first molecules | 2.3 s | 0.30 s | 0.13 s |
| gene → settled view | **494 s** | 31 s | 4.2 s |
| gene bytes | 54.4 MB | 27.1 MB | 3.2 MB |
| molecules drawn / resident after the gene | 11.0 M / 22.1 M | 4.5 M / 4.5 M | 0.43 M / 0.43 M |
| zoom to ~1.5 mm (one wheel event) → settled | **266 s** | 0 (all resident) | 0 (all resident) |
| zoom bytes / requests | 86.8 MB / 640 | 0 | 0 |
| molecules drawn / resident after the zoom | 18.6 M / 30.0 M | — | — |
| drawCalls | 887 | 185 | 254 |
| worstFrameMs (SwiftShader) | 4.7 s | 1.3 s | 0.22 s |
| jsHeapMB (after the zoom on all_sets) | 242 | 61 | 82 |

Note on the zoom metric: a 12-notch zoom never made the planner choose a
finer level on any store (`viewSettled` 0, no bytes), so the suite zooms
32 notches (a fifth of the width) and the real-store script 64 (about
1.5 mm from a 40 mm atlas view). Done as a hand does it, in several wheel
events, **the zoomed view did not settle at all**: not within 120 s on the
medium store (4 events) nor within 600 s on the ten-section store
(7 events), because every intermediate view's pieces stay in the queue and
are fetched, coarse first, before the final view's (reads are never
cancelled — `docs/specs/AUDIT-2026-10-17.md` §1). The metric is therefore
taken on a single wheel event (one replan): 132 s and 59 MB on the medium
store, 266 s and 87 MB on the atlas, against a 1 s target. The stepped zoom
is the finding that step 3 (cancellation) must fix; the one-event numbers
show that even one view costs five to seven times the byte budget today.

What the numbers already say: the ten-section store's gene load reads
54 MB and holds 22 M molecules for an 11 M-molecule view, with 887 draw
calls; a slider step on a heavy scene is seconds, not 100 ms; the time to
first molecules is within budget on a single section (0.3 s) and three
times over on the atlas (2.3 s); frames on a software rasterizer are
seconds, which is why uploads paced per frame stall there.

Numbers the writer measured before this harness existed (three genes on
the full store: Slc17a7 589 M, Slc1a3 286 M, Gfap 284 M):

| View | Levels | Drawn | Read |
|---|---|---|---|
| 40 mm (all ten sections) | 5, 4, 4 | 2.8 M | 29 MB |
| 9 mm (one section) | 3, 2, 2 | 3.6 M | 58 MB |
| 1.5 mm | 0, 0, 0 | 4.8 M | 55 MB |
| zoom further, or back out | | | nothing new |

At 12 MB/s those views settle in 2.4 s, 4.8 s and 4.6 s. The target is 1 s.

## After step 3 (2026-10-17, one store, cancellation, byte budget)

Same fixtures, same server at 12 MB/s, same headless Chromium
(`perf/results/latest.json` at the step's commit). The Auto detail is now
8 M molecules and 12 MB a view (was 15 M, unbounded in bytes).

| Metric | spatialdata-tiny | spatialdata-medium | before |
|---|---|---|---|
| openToCells | 2.8 s | 4.8 s | 3.6 s |
| firstMolecules (default genes) | 3.1 s | 5.4 s | 4.1 s |
| geneToSettled (18 M-molecule gene) | — | **2.7 s** | 14.9 s |
| geneBytesKB / geneRequests | — | 6.8 MB / 25 | 24 MB / 148 |
| viewSettled (one wheel event) | 0 | **13.5 s** | 140 s |
| bytesPerViewKB / requestsPerView | 0 | 16.6 MB / 95 | 59 MB / 364 |
| viewSettledStepped (four wheel events) | 0.3 s | **22.6 s** | never |
| sliderToFrame (light) | 71 ms | 104 ms | 100 ms |
| sliderToFrameHeavy | — | **3.0 s** | 32 s |
| avgRenderMs | 0.3 ms | 1.1 ms | 2.6 ms |
| drawCalls | 10 | 79 | 282 |

What changed the numbers: a view reads the level its byte share allows,
so the heavy gene costs 6.8 MB instead of 24; a replan cancels every read
the new view does not want, so a zoom in several wheel events settles at
all; the coarsest level of a fresh gene is read first. What did not: the
one-event zoom still reads 16.6 MB against a 12 MB target, because a
piece covers a power-of-two span of tiles with a row of padding, so the
molecules read exceed those in view (step 4 narrows pieces); times on
SwiftShader include frames of hundreds of milliseconds, so the
13.5 s to settle is mostly drawing 8 M molecules on a software rasterizer,
not the 1.4 s the bytes take at 12 MB/s. `openToCells` on the medium store
is latency-bound: its table is 1,092 ranged reads of 36 bytes (the
generator's sharding), and the store's queue of 32 in flight is what
keeps it near the unbounded number.

## After step 4 (2026-10-17, rendering)

Same fixtures and server. Molecule clouds now carry positions only, a
cell tick reuses its buffers, hover picks through a grid.

| Metric | spatialdata-tiny | spatialdata-medium | after step 3 |
|---|---|---|---|
| sliderToFrame (light) | 64 ms | 95 ms | 104 ms |
| sliderToFrameHeavy | — | 2.6 s | 3.0 s |
| hoverFrameMs | 50 ms | 49 ms | 50 ms |
| avgRenderMs | 0.2 ms | 1.1 ms | 1.1 ms |
| jsHeapMB | 58 | 61 | 82 |
| geneToSettled / viewSettled / stepped | — | 2.9 s / 16.1 s / 22.5 s | 2.7 / 13.5 / 22.6 s |
| drawCalls | 10 | 79 | 79 |

The heavy slider step is no longer a repaint of every molecule's colour
(nothing is uploaded now: `jsHeapMB` fell by a quarter) but the time
SwiftShader takes to draw 8 M molecules once; on a GPU that frame is the
16 ms budget's concern, measured by `scripts/perf/real-store.ts` on a
machine with one. The medium store's `hoverFrameMs` is the same because
the pointer sweep's worst frame is a draw, not the pick.

## After SDDB step 2 (2026-10-18, coordinate systems)

Same fixtures and server. Images are now placed through the store's
coordinate systems (`lib/spatialdata/coordinates.ts`); read once from the
consolidated metadata, so nothing is fetched for it.

| Metric | spatialdata-tiny | spatialdata-medium | after step 4 |
|---|---|---|---|
| openToCells | — | 4.7 s | 4.8 s |
| sliderToFrame (light) | — | 109 ms | 95 ms |
| sliderToFrameHeavy | — | 2.7 s | 2.6 s |
| hoverFrameMs | — | 49 ms | 49 ms |
| geneToSettled / viewSettled | — | 2.3 s / 14.4 s | 2.9 s / 16.1 s |
| drawCalls | 10 | 79 | 79 |

Unchanged within the run-to-run noise. The real Visium stores were opened
through Lamin's proxy for the placement only (the cerebellum
`fJSOO03VzCTYtD1o`, the ovary `haOzFbJ7Bc9x7NW6`: the image lies under the
spots); their load times are not measured yet (step 4's nightly).

## How the numbers are read

- Layer hooks (`window.__merfish.layers[kind]`, `lib/utils/test-hooks.ts`):
  `firstContentAt`, `changedAt` (the view asked for something new),
  `settledAt` (everything wanted is drawn), `drawn`, `wanted`, `resident`.
- Frames (`__merfish.frames`): intervals between ticks of the scene loop
  (every animation frame, drawn or not, so a still on-demand scene is not a
  slow one); `over16ms` counts intervals above 20 ms (one dropped 60 Hz
  frame); `drawn` and `renderMs` cover only the frames actually drawn.
- Bytes: the Resource Timing API when the server sends
  `Timing-Allow-Origin` (the fixture server does); the proxy's log for the
  real bucket.
- Memory: `performance.memory` (Chrome) and three.js `renderer.info`.
- **Headless Chromium draws with SwiftShader (a software rasterizer)**, so
  its frame intervals measure CPU rasterization of millions of points, not a
  GPU: `frameOver16Pct` and `worstFrameMs` from the suite are only
  comparable with themselves, and anything paced per frame (uploads) is
  slower there. `avgRenderMs` (CPU inside `render`) is the gate on CI; frame
  intervals are read from the real-store script on a machine with a GPU.
- The medium fixture is 30 M molecules (`make-medium-store.py`): its
  heaviest gene, 18 M, exceeds the Auto budget, so opening it whole draws a
  coarser level and zooming in reads the finer one.
