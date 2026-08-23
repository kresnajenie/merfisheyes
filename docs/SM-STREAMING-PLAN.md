# Single-molecule streaming (plan) — branch `feat/sm-streaming`

Problem: a single gene can be 50M molecules (~600 MB of float32). The viewer
downloads the whole `genes/{gene}.bin.gz`, gunzips it in one go, then builds the
point cloud — minutes of spinner. The file format already allows progressive
consumption (flat `[x,y,z,…]` float32 in one gzip stream), so stream it.

## Decisions (user-approved 2026-08-23)
- **Threshold: per gene file.** Stream when `manifest.genes.molecule_counts[gene].assigned`
  (or `.unassigned` for the `_uuuuuuuuuu` file) ≥ `SM_STREAM_THRESHOLD = 2_000_000`
  (in `lib/config/visualization.config.ts`); smaller files keep the one-shot load.
  Counts are always present after `feat/pipeline-parity`; if absent (legacy
  manifest) → stream when the response `Content-Length` ≥ ~8 MB (gzipped), else one-shot.
- **Progressive order in the files.** Both processors write each gene's molecules in
  **bit-reversal index order** so the first N bytes of any stream are an evenly spread
  subsample of the whole slide (the picture densifies instead of sweeping across).
  Deterministic, no RNG → byte-identical JS/Python; old files still stream in file order.
- Branch stacked on `feat/pipeline-parity` (writers change); merge after parity.

## 1. Writers: progressive molecule order
Permutation for a gene with `n` molecules: sort indices `0..n-1` by
`bitrev(i, ceil(log2 n))`, skipping values ≥ n (the standard "bit-reversal permutation
of the next power of two, dropping out-of-range"). Implemented once per language:
- `scripts/process_single_molecule.py` → `progressive_order(n) -> np.ndarray[uint32]`
  (vectorised: bit-reverse with numpy on `np.arange(2**k)`, filter `< n`), applied to
  each gene's index array before slicing coords (assigned and unassigned separately).
- `lib/SingleMoleculeDataset.ts` (browser parse) → `progressiveOrder(n): Uint32Array`,
  applied when filling each gene's Float32Array (pass 2 writes into the permuted slot;
  or permute after). `lib/utils/SingleMoleculeProcessor.ts` unchanged (writes what the
  dataset holds).
- Unit test: `tests/unit/progressive-order.test.ts` — permutation is a bijection, first
  k entries of a spatially sorted input cover the range evenly; Python twin via a
  golden file decoded in TS (same pattern as `tests/fixtures/obs-binary/`).
- No manifest change needed; optionally `processing.molecule_order: "bitrev-v1"` so
  the viewer can say "progressive preview" vs "loading in file order".

## 2. Viewer: streamed gene loading
- `lib/SingleMoleculeDataset.ts` (fromS3 + fromCustomS3 overrides): add
  `streamCoordinatesByGene(gene, { onChunk(coords: Float32Array, loaded, total), signal })`
  and the unassigned twin:
  1. resolve URL as today (`/api/single-molecule/{id}/gene/{gene}`; custom-S3 direct URL);
  2. `fetch(url, { signal })` → `response.body.pipeThrough(new DecompressionStream("gzip"))`
     → reader loop; keep a carry buffer so only complete 12-byte triples are emitted;
  3. emit in batches of ≥ 250k molecules or every ~120 ms (whichever first) to keep
     `requestAnimationFrame` smooth; apply `denormFactor` per batch (legacy);
  4. on completion, cache the full Float32Array in `geneIndex` (so re-selection is the
     instant path), mark `moleculeCounts`, resolve. On abort, discard.
  - Keep `getCoordinatesByGene` as is (one-shot, cached) for small files and local datasets.
  - Memory: stream into a preallocated `Float32Array(count*3)` when the count is known
    (manifest), else grow by doubling and `slice` at the end.
- `components/single-molecule-three-scene.tsx`:
  - new `createStreamingPointCloud(count, geneViz, shape)`: `BufferGeometry` with a
    preallocated position attribute (`count*3`), `setDrawRange(0, 0)`,
    `attribute.setUsage(THREE.DynamicDrawUsage)`; per batch: copy into the attribute,
    `addUpdateRange(offset, length)` + `needsUpdate = true`, `setDrawRange(0, loaded)`,
    recompute bounding sphere only every N batches (or from the manifest's coordinate
    extent if we add it; otherwise from the first batch — bit-reversal order makes the
    first batch representative of the full extent);
  - decide per gene: `count >= SM_STREAM_THRESHOLD` → streaming path, else current path;
  - cancellation: `AbortController` per gene load; deselecting / changing view mode /
    unmount aborts; already-rendered partial cloud is disposed;
  - toast: `Loading Gad1 — 12.3M / 50.0M (25 %)` updated per batch, then dismissed;
  - unassigned (`_uuuuuuuuuu`) cloud: same machinery, same threshold.
- Colour/size/visibility updates during a stream must not restart it (they mutate
  the material only — already the case).
- SC viewer overlay (`components/three-scene.tsx` loads linked SM genes via the same
  dataset API): adopt `streamCoordinatesByGene` there too, or leave one-shot for v1 (decide
  during implementation; overlay genes are usually smaller).

## 3. API / infra
- No new endpoint: presigned GET URLs work with streaming fetch. Confirm S3 objects
  are served with `Content-Encoding` unset (raw .gz bytes) — they are today.
- CORS: streaming reads need the presigned-URL host to allow the origin (already
  required for the current fetch).

## 4. Tests / verification
- Unit: progressive-order bijection + coverage; a `streamFloat32Triples(readable)` helper
  tested with a chunked fake stream that splits mid-triple.
- Manual: a ≥ 2M-molecule gene on dev (e.g. the 46M-molecule spots_set4 file once re-downloaded
  intact) — first points within ~1 s, smooth densification, deselect mid-stream aborts,
  re-select instant from cache; memory after completion equals today's.
- Parity harness (docs/PARITY-TESTING.md): add the order to the acceptance list — gene
  files from both pipelines must be byte-identical including order.

## 5. Rollout
1. Viewer streaming (works on existing files, file order).
2. Writer progressive order (both pipelines) + worker image rebuild.
3. Optional later: LOD/tiles (spatial streaming) — separate design.
