# Reading the stitcher's images, and making one image of five channels

Two jobs that look like one. **A:** merfisheyes refuses the CoQ1s stitcher's
`.ome.zarr` output, which is otherwise fine. **B:** that output is five
one-channel images, and a user wants one five-channel image. **C** follows
from A: stop the stitcher writing it that way again.

Nothing here is written yet. Every number below is measured from
`/media/stephenSSD/20260710T132002_B6.3_G93A126_135_126/Coq1s_Output/W1/A1`.

## 1. What the data is

Five folders, `W1A1_C1.ome.zarr` … `W1A1_C5.ome.zarr`, identical in shape:

| | |
|---|---|
| Container | zarr v3 (`zarr.json`, no `.zattrs`) |
| NGFF metadata | `attributes.multiscales`, **not** `attributes.ome.multiscales`; `multiscales[0].version` is `"0.4"` |
| Axes | `c, y, x` — one channel each, `omero.channels[0].label` = `"C1"` … `"C5"` |
| Level 0 | `[1, 73519, 154331]`, `uint16`, zstd |
| Levels | 10, each half the last, down to `[1, 143, 301]` |
| Scale | 0.65 µm/px at level 0 |
| Chunks | `[1, 2048, 2048]` at every level |
| Chunk keys | `<level>/c/<ci>/<yi>/<xi>`, `/` separator; `ci` is always `0` |
| Size | **30 GB per channel, ~150 GB for A1** |
| Level-0 chunk bytes | median **2.98 MB**, max 5.69 MB compressed (8 MB raw) |
| Coarsest level | one tile, 44 KB |

## 2. Why it is refused, and where it came from

OME-NGFF 0.5's own schema
(`https://ngff.openmicroscopy.org/0.5/schemas/image.schema`) has **`ome` as its
only required top-level property**, with `multiscales` and `version` required
inside it. 0.4 is a zarr **v2** layout. This data is neither: a v3 container
carrying 0.4-style attributes.

Three readers in this repo require `attributes.ome.multiscales` and so reject
it — `lib/upload/uploadImageToS3.ts:25` (`isOmeZarrImage`, the gate for
drag-and-drop and for Upload & Save), `lib/spatialdata/image.ts:203`
(`readMultiscale`), `lib/spatialdata/resolveUrl.ts:32` (`zarrKind`).

**Four others already accept it**: `lib/spatialdata/coordinates.ts:389`,
`scripts/validate_store.py:333`, `scripts/sddb_census.py:199`,
`scripts/prepare_folder.py:296`, each with its own inline `(attrs.ome ?? attrs)`.
So a store can load halfway. *That split is the actual defect*; the stitcher
just exposed it.

**Where the files came from** (reproduced, `scratchpad/repro_hybrid.py`):

- `ome-zarr 0.21.0` + `zarr 3.4.1` → `CurrentFormat()` is `FormatV05`, and the
  output is **valid 0.5**, `attributes.ome` and all.
- Forcing `FormatV04` under zarr 3 is **refused by ome-zarr-py itself**:
  `ValueError: Group is zarr_format: 3 but OME-Zarr 0.4 is 2`.

So this data was written by an **older ome-zarr-py**, whose `CurrentFormat()`
was 0.4 and which had no such guard, running against zarr-python 3. Today's
libraries cannot produce it.

**But the reference reader reads it.** `ome_zarr.reader.Reader` opens
`W1A1_C1.ome.zarr` and returns all ten levels, warning only
`version mismatch: detected: FormatV04, requested: FormatV05`. So upstream
refuses to *write* this shape yet still *reads* it — which is exactly the
position this plan takes for merfisheyes, and the precedent for the MAY row
below. Note `CoQ1s/numpy_zarr.py`'s
`require_version("zarr", max_major=2)` is dated 2026-05-11 and the data
2026-07-24, so the run did not use this checkout or this environment — worth
confirming on the imager before changing anything there.

## 3. Part A — accept it

**Contract first.** `docs/specs/image-profile.md` currently makes 0.5 under
`attributes.ome` a MUST. Replace that row with two:

- 0.5 under `attributes.ome` — MUST for anything we write.
- 0.4-style `attributes.multiscales` in a v3 container — **read** MAY,
  transitional, never written. Name the producer and the upstream error above,
  so the next reader knows it is deliberate and not a spec misreading.

Change the row, `docs/specs/schemas/`, and `scripts/validate_store.py` in one
commit, as this repo requires of a format change.

**Then one helper, not three patches.**

```ts
/** NGFF metadata, wherever the writer put it (image-profile.md). */
export const omeAttrs = (attrs: Record<string, unknown>) =>
  (attrs.ome ?? attrs) as OmeAttrs;
```

One expression covers `multiscales` **and** `omero`, because in the 0.4 layout
they are siblings. Call sites: `image.ts:203`, `uploadImageToS3.ts:25`,
`resolveUrl.ts:32`, and `coordinates.ts:389` switched off its inline copy so
there is one implementation.

**Tests — the inconsistency is what must not come back.**

- `tests/unit/ome-layout.test.ts`: a table over *every* reader entry point ×
  both layouts, in the style of `split-sync-capabilities.test.ts`'s pair
  matrix, so adding a reader that forgets the helper fails here.
- Fixture `tests/data/tiny/zarr/image_v04_in_v3.zarr`, built by
  `scripts/testdata/make-tiny-ome-04.py` — hand-written metadata, since
  ome-zarr-py now refuses to emit it.
- `@quick` browser test through `e2e/helpers/fixture-server.ts`:
  `/image-viewer/from-s3` on that fixture draws.
- Python side (`prepare_folder.py`, 5 spots) **only** if these images are to go
  through `prepare.py`; drag-and-drop and save-as-is do not touch it.

## 4. Part B — one image of five channels

The chunk keys decide this. A combined array of shape `[5, y, x]` with chunks
`[1, 2048, 2048]` addresses channel *k* at `<level>/c/<k>/<yi>/<xi>`; each
source channel already holds exactly those bytes at `<level>/c/0/<yi>/<xi>`.

**B1 — relink. Verified on this data** (`scratchpad/prove_combine.py`): a
two-level, five-channel store built from the real files as ten symlinks and no
copied pixels reads back `(5, 287, 602)` and `(5, 143, 301)`, every channel
**pixel-identical** to its source array. At full depth that is fifty symlinks
and seconds of work. Write the root as proper 0.5 while you are at it, so the
combined store is valid even though its inputs are not.

Two caveats: it is local-only (an upload follows the links and copies the
150 GB anyway), and whether Chrome's directory picker traverses symlinks for
drag-and-drop is **not yet tested** — check that before relying on it.

**B2 — rewrite.** One real `[5, y, x]` store. 150 GB read and written. Needed
for S3, Explore, or anything shared.

**The chunking decision, which is the real content of B2.** Measured:
`[1, 2048, 2048]` is 8 MB raw, median 2.98 MB compressed, against
`image-profile.md`'s SHOULD of ≤1024² so a tile is ≤2 MB. With five channels
visible, `ImageLayer`'s 12-tile view costs **5 × 12 × ~3 MB ≈ 180 MB** at level
0. Rewriting at `[1, 1024, 1024]` quarters that to **~45 MB**. Keep zstd, keep
10 levels. Measure both with the `e2e/perf` harness against the tile budget
rather than assuming the arithmetic holds.

Also: the channel labels are literally `C1`…`C5`. Ask what the stains are and
write real labels into `omero`, or the channel list in the viewer says nothing.

## 5. Part C — stop writing it

`ome-zarr 0.21 + zarr 3` already writes valid 0.5, so the fix is to **drop the
`zarr<3` pin**, not to hand-write metadata.

What that pin is holding back, in `CoQ1s/numpy_zarr.py`: `zarr.GZip` / `LZ4` /
`Zstd` / `Blosc` (≈ lines 223–229), `zarr.group(store=…)`, and
`parse_url(...).store` — all zarr-v2 APIs that moved in v3. Port those, set
`zarr>=3` and `ome-zarr>=0.21` in `requirements.txt` and `environment.yml`, and
turn `require_version` around: require ome-zarr ≥ 0.21, so the broken pairing
cannot recur silently. Confirm with a tiny write that `attributes.ome` and
`version: "0.5"` come out.

## 6. Order

1. **A** — the viewer accepts today's 150 GB as it stands. Small.
2. **C** — new runs are valid. Medium, gated on porting the compressor calls.
3. **B** — decide B1 vs B2 after A lands and one channel has been looked at.

## 7. What says it worked, other than us

- **`napari-ome-zarr` / `ome-zarr-py`** opening the combined output — already in
  the stitcher's `requirements.txt`, so it is their own reader, not ours.
  (It already reads the un-combined hybrid; §2.)
- The **ome-ngff validator** and `python scripts/validate_store.py`.
- A **round-trip**: read the combined store back and compare sampled tiles
  pixel-for-pixel against the five source arrays.
- The **perf harness on the real A1**, not the tiny fixture: tiles and bytes per
  view, held under the profile's budget.
- `image-profile.md` alone being enough for someone else to write an image this
  viewer accepts.

## 8. Open questions

- What are the five stains? (channel labels)
- Do these need to reach S3 / Explore, or stay local previews? That decides B1
  vs B2.
- Which `ome-zarr` version is on the imager that produced this?
