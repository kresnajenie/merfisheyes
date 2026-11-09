# Image profile

What an image must be for `lib/spatialdata/image.ts` and `ImageLayer` to
draw it, inside a SpatialData store (`images/<name>`) or on its own (an
OME-Zarr group at a URL, an image dataset in the app's bucket).

| Requirement | Level |
|---|---|
| OME-NGFF 0.5 metadata under `attributes.ome`: `multiscales[0]` with `axes` and `datasets` | MUST for anything written here |
| Read also: 0.4-style metadata (`multiscales` and `omero` at the top of `attributes`) inside a Zarr v3 container, as ome-zarr-py wrote before it guarded against the pair and still reads, warning `version mismatch: detected: FormatV04`. Transitional, never written. Every reader goes through `omeAttrs` (`lib/spatialdata/ome.ts`); `tests/unit/ome-layout.test.ts` holds them to it, so a store cannot load halfway | MAY |
| Axes exactly `c, y, x`, in that order (a `z` or `t` axis is not read) | MUST |
| Data type `uint8` or `uint16` | MUST |
| Zarr v3 arrays; chunks at most 1024 × 1024 pixels per channel, so one tile is at most 2 MB | SHOULD (larger chunks are read whole: the skin store's 4096² chunks are 20 MB each) |
| Every level's `coordinateTransformations`: a `scale` and a `translation` (pixel → shared coordinates), or a `sequence` of them as the `spatialdata` library writes | MUST |
| Levels from finest (`datasets[0]`) to coarsest, each about half the previous | SHOULD; the coarsest is read whole for the first view, so it should be a few tiles |
| `omero.channels[].label` for the channel names | SHOULD |
| Any number of channels, one per chunk when there are many (a CODEX image: 84 channels, chunks `[1, 1024, 1024]`); the layer shows at most 8 at once (`MAX_VISIBLE_CHANNELS`: WebGL guarantees 16 texture units), the first by default, and turning a ninth on turns the earliest off | MAY |
| A colour photograph (an H&E stain): `uint8`, channels labelled `r, g, b` (or three unlabelled / index-labelled channels) — shown as red, green and blue, all on, over 0–255; fluorescence shows its first channel, windowed by the data | MAY |
| Channels may share a chunk (`[3, y, x]`, as every SpatialData-DB H&E image is written): one read fills every channel the chunk holds | MAY |
| Inside a store: an entry of `multiscales[0].coordinateTransformations` into the target system (`store-layout.md` §8: `global`, else the system the cells' region is identity into), followed by the inverse of the region's entry, so the image lands in the cells' intrinsic coordinates. Of several images, the finest that reaches the target is drawn | MUST when it is to line up with the cells |
| The transform into the target: `identity`, `scale`, `translation`, `affine` (its x, y part, which may rotate: a MALDI optical image sits at an angle to its spots) or a `sequence` of them. The levels stay in the image's own frame and the placement is one matrix on the layer, so a rotated image is drawn as is | MUST |
| Under a merfisheyes cell folder (`mapping.json` `image`): transforms must place the image directly in the cells' raw coordinates (microns); normalized cells draw no image | MUST |

How it is drawn: the finest level whose pixel is no smaller than 0.75 of
a screen pixel, within 12 tiles in view, with the coarsest level's tiles
always resident underneath (`ImageLayer`, `docs/specs/RENDERING.md`). Each
channel is one texture per tile; its window (min, max) defaults to the
99.9th percentile of the coarsest tile and is kept in the URL (`img=`).
