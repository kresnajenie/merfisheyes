# Image profile

What an image must be for `lib/spatialdata/image.ts` and `ImageLayer` to
draw it, inside a SpatialData store (`images/<name>`) or on its own (an
OME-Zarr group at a URL, an image dataset in the app's bucket).

| Requirement | Level |
|---|---|
| OME-NGFF 0.5 metadata under `attributes.ome`: `multiscales[0]` with `axes` and `datasets` | MUST |
| Axes exactly `c, y, x`, in that order (a `z` or `t` axis is not read) | MUST |
| Data type `uint8` or `uint16` | MUST |
| Zarr v3 arrays; chunks at most 1024 × 1024 pixels per channel, so one tile is at most 2 MB | SHOULD (larger chunks are read whole: the skin store's 4096² chunks are 20 MB each) |
| Every level's `coordinateTransformations`: a `scale` and a `translation` (pixel → shared coordinates), or a `sequence` of them as the `spatialdata` library writes | MUST |
| Levels from finest (`datasets[0]`) to coarsest, each about half the previous | SHOULD; the coarsest is read whole for the first view, so it should be a few tiles |
| `omero.channels[].label` for the channel names | SHOULD |
| Inside a store: the image's `multiscales[0].coordinateTransformations` to the shared system, undone by the cells' element transform (shapes first, else points) so the image lands in the cells' intrinsic coordinates | MUST when it is to line up with the cells |
| Under a merfisheyes cell folder (`mapping.json` `image`): transforms must place the image directly in the cells' raw coordinates (microns); normalized cells draw no image | MUST |

How it is drawn: the finest level whose pixel is no smaller than 0.75 of
a screen pixel, within 12 tiles in view, with the coarsest level's tiles
always resident underneath (`ImageLayer`, `docs/specs/RENDERING.md`). Each
channel is one texture per tile; its window (min, max) defaults to the
99.9th percentile of the coarsest tile and is kept in the URL (`img=`).
