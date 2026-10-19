# Shapes profile

What a shapes element must be for `lib/spatialdata/shapes.ts` and
`ShapesLayer` to draw cell polygons in each cell's colour.

| Requirement | Level |
|---|---|
| `shapes/<name>/shapes.parquet`: **one** GeoParquet file. A directory of parquet parts (as dask-geopandas writes) is not read: the viewer asks for that one key and gets a 403 | MUST |
| `<name>` one of `cell_boundaries`, `cell_polygons`, `cells` (the first found, in that order) | MUST |
| GeoParquet `geo` metadata naming the `primary_column`; the column holds WKB polygons (or GeoJSON-like objects); the outer ring is drawn | MUST |
| An index column (the `pandas` metadata's first `index_columns` entry, else the first non-geometry column) holding each polygon's instance id, equal to the table's `instance_key` column value for that cell | MUST for the polygons to take their cell's colour; unmatched polygons are drawn grey |
| `coordinateTransformations` of the element: `identity`, or a `scale` (and `translation`) to the shared system; the image is placed by undoing this transform (`image-profile.md`) | MUST |
| Row groups ordered by position (a tile-major order), so a region can be read by row group | SHOULD (today the whole file is read on the first zoom-in: 31 MB on atera); the planned ranged reads need it |
| At most 60,000 polygons drawn per view; a 100 µm grid indexes them | behaviour |

Reading: the whole file once, in a worker (`lib/workers/shapes.worker.ts`,
hyparquet), decoded to the outer rings and ids; polygons are built for the
view when it is narrower than 1 mm, filled or outlined, under the molecules
and the cell dots (`docs/specs/RENDERING.md`).

Known stores: atera has one file (31 MB) and matches its table; the skin
store has a directory of parts and is not drawn. `scripts/validate_store.py`
reports both.

## Spots and bins: the cells' footprint

When the table's region (the shapes element it annotates) is a Visium
spot or Visium HD bin element, the shapes are not read: each cell is drawn
at its real size and shape by the cell cloud itself
(`lib/spatialdata/footprint.ts`, `docs/specs/RENDERING.md`), at any zoom,
in its own colour, picked and selected as a cell.

| Requirement | Level |
|---|---|
| The region element's GeoParquet `geo` metadata gives one geometry type | MUST for a footprint |
| Circles: `Point` geometry and a `radius` column whose statistics (min, max) are in the footer and equal — one radius for every spot (Visium). Xenium's `cell_circles` vary and are drawn as dots | MUST for circles |
| Squares: `Polygon` geometry, an element name containing `square`, the bins tiling a grid; the side is taken as the grid's pitch, the median nearest-neighbour distance of the table's coordinates (a Visium HD grid is rotated by a fraction of a degree in the image's pixels, so its columns do not share an x) | MUST for squares |
| The region in the cells' frame (`store-layout.md` §8) | MUST |

Read: one suffix range of the region's `shapes.parquet` (its footer); the
2 µm bins of a Visium HD store are 348 MB and are never fetched.
