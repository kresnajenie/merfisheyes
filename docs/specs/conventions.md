# Store conventions

The viewer and `scripts/validate_store.py` meet three kinds of rule. Each
has one home:

| Kind | Example | Where it lives |
|---|---|---|
| A rule of the format | an element's transforms are alternatives, one per system; a placement may rotate | code, defined by the profiles in this folder (`store-layout.md`, `table-profile.md`, `image-profile.md`, …) and tested |
| A name a producer uses | SpatialData-DB's `<region>_for_vitessce`, Visium's `array_row`, Visium HD's `square_` bins, MALDI's `ion_image_pixel_x` | [`conventions.json`](conventions.json), schema [`schemas/conventions.schema.json`](schemas/conventions.schema.json) |
| A choice for one dataset | its default column, flips, rotation | the dataset's `viewerConfig` in the database, set by its owner |

`conventions.json` is read by the viewer through
`lib/spatialdata/conventions.ts` and by the validator at import, so the two
cannot disagree about a name.

## What is in it

| Key | Read by | What it decides |
|---|---|---|
| `tables.curatedDefaultTable`, `tables.fallbackTableName` | `defaultTableName` (`lib/spatialdata/store.ts`), `default_table` (validator) | which table a store's cells are read from |
| `tables.featureIndexColumns` | the SpatialData, h5ad and h5ad-zarr adapters | `var` columns tried for feature names after the one `var` declares |
| `tables.defaultColumnPriority` | `selectBestClusterColumn*` (`lib/utils/dataset-utils.ts`), `de-stats.ts` | obs columns chosen first as the default colouring |
| `tables.notDefaultColumns` | the same | positional columns never chosen as the default colouring |
| `shapes.cellOutlineElements` | `cellShapesName` (`lib/spatialdata/shapes.ts`), validator | shapes elements the Shapes layer draws as outlines |
| `shapes.pointSiblingSuffix` | `centroids.ts`, `footprint.ts`, validator | the Point element beside a region with each cell's centre and radius |
| `shapes.squareBinPattern` | `footprint.ts` | a polygon region drawn as square bins of the grid's pitch |
| `catalog.bucketRegions`, `catalog.species`, `catalog.platforms` | `lib/sddb/catalog-row.ts` | a SpatialData-DB catalog row's URL and spellings on Explore |

## Adding one

**A new name for an existing rule** (a producer calls its positional column
`grid_x`, its outlines `segmentation`, a new technology needs a spelling):
add it to `conventions.json`. The schema check (`tests/python/test_specs.py`)
and `tests/unit/conventions.test.ts` cover the file; add a case to the test
of the reader that uses it when the name changes what a real store does.

**A new kind of rule** (something no key above decides): add a key to the
schema and the file, read it in `lib/spatialdata/conventions.ts` and, when
the validator checks it, in `validate_store.py`; describe it in the table
above. If it is part of the format rather than one producer's habit, it
belongs in a profile and in code instead.

`tests/unit/conventions.test.ts` fails when a reader in `lib/` or
`components/` writes one of the file's distinctive names as a string
literal, so a name added here is not also hard-coded elsewhere.
