# Points in Zarr: a multiscale, doubly ordered points element (draft spec)

Status: **draft 0.1, proposed by merfisheyes, 2026-10-16.** A proposal for
storing a SpatialData points element (molecules, transcripts, any point
cloud with a categorical feature) as Zarr arrays instead of Parquet, so
it can be viewed at any zoom and queried by feature or by region over
HTTP. Background, measurements and the alternatives considered are in
[molecule-pyramid.md](molecule-pyramid.md); this file is the format itself.

The key words MUST, SHOULD and MAY are used as in RFC 2119.

## 1. Overview

A points element is a Zarr v3 group under `points/<name>/` in a SpatialData
store. It holds:

- **levels**: level `0` is every point; each further level is a uniform
  random sample of the one below it (a quarter of it, by default);
- in each level, the same points **in two row orders**:
  - `tile_major`, grouped by spatial tile — for "everything in this region";
  - `feature_major`, grouped by feature (e.g. gene) — for "this feature
    everywhere";
- in each level, an **index** of row numbers, small enough to read before
  any points, that says where each tile and each feature's rows are and how
  many points each feature has in each tile.

Each order is a dataframe of 1-D columns of equal length.

```
points/<name>/                      zarr.json: element metadata (§2)
├── 0/                              level 0: every point
│   ├── tile_major/                 dataframe (§4), every column
│   ├── feature_major/              dataframe (§4): coordinates, feature, id
│   └── index/                      §5
├── 1/                              a uniform sample of 0
└── …
```

It builds on three Zarr conventions:

| Convention | UUID | Used for |
|---|---|---|
| [multiscales](https://github.com/zarr-conventions/multiscales) v0.1 | `d35379db-88df-4056-af3a-620245f8e347` | declaring the levels |
| [dataframe](https://github.com/ilan-gold/zarr-convention-dataframe) (`df:`) | `1f1b36d2-3185-6da8-91f6-9a0823af2d33` | each row order |
| [enum](https://github.com/ilan-gold/zarr-convention-enum) (`enum:`) | `c906c423-56ed-413c-9943-7b2ff52d18f2` | categorical columns (feature, sample, cell) |

Properties this spec adds are namespaced `points:`.

## 2. Element metadata

The element group's attributes MUST contain:

| Attribute | Value |
|---|---|
| `encoding-type` | `"ngff:points"` |
| `axes` | `["x", "y"]` or `["x", "y", "z"]` |
| `coordinateTransformations` | as for any SpatialData element: intrinsic → named coordinate systems |
| `spatialdata_attrs` | `{"feature_key": <feature column>, "version": "0.3-zarr"}`, plus `instance_key` if a column assigns points to instances |
| `zarr_conventions` | entries for `multiscales` and for this spec (§9) |
| `multiscales` | the levels (§3) |
| `points:version` | `"0.1"` |
| `points:units` | unit of the intrinsic coordinates, e.g. `"micrometer"` |
| `points:n_points` | number of points in level 0 |

and SHOULD contain `points:extent`: `[[xmin, ymin(, zmin)], [xmax, ymax(, zmax)]]`
of level 0.

## 3. Levels

`multiscales.layout` lists the levels from level `0` (every point) to the
coarsest. Each entry:

| Key | Value |
|---|---|
| `asset` | the level group's name: `"0"`, `"1"`, … |
| `derived_from` | the level below (absent on level `0`) |
| `transform` | `{"scale": [1, 1(, 1)]}` — points keep their coordinates across levels |
| `resampling_method` | `"random"` |
| `points:n_points` | points in this level |
| `points:grid` | the level's tile grid (§3.1) |

Rules:

- Level `0` MUST contain every point of the element exactly once.
- Every point of level *k* + 1 MUST be in level *k* (levels are nested).
- Each level SHOULD be a uniform random sample of the level below,
  independent of feature and position. The reference writer draws one
  uniform random number `r` per point and puts it in every level *k* with
  `r < 4^-k`.
- Levels SHOULD be added until a level holds about a million points.

### 3.1 Tile grid

```json
"points:grid": { "origin": [x0, y0], "tile_size": 100.0, "shape": [ny, nx] }
```

A point at (x, y) is in tile `t = ty * nx + tx`, where
`tx = floor((x - x0) / tile_size)` and `ty = floor((y - y0) / tile_size)`,
both clamped to the grid. Tiles are 2-D even for 3-D points. Each level
has its own grid; the reference writer doubles `tile_size` per level, so
tiles hold similar numbers of points at every level.

## 4. The two row orders

Both are groups following the dataframe convention: a `columns` array of
names and members `column_0` … `column_{n-1}`, all of length
`points:n_points` of the level. Column names:

| Name | Type | In `tile_major` | In `feature_major` |
|---|---|---|---|
| `x`, `y` (`z`) | float32 or float64, intrinsic coordinates | MUST | MUST |
| the feature column (`spatialdata_attrs.feature_key`) | enum, unsigned integer codes | MUST | MUST |
| `id` | uint64 | — (implicit: the row number) | MUST |
| any other column | any 1-D type or enum | MAY | MAY, if listed in `points:copied_columns` |

**Row order.**

- `tile_major` MUST be sorted by (tile, feature code). Within a (tile,
  feature) run, level 0 SHOULD be sorted by `r`, so any prefix of a run is
  itself a uniform sample; levels above 0 SHOULD be sorted by `id`, so the
  ids rise and delta-code to almost nothing (those levels are samples
  already).
- `feature_major` MUST be sorted by (feature code, tile), then in the same
  order as `tile_major` within each (tile, feature) run.

**Ids.** A point's id is its row number in **level 0** `tile_major`. The `id`
column of every `feature_major`, at every level, holds those level-0 row
numbers, so a point drawn from any level can be looked up in the full table.
A `tile_major` of level *k* > 0 MAY also carry an `id` column.

**Group attributes** of each order: `df:shape` (dataframe convention),
`points:order`: `"tile_major"` or `"feature_major"`, and on `feature_major`
`points:copied_columns`: the names of columns beyond coordinates, feature
and id.

The feature column's `enum_values` MUST be identical in both orders and at
every level, so a code means the same feature everywhere.

## 5. Index

Each level has an `index/` group of 1-D unsigned integer arrays, all in
**row numbers** of that level (never in chunks, so chunking can change
freely):

| Array | Type | Shape | Meaning |
|---|---|---|---|
| `tile_rows` | uint64 | `nx·ny + 1` | tile *t* is `tile_major` rows `[tile_rows[t], tile_rows[t+1])` |
| `feature_rows` | uint64 | `n_features + 1` | feature *f* is `feature_major` rows `[feature_rows[f], feature_rows[f+1])` |
| `feature_tile_ptr` | uint64 | `n_features + 1` | feature *f*'s entries are `[feature_tile_ptr[f], feature_tile_ptr[f+1])` of the next two |
| `feature_tile` | uint32 | R | the tiles holding feature *f*, ascending |
| `feature_tile_count` | uint32 | R | how many of *f*'s points each of those tiles holds |

where R is the number of non-empty (feature, tile) pairs.

Derived positions:

- Feature *f* in tile *t*, in `feature_major`: starts at `feature_rows[f]`
  plus the counts of *f*'s tiles before *t*.
- In `tile_major`: starts at `tile_rows[t]` plus the counts, in tile *t*, of
  the features with a lower code.

## 6. Storage

Not normative; what the reference writer does:

- Point columns: chunks of 65,536 rows, shards of 64 chunks.
- Coordinates: `bytes` → `blosc` (zstd, byte shuffle).
- `id`: `numcodecs.delta` → `blosc` (ids increase within a (feature, tile) run).
- Codes: `blosc` (zstd, byte shuffle).
- Index arrays: chunks of 65,536 entries.
- Consolidated metadata at the store root.

All of these are decodable by zarr-python ≥ 3 and zarrita ≥ 0.7.

## 7. Reading

**A viewer**, for a set of selected features and a view rectangle:

1. Read the element metadata and, for each level, `feature_rows` and
   `feature_tile_ptr` (`n_features + 1` entries each).
2. For each selected feature, read its slices of `feature_tile` and
   `feature_tile_count` at each level.
3. From the counts, pick per feature the finest level whose points in view
   fit a point budget.
4. Read `x`, `y` (`z`) of `feature_major` for that feature's tiles in view:
   contiguous row ranges, one per run of adjacent tiles.

**Everything in a region**: `tile_rows` of the tiles overlapping it →
`tile_major` row ranges → exact filter on coordinates.

**One feature, all points**: `feature_rows` at level 0 → one row range of
`feature_major`.

**As a dataframe**: level 0 `tile_major` is the whole element; its columns
read straight into pandas, dask or polars.

## 8. Adding a column

Write it as `column_{n}` of level 0 `tile_major`, aligned to its rows, and
append its name to `columns` and to `df:shape`. To make it available at other
levels or in `feature_major`, gather it by `id`.

## 9. Convention entry

Until this spec has a home, the element's `zarr_conventions` lists it as:

```json
{
  "name": "points:",
  "description": "Multiscale points: levels, tile- and feature-major orders, row index",
  "spec_url": "https://github.com/kresnajenie/merfisheyes/blob/develop/docs/points-zarr-spec.md"
}
```

## 10. Compatibility

- spatialdata 0.8 reads points only from `points.parquet` and fails on this
  element. Until it understands this format, open such a store with
  `sd.read_zarr(path, selection=("tables", "images", "labels", "shapes"))` and
  read the points with the reference reader.
- The element doesn't use `points.parquet`, so writing one beside it (for
  today's readers) is possible and doesn't conflict.

## 11. Reference implementation

- Writer, reader and validator (Python): `scripts/points_zarr.py`
  (`python points_zarr.py validate <element>`).
- Viewer: merfisheyes, `lib/spatialdata/points.ts` (reader) and
  `lib/spatialdata/PointsLayer.ts` (drawing a view).
- Tiny example: `tests/data/tiny/zarr/points_v1.zarr`, made by
  `scripts/testdata/make-tiny-points.py`.

## 12. Proposals for 0.2

Kept apart from 0.1, which every store on the bucket is written in. A
reader of 0.1 must ignore what it does not know; a writer of 0.2 sets
`points:version` to `"0.2"`. The JSON Schema for the metadata of both is
[specs/schemas/points.schema.json](specs/schemas/points.schema.json).

- **`points:bytes_per_row` per level** (`multiscales.layout[k]`): the
  compressed bytes one `feature_major` row costs at that level, measured by
  the writer (2.1–3.2 B on the stores measured in
  [specs/AUDIT-2026-10-17.md](specs/AUDIT-2026-10-17.md)). A viewer with a
  byte budget per view ([specs/BUDGETS.md](specs/BUDGETS.md)) can then pick
  a level before reading anything; today it learns the figure from its own
  reads and assumes 3 B before the first.
- **Tile-relative `uint16` coordinates**: `x`, `y` as `uint16` offsets
  within the tile, in units of `tile_size / 65536`, with the tile known
  from the index. Halves the bytes on the wire against `float32` and
  compresses better; the viewer adds the tile origin on upload. Level 0 of
  a 100 µm tile resolves 1.5 nm.
- **`id` optional at level 0**: level 0 is every point in `tile_major`
  order, so its row number is its identity; `id` there costs 0.1–0.7 B a
  row and is only needed to join another column written later.
- **Index chunking**: `feature_tile` and `feature_tile_count` chunked by
  feature (`feature_tile_ptr` gives the boundaries), so one gene's slice is
  one or two reads on a store with 23.9 M pairs.
