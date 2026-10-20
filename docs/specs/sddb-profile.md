# SpatialData-DB stores against the MERFISHEYES profiles

Status: measured 2026-10-18. Decision record for `docs/plans/SDDB.md`
(step 1) and the written profile for the SpatialData-DB maintainers.

Every number below comes from the census: `scripts/sddb_census.py`, which
reads every latest SpatialData artifact of the public Lamin instance
`scverse/spatialdata-db` through the hub REST, each store's root
`zarr.json` and each parquet footer through the storage proxy, and runs
`scripts/validate_store.py --root-json` on each. Per-store records:
[`perf/sddb/census.json`](../../perf/sddb/census.json); counts:
[`perf/sddb/CENSUS.md`](../../perf/sddb/CENSUS.md). Re-run the script to
refresh both.

## The population

| | Stores |
|---|---|
| Latest SpatialData artifacts in the instance | 1,884 |
| … converted by SDDB (root has `sddb_provenance`) | 1,761 |
| … older stores at the bucket's top level (no provenance): 65 Visium, 44 Xenium; 12 earlier Visium HD versions; 2 `_csctest` | 123 |
| Listed on the website | 1,778 (17 more than the converted stores; not reconciled) |
| Total size | 3.1 TB |

By folder: Visium 935, spatial transcriptomics 552, Xenium 147, Visium HD
137, Atera 2 (plus the top-level stores). Tier, from the artifact
description: silver 1,481, bronze 267, none 136 (15 of them converted).
Licence, a Lamin label (not in the store): CC BY-NC-SA 4.0 1,090; CC BY 4.0
350; Creative Commons Attribution 4.0 International 114; NCBI GEO policy
224; none 106 (104 converted: 63 Visium HD, 39 Visium, 2 Xenium).

The root of every one of the 1,884 was read; each is classified
(`profile` in the census), none unknown.

## Store layout (`store-layout.md`)

| Requirement | Their stores |
|---|---|
| 1. Zarr v3 | all 1,884 |
| 2. Consolidated metadata at the root | all 1,884 |
| 3. `spatialdata_attrs` in the root | all 1,884 (`version 0.2`) |
| 5. Byte ranges | the proxy honours them (206), including suffix ranges (`bytes=-N`), which is how the census reads parquet footers |
| 6. A missing key is 404 | yes (a JSON body). Python's default user agent gets **403** for every key; browsers and curl do not |
| 7. CORS | `*`, `Range` allowed, `Content-Range` exposed. The bucket itself has no CORS: the proxy is the only path from a browser |
| 8. Every element into one shared system (`global`) | **1,283 stores**. The other **601** (Visium, Visium HD) name one system per image and none is common to all elements; elements there declare several transforms each, one per system |

Requirement 8 is the one the viewer breaks on: it multiplies all of an
element's transforms (`docs/plans/SDDB.md` §1), so in those 601 stores the
image is placed 148× (Visium) to 198× (Visium HD) too large. The rule the
viewer will follow instead (plan §2), as the census applies it: the target
is `global` in 1,283 stores and the system the table's region is identity
into in 601; every store has exactly **one sample** (no store's default
table annotates regions in different systems).

## Table profile (`table-profile.md`)

| Requirement | Their stores |
|---|---|
| One table is read, the first by name | 1,785 stores hold one table; 99 hold 2 to 5. In **90** Visium HD stores the first by name (`square_002um`, 7.5 M+ bins) is not the curated `Default Table` (`square_008um` in 127, `square_016um` in 1, none named in 9) |
| `obsm/spatial` n × 2, float | n × 2 in every table that has it; **int64 in 654 tables** (float64 in 1,421). The reader converts with `Number()`, so int64 works; the profile's "float" is to be relaxed to "numeric" |
| `obsm/spatial` present | missing in 10 tables, all `cell_segmentations` / `nucleus_segmentations` of 5 Visium HD stores (validator FAIL) |
| `X` CSC, or a CSC copy in `layers/` | `X` is CSR (dense in the spatial-transcriptomics stores) with a CSC copy `layers/csc_counts` in **1,845 tables**; **240 tables have none**: every `square_002um` (99), `square_008um` in 82 stores, the segmentation tables (10), `square_016um_hvg` (1), and 9 earlier versions' `table`. The reader takes any CSC layer of X's shape; the validator now does too (it used to want the name `X_csc`) |
| The default table has a CSC layer | 1,682 stores yes; **82 no** (all Visium HD, the `square_008um` default); 120 name no default (the 109 top-level stores, 9 earlier Visium HD versions, 2 tests) |
| `instance_key`, `region_key` | set on every table |

## Image profile (`image-profile.md`)

| Requirement | Their stores |
|---|---|
| Axes `c, y, x`; `uint8` or `uint16` | every image |
| Chunks ≤ 1024² | SDDB-converted stores: 1024² or smaller, except **13 Visium HD hires images written as one 6000-pixel chunk** (35 MB). Larger chunks are in the older top-level stores only: 4096² uint16 in the 44 Xenium stores (88 images, 32 MB a tile) and 2000² in the 65 Visium stores |
| Levels about halving | the full image has 4 levels in 1,090 stores, 5 in 836, 1 in 8; hires and lowres images are one level each |
| A full-resolution image | **77 stores have none**, only the hires and lowres images (66 Visium HD, 11 Visium): the viewer draws the hires image there, scaled into the target (plan §2) |
| Transforms | `identity` and `scale`, plus **`affine` in 5 Visium HD stores**: a 90° rotation that places the CytAssist image. The viewer rejects it until needed; each of those stores also has a hires image |

## Shapes profile (`shapes-profile.md`)

| Requirement | Their stores |
|---|---|
| One `shapes.parquet` file | every element (2,607), none a directory |
| GeoParquet `geo` metadata | every element |
| Name `cell_boundaries` / `cell_polygons` / `cells` | polygons named `cell_boundaries` in 193 stores (Xenium, Atera; `nucleus_boundaries` beside them); 5 Visium HD stores name theirs `<sample>_cell_segmentations` and are not found |
| Polygons only | **circles** (Point + `radius`) in every store, 1,884 elements: Visium spots, Xenium and Atera `cell_circles`, Visium HD `_for_vitessce`; **squares** (Polygon, the Visium HD bins) 328 elements; polygons 395. Circles and squares are not drawn today (plan step 3) |
| Rows per element | circles: 1,485 elements ≤ 10,000, but 313 over 60,000 (the per-view cap), up to 11.2 M; among them 35 in the `visium/` folder that hold 2 µm bins as circles. Squares: up to 28.1 M; 204 over 250,000 |
| An index column holding the instance id | a named index on 2,160 elements; a `RangeIndex` (instance id = row position) on 313 squares and 134 circles |
| Row groups ordered by position | 125 elements have more than one row group; **699 elements over 100,000 rows are one row group**, so a region cannot be read by range |

## Points (`points-zarr-spec.md`)

193 points elements, all spatialdata's parquet **directory**
(`points.parquet/part.N.parquet`, about 4 M rows a part), in the Xenium and
Atera stores: about **20.5 billion molecules**, up to 2.8 billion in one
store. None is in our tiled points format, so none is streamed. The plan
leaves them for later (conversion with `scripts/points_zarr.py`, or a
row-group reader once a spatial index exists).

## Labels

193 stores hold labels (`cell_labels`, `nucleus_labels`): uint32, 4096²
chunks. Not read (plan: later, on a `TileLayer` base).

## What we ask the maintainers to write differently

To be filed at `theislab/spatialdata-db-feedback` once agreed (not filed
yet), each with the counts above:

1. **A CSC layer on every table the website offers**, the default first: 82
   Visium HD stores' default `square_008um` has none, nor any
   `square_002um` (343 M non-zeros per table, scanned per gene without it).
2. **Image chunks at most 1024²**: the 13 Visium HD hires images written as
   one 6000-pixel chunk; and, if the older top-level stores stay
   published, their 4096² and 2000² chunks.
3. **A spatial index for points** (our tiled points format, or parquet row
   groups by tile with a manifest), so a view can be read without the
   whole element: 193 elements, 20.5 billion molecules.
4. **Shapes row groups ordered by position** for elements over 100,000
   rows (699 elements).
5. **One declared default coordinate system**: `global` in the 601 stores
   that have none, or say that "the system the table's region is identity
   into" is the rule.
6. **The default table first, or named in the store** (it is, in
   `metadata_general["Default Table"]`, but not in the 109 top-level stores
   and 9 earlier versions), and the licence in the store's metadata rather
   than only as a Lamin label (106 stores have neither).
7. **The proxy**: is the 403 for non-browser user agents intended, and is
   there a rate limit? (The census ran 6 requests in flight for 9.5 minutes
   with no 429.)
8. **The 17 stores the website lists beyond the 1,761 converted ones**, and
   whether the 123 older top-level stores are meant to stay listed.
