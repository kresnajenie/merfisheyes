# Store layout

What a folder must look like for MERFISHEYES to open it, wherever it is
served from. The profiles for each element are their own documents:
[table](table-profile.md), [points](../points-zarr-spec.md),
[shapes](shapes-profile.md), [image](image-profile.md); the file beside
the store that links them is [mapping.md](mapping.md). The validator that
checks all of this is `scripts/validate_store.py`.

## The folder

```
<folder>/
  <name>.zarr/            the SpatialData store (Zarr v3)
    zarr.json             root group, with consolidated metadata (required)
    tables/<table>/       the cells (table profile)
    points/<element>/     molecules (points spec)
    shapes/<element>/     cell polygons (shapes profile)
    images/<image>/       a multiscale image (image profile)
    labels/<image>/       (read by nothing yet)
  mapping.json            optional: molecules linked by cell, an image (mapping.md)
  annotations/            optional: cell annotations sidecar
```

A folder opened by URL is not listable, so its store is found by name: the
`cells` entry of `mapping.json`, else `<folder>/<folder>.zarr`
(`lib/spatialdata/resolveUrl.ts`). A dataset in the app's bucket has the
store at its root or in one folder under it.

## Requirements

1. **Zarr v3** (`zarr_format: 3`), written by zarr-python ≥ 3 or anything
   zarrita 0.7 reads. Zarr v2 stores are not opened.
2. **Consolidated metadata at the root** (`consolidated_metadata.metadata`
   in the root `zarr.json`, as zarr-python's `consolidate_metadata` writes
   it). A public URL cannot be listed; the consolidated metadata is the
   only way the viewer learns what the store holds, and it is read once
   (`docs/specs/LOADING.md`). Write it last: it must name every node.
3. **Root attributes** `spatialdata_attrs: {version}` (what tells a store
   from an image: `lib/spatialdata/resolveUrl.ts`).
4. **Element names**: `tables/<name>`, `points/<name>`, `shapes/<name>`,
   `images/<name>`; a store may hold any combination. The viewer shows the
   cells' scene when there is a table, else the molecules' (a linked folder
   named by `mapping.json`), else the image alone.
5. **Sharded arrays** are read by byte range: the server must honour
   `Range` requests (S3 does). A server that ignores them still works, at
   the cost of whole objects.
6. **A missing key** answers 404, or 403 from a public bucket that cannot
   be listed; the viewer treats both as "not there". Any other status is an
   error.
7. **CORS** for a public bucket: `GET` and `HEAD` from the viewer's origin,
   with `Range` allowed in `AllowedHeaders` and `Content-Range`,
   `Content-Length` exposed; `Timing-Allow-Origin` lets the byte counts be
   measured. The merfisheyes-bil bucket allows the production origin and
   the standard local ports only.
8. **Transforms**: every element declares `coordinateTransformations` into
   one shared coordinate system (`global`). The viewer draws in the cells'
   intrinsic coordinates (the table's `obsm/spatial`); an image is placed
   by undoing the shapes element's transform, else the points element's
   (`lib/spatialdata/image.ts`). Where the points and shapes declare
   different transforms to `global` (the atera store does), the shapes'
   wins for the image and the points are drawn as written.

## What the validator checks

`python scripts/validate_store.py <store.zarr> [--mapping]` reads the store
from disk and reports, per element, the profile's checks; it exits non-zero
on the first failure of a MUST. It is run in CI on the tiny fixtures
(`tests/python/`) and by hand on the real stores.
