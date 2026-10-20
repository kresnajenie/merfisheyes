# mapping.json

A JSON file beside a cell dataset (a merfisheyes cell folder or a
SpatialData store) that names what is drawn with the cells: which molecule
dataset each cell links to, and an image under them. Schema:
[schemas/mapping.schema.json](schemas/mapping.schema.json). Read by
`lib/viewer/resolve-overlays.ts` (public folders), the app's
`/api/datasets/{id}/mapping` route (uploaded datasets) and
`lib/spatialdata/resolveUrl.ts` (to find a folder's store).

```json
{
  "cells": "skin.zarr",
  "linkColumn": "__all__",
  "links": { "__all__": "https://bucket/path/molecules" },
  "source": "s3",
  "image": "https://bucket/path/dapi.zarr",
  "imageSource": "s3"
}
```

| Key | Meaning |
|---|---|
| `cells` | The store in this folder (a name inside it, or a full URL). Only needed when the store is not named after the folder. |
| `linkColumn` | `"__all__"`: every cell links to the one dataset in `links.__all__`, drawn over the cells. Any other value: the name of a cluster column; a cell's value picks its entry in `links`, opened in a split panel on right-click. |
| `links` | Column value → a molecule dataset: a public folder's URL or a name inside this folder (`source: "s3"`, the default), or an app dataset id (`source: "app"`). |
| `source` | How `links` are opened: `"s3"` (public URL, `SingleMoleculeDataset.fromCustomS3`) or `"app"` (presigned API, `fromS3`). |
| `image` | An OME-Zarr image drawn under the cells: a public URL or a name inside this folder, or with `imageSource: "app"` an image dataset's id. Its transforms must place it in the cells' coordinates; it is drawn only for cells with raw (un-normalized) coordinates. |
| `imageSource` | `"s3"` (default) or `"app"`. |
| `imageDatasetId` | The image's dataset id when it is registered (for the overlay picker). |

Rules:

- A mapping with `links` MUST have `linkColumn`; one with neither links no
  molecules and may name an image alone.
- A `links` value or `image` without a scheme (`molecules`, `dapi.zarr`,
  `./dapi.zarr`) is a name inside the folder the mapping was fetched from
  (`resolveMappingPaths` in `lib/config/dataset-links.ts`), so a folder
  prepared once by `scripts/prepare_folder.py` opens from any bucket; app
  ids are never resolved.
- A SpatialData store whose own points element is drawn over the cells
  does not need `links`: the store's points stand in for a linked folder
  (`lib/scene/registry.tsx`, `openElements`).
- The linked folder is the legacy per-gene format (`manifest.json.gz`,
  `genes/<gene>.bin.gz`); `scripts/legacy_to_points.py` converts it into a
  points element inside the store, after which the link is not needed.
