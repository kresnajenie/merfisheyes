# SpatialData stores

MERFISH Eyes can open a [SpatialData](https://spatialdata.scverse.org) zarr
store directly from S3 at `/spatialdata-viewer`, reading it lazily and **as
written** — the store is never converted or modified. Anything the viewer
needs that the store doesn't have (per-gene molecule files, clusters) is
written **beside** it by the scripts below.

It was built against the Xenium skin store from the Celldega
[spatial-tiling proof of concept](https://github.com/cornhundred/POC_SpatialData_spatial-tile_viz_with_Celldega)
(`skin_adapt_v2.zarr`: 112,551 cells × 5,006 genes, 74M transcripts, a
4-channel morphology image).

## What is shown

| Layer | Read from | How |
| --- | --- | --- |
| Cells | `tables/<name>` in the store | Coordinates, gene list and one cell column up front; expression per gene on demand |
| Image | `images/<name>` in the store | The pyramid level matching the zoom, only the chunks in view |
| Molecules | `molecules/` beside the store | Per-gene files, converted once from `points/<name>` |
| Clusters + DEGs | `annotations/` beside the store | Leiden clusters and their DE stats, computed once from the table |

Cell boundaries (`shapes/<name>`) are not drawn yet.

A store does not need all three. Whatever it has is shown:

| The store has | Shown as |
| --- | --- |
| Cells (with or without molecules and an image) | The single cell scene: molecules over the cells, image under them |
| Molecules but no cells (with or without an image) | The single molecule scene, image under the molecules |
| Only an image | The image on its own |
| None of them | An error |

"Has molecules" means the converted `molecules/` folder is linked by a
`mapping.json` beside the store (see below); the store's own `points/` are
not read directly.

## Opening a store

**By public URL** — a store in a public-read bucket:

```
/spatialdata-viewer/from-s3?url=https://<bucket>.s3.<region>.amazonaws.com/<folder>/<name>.zarr
```

The URL may also be the **folder that holds the store**
(`…/<folder>/`). A public URL can't be listed, so the store is found by
name: the `cells` entry of the folder's `mapping.json`, else a store named
after the folder (`<folder>/<folder>.zarr`). Either way everything beside the
store is overlaid, and the link reopens all of it.

The homepage's **Load from S3** box takes the same URLs — a store, the folder
holding one, or an OME-Zarr image — as well as merfisheyes folders.

A store opened by URL that nobody owns shows the **Add to my account**
banner, like any other dataset opened by URL. Nothing is copied: the row
(`datasetType: "spatialdata"`, `ingestSource: "s3_registered"`) remembers the
store's URL, it is listed under **Your datasets**, and it opens at the same
`from-s3` link.

**By dataset id** — a store under `datasets/{id}/` in the app's own bucket,
with a `Dataset` row (`COMPLETE`):

```
/spatialdata-viewer/{id}
```

The store is found by its root `zarr.json`, so it can sit at
`datasets/{id}/` or in a subfolder such as `datasets/{id}/data.zarr/`.

**From your computer** — drop the store, or a folder that holds it, on the
homepage's single-cell **Folder** card in *Preview in browser*. See the next
section.

## Opening already-processed data from disk

A folder whose parts are already in a format the viewer reads opens as it is:
nothing is converted, and each part is read from the dropped files a piece at
a time — a table column, one gene, the image tiles in view — so its size does
not matter. Two layouts are accepted, and any part may be missing:

```
sample/                          sample/
├── <store>.zarr/   cells+image  ├── <folder>/        cells (merfisheyes output)
├── <folder>/       molecules    ├── <folder>/        molecules (merfisheyes output)
└── annotations/    optional     └── <image>.zarr/    an OME-Zarr image
```

| Part | Recognised by |
| --- | --- |
| Cells (SpatialData) | a `zarr.json` plus `tables/`, `images/` or another element group |
| Cells (table alone) | a `zarr.json` with `obs/` and `var/` beside it: `tables/<table>/` or `tables/` out of a store, or any AnnData written to Zarr v3 |
| Cells (merfisheyes) | `manifest.json`, `expr/index.json`, `coords/spatial.bin.gz` |
| Molecules | `manifest.json.gz` (or `manifest.json`) and `genes/*.bin.gz` |
| Image | a `zarr.json` describing an OME-Zarr image with axes `c, y, x` |

- Names don't matter, except `annotations/`. A part is a folder directly
  inside the dropped one, or the dropped folder itself (an image folder
  dropped alone opens in the image viewer).
- With cells, the molecules go over them and the image under them (a store's
  own image, else the standalone one). Without cells the molecules are the
  scene, image under them. With several candidates for a part the first by
  name is used and a message names the rest.
- The preview runs in the ordinary viewers — `/viewer/from-local`,
  `/sm-viewer/from-local`, `/image-viewer/from-local` — with the extras handed
  over through `lib/stores/localBundleStore.ts`.
- Raw cells (H5AD, Xenium, MERSCOPE) beside processed parts are shown on
  their own, with a note: they have to be parsed first.

**Upload & Save** from such a preview sends every part as it is, each as its
own dataset in one project: a merfisheyes cell folder through the pre-chunked
upload, the molecules through the pre-chunked molecule upload, an image as an
image dataset, and a SpatialData store whole (`datasetType: "spatialdata"`,
its `annotations/` beside it), which then opens at `/spatialdata-viewer/{id}`.
The cells' `mapping.json` is written with both links, so the saved link opens
overlaid. Molecules and an image saved without cells are grouped in a project
but not overlaid: only cells carry overlays.

## Layout on S3

```
<folder>/
├── <name>.zarr/          the SpatialData store, exactly as written
├── molecules/            from scripts/spatialdata_points_to_sm.py
│   ├── manifest.json.gz
│   └── genes/{gene}.bin.gz
├── annotations/          from scripts/spatialdata_cluster.py
│   ├── annotations.json
│   ├── leiden.bin
│   └── de/leiden.bin.gz
└── mapping.json          names the store ("cells") and links it to molecules/
```

`molecules/`, `annotations/` and `mapping.json` are all optional; without
them the viewer shows the cells and the image only.

## Preparing a store

Both scripts need the Python environment of the other `scripts/` plus
`scanpy`, `leidenalg` and `igraph` for clustering.

```bash
STORE=/path/to/name.zarr
DEST=s3://<bucket>/<folder>
URL=https://<bucket>.s3.<region>.amazonaws.com/<folder>

# 1. The store itself, unchanged
aws s3 sync $STORE $DEST/name.zarr/ --exclude '*.DS_Store'

# 2. Molecules: points/<name> → per-gene files (74M transcripts: ~1 min)
python scripts/spatialdata_points_to_sm.py $STORE molecules/ --workers 16 --drop-z
aws s3 sync molecules/ $DEST/molecules/

# 3. Name the store and link its cells to the molecules. "cells" is what lets
#    the folder URL be opened on its own; it is not needed when the store is
#    named after the folder.
echo "{\"cells\":\"name.zarr\",\"linkColumn\":\"__all__\",\"links\":{\"__all__\":\"$URL/molecules\"}}" > mapping.json
aws s3 cp mapping.json $DEST/mapping.json

# 4. Clusters and DEGs (112k cells: ~1 min)
python scripts/spatialdata_cluster.py $STORE annotations/ --resolution 0.5
aws s3 sync annotations/ $DEST/annotations/
```

For a store opened by dataset id, the same three sidecars go in the folder
that holds the store, and `mapping.json` links to the molecule dataset's id
with `"source": "app"` (see `lib/ingest/overlay.ts`).

### `--drop-z`

Without it molecules keep their own z (Xenium: roughly 10–44 µm) while the
cells and the image sit at z = 0. The scene uses a perspective camera, so
molecules above the plane slide against it as you pan. `--drop-z` writes
them flat on that plane; leave it off only for a 3D view of the molecules
on their own.

### Clustering

`spatialdata_cluster.py` runs normalize_total → log1p → PCA → neighbours →
Leiden (`--resolution`, `--n-pcs`, `--n-neighbors`) and writes one category
code per cell, in table order. The DE stats are per-cluster means of the
**log-normalized** matrix plus the fraction of cells expressing each gene —
the DEG panel un-logs the means for fold change. (Raw Xenium counts are low
enough that the panel would mistake their means for log values.)

Clusters are numbered, not named, and no UMAP is written.

The sidecar can also carry DEG stats for columns the table already has
(e.g. its own `leiden` or cell type labels): list them in `de_stats`, leave
them out of `columns`, and write `de/<column>.bin.gz` in the same format —
per group, each gene's mean of the log-normalized matrix and the fraction of
cells expressing it, with the groups named as the column's values.

```json
{ "num_cells": 147402, "num_genes": 275, "columns": {}, "de_stats": ["leiden", "celltype"] }
```

## What a store needs

- **Zarr v3** (`zarr.json`). Zarr v2 SpatialData stores are not read.
- **A table** (optional) under `tables/` with `obsm/spatial` (or
  `obsm/X_spatial`). With several tables, the one named `table` is used,
  else the first. It is read as anndata writes it, old or new: plain or
  sharded arrays, and plain or nullable string columns (a missing string
  reads as empty).
- **Expression that can be sliced by gene.** `X` is CSR in these stores, so
  genes are read from a CSC copy in `layers/` (e.g. `layers/X_csc`): one gene
  is ~17 small chunk reads. A dense or CSC `X` also works; a dense `X` written
  one gene per chunk (`chunks=(n_cells, 1)`) is read in one request per gene.
  A CSR `X` with no CSC copy still works but reads the whole matrix for every
  gene.
- **Consolidated metadata** in the root `zarr.json`, for stores opened by
  URL: a plain URL can't be listed, so this is how the viewer finds the
  store's elements and columns. Stores opened by dataset id are listed
  through S3 instead.
- **An image** (optional) with axes `c, y, x` and `uint8` or `uint16`
  pixels. Its transforms and those of the shapes/points element may be
  `identity`, `scale`, `translation` and a `sequence` of those (how the
  `spatialdata` library writes scale-then-translation); others are ignored with a console
  warning.
- **CORS** on the bucket for the origin the app is served from, allowing
  `GET`. A store with sharded arrays (anndata ≥ 0.12 writes them by default)
  is read by byte range, and the browser asks permission for that first, so
  the bucket must also allow the `Range` request header:

  ```json
  [
    {
      "AllowedOrigins": ["https://<site>"],
      "AllowedMethods": ["GET", "HEAD"],
      "AllowedHeaders": ["*"]
    }
  ]
  ```

  Without it a sharded store fails with "Failed to fetch". A server that
  ignores `Range` and returns the whole object still works; the part asked
  for is cut out in the browser.

The SpatialData join-key columns (`instance_key`, `region_key`) are hidden
from the cell-column picker.

## Coordinates

The scene is in the units of the cells: `obsm/spatial`, which for Xenium is
microns — the intrinsic coordinates of the shapes and points elements. The
image is placed in those units by undoing the shapes (or points) element's
transform to the shared coordinate system (for Xenium, pixels = microns ×
4.7059).

## Image

- One pyramid level is chosen per view: the finest that is no finer than
  the screen. The coarsest level stays loaded as a backdrop, so the view
  sharpens as chunks arrive rather than blanking.
- The **Image** button (top left) shows or hides the image and sets each
  channel's visibility, colour and display window: a lower bound (`min`,
  values at or below it draw black) and an upper bound (`max`, values at or
  above it draw full colour). Only the first channel is on to start with;
  its window starts at 0 to the 99.9th percentile of the coarsest level.
- What is set there is kept in the link (`img=`, or `rimg=` in the right
  split panel) once something is changed, so a copied link opens looking the
  same. An untouched image leaves the link alone.
- Chunks are whatever the store was written with. The test store uses
  4096 × 4096 chunks of ~19 MB at full resolution, and a zoomed-in view
  needs two to four per visible channel, so that is the slow step on a
  slow connection.

## An image is a dataset

An OME-Zarr image is a dataset of its own (`datasetType: "image"`), owned,
listed and deleted like a single cell or single molecule one. There are two
ways to have one:

- **By URL.** Paste its URL into the homepage's **Load from S3**, or open
  `/image-viewer/from-s3?url=…`. If nobody
  owns it yet a banner offers **Add to my account** — the same banner the
  other viewers show for a dataset opened by URL. Nothing is copied; the row
  remembers the URL (`ingestSource: "s3_registered"`).
- **By upload.** As part of a combined upload (below). The image's files go
  into the app's own bucket under `datasets/{id}/` and it opens at
  `/image-viewer/{id}`.

Either way it then shows under **Your datasets** (filter *Type → Image*) and
in the overlay dialog.

## An image under a dataset's cells

A single cell dataset that is not a SpatialData store can still have an image
under its cells. Its owner opens the viewer's overlay button (or the dataset's
page) and, under **Image overlay**, picks one of their images from the list.
**Add by URL** takes an image that is not in the account yet: it is opened
first, so one the viewer could not read is refused, then added to the account
and attached. (An image someone else already owns is attached by its URL
without being added.)

The image must be an OME-Zarr multiscale image group (the group itself, not a
SpatialData store): Zarr v3, axes `c, y, x`, `uint8` or `uint16`. Its own
`scale` / `translation` must place it in the cells' coordinates, e.g.
`scale: [1, 0.10833, 0.10833]` for an image at 0.10833 µm per pixel whose
pixel (0, 0) is at the cells' origin. It is drawn only for datasets with raw
(not normalized) coordinates, and gets the same **Image** button as a store's
own image. Rotations and flips are not applied, so they must already be in
the pixels.

### In `mapping.json`

The choice is stored in the dataset's `mapping.json`, next to the molecule
overlay and independent of it: either may be there without the other, and
changing one keeps the other.

```json
{
  "linkColumn": "__all__",
  "links": { "__all__": "sm_abc123" },
  "source": "app",
  "image": "ds_def456",
  "imageSource": "app",
  "imageDatasetId": "ds_def456"
}
```

| Field | Meaning |
| --- | --- |
| `image` | The public URL of the image, or — with `imageSource: "app"` — the id of an image dataset in the app's bucket |
| `imageSource` | `"app"`: open `image` by dataset id through the presigned API. `"s3"` or absent: `image` is a URL |
| `imageDatasetId` | The image's dataset id whenever it is a registered dataset, so the dialog can name it |

A `mapping.json` beside a dataset on a public bucket (opened with
`/viewer/from-s3`) only needs `image` with a URL.

### Overlay API

Both overlays of an app dataset are edited through one route. Reading needs
no sign-in; writing is for the dataset's owner (`401` signed out, `404` for
anyone else).

| Request | Does |
| --- | --- |
| `GET /api/datasets/{id}/overlay` | `{ smDatasetId, image, imageDatasetId }`, each `null` when unset |
| `PATCH …/overlay` with `{ "imageDatasetId": "<id>" }` | Attaches one of your own image datasets (`400` if it is not yours or not an image) |
| `PATCH …/overlay` with `{ "image": "<http(s) URL>" }` | Attaches an image by URL, without it being a dataset |
| `PATCH …/overlay` with `{ "image": null }` | Removes the image |
| `POST` / `DELETE …/overlay` | Attach / remove the molecule overlay |

Deleting an image dataset removes it from every dataset of yours that drew it.

## Cells, molecules and an image from one folder

With **Upload & process on server** chosen on the homepage, a picked folder
may hold all three:

```
sample/
├── cells.h5ad            or a Xenium / MERSCOPE export
├── transcripts.parquet   (optional) the molecules
└── dapi.zarr/            (optional) an OME-Zarr image, any name ending .zarr
```

The dialog that names the upload offers each extra it found as a tick box.
On confirm:

1. The molecules and the image each start uploading as their own dataset
   (“… (molecules)”, “… (image)”). The image goes up as it is — it needs no
   processing — and is viewable as soon as its last file lands.
2. Once their ids are known the cells start (“… (cells)”), carrying those ids
   (`processingParams.linkedSmDatasetId` / `linkedImageDatasetId`).
3. All of them are grouped in one project named after the upload.
4. When the cells finish processing, `mapping.json` is written with both
   links. Until then the mapping route answers from the ids on the row.

The cells' link — on their upload bar and in their “ready” email — opens
their viewer, which draws the molecules over the cells and the image under
them with no further step.

A `.zarr` folder is only offered when its root `zarr.json` describes an image
the viewer can draw (axes `c, y, x`). The image is uploaded under one
refreshable POST policy, eight files at a time, so a multi-gigabyte image is
not cut off when a signature expires; no single file may exceed 5 GB (what
S3 accepts in one POST — a sharded array is one file however large it is).

### Writing an image

Any tool that writes OME-Zarr 0.5 (Zarr v3) with `c, y, x` axes will do. With
zarr-python ≥ 3, for one stitched plane:

```python
import zarr

plane = ...            # 2D uint8/uint16 array: rows are y, columns are x
um_per_pixel = 0.10833
x0, y0 = 0.0, 0.0      # where the centre of pixel (0, 0) sits among the cells

root = zarr.open_group("dapi.zarr", mode="w", zarr_format=3)
datasets, level = [], plane
for i in range(12):
    root.create_array(str(i), data=level[None], chunks=(1, 1024, 1024))
    size = um_per_pixel * 2**i
    shift = (size - um_per_pixel) / 2   # a coarser pixel's centre moves out
    datasets.append({
        "path": str(i),
        "coordinateTransformations": [
            {"type": "scale", "scale": [1, size, size]},
            {"type": "translation", "translation": [0, y0 + shift, x0 + shift]},
        ],
    })
    if max(level.shape) <= 1024:
        break
    h, w = level.shape[0] // 2 * 2, level.shape[1] // 2 * 2
    level = level[:h, :w].reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3)).astype(plane.dtype)

root.attrs["ome"] = {
    "version": "0.5",
    "multiscales": [{
        "axes": [{"name": "c", "type": "channel"},
                 {"name": "y", "type": "space", "unit": "micrometer"},
                 {"name": "x", "type": "space", "unit": "micrometer"}],
        "datasets": datasets,
    }],
    "omero": {"channels": [{"label": "DAPI"}]},
}
```

Keep halving until a level fits in one or two chunks: the coarsest level is
loaded whole as the backdrop. Upload the folder to a public-read bucket whose
CORS allows the site, and use the folder's URL.

## An image on its own

An image opens without any cells or molecules, by URL or — when it was
uploaded — by dataset id:

```
/image-viewer/from-s3?url=https://<bucket>.s3.<region>.amazonaws.com/<folder>/dapi.zarr
/image-viewer/{id}
```

It is drawn in its own coordinates with the same layer and **Image** button,
so it looks the same here as it does under a dataset's cells. From the
homepage, drop the image folder on the single-cell **Folder** card (*Preview
in browser*) or paste its URL into **Load from S3**. The URL may also be an
image inside a SpatialData store (`…/<name>.zarr/images/<image>`).

## Not supported yet

- Cell boundaries (`shapes/`).
- Reading transcripts straight from the store's spatially tiled Parquet
  (they are converted to per-gene files instead).
- Labels, multiple images, multiple coordinate systems.
- Uploading a store by itself from the homepage's server mode; a store is
  uploaded with *Upload & Save* from its local preview, or synced to S3 by
  hand.

## Backlog

- **View the parts that are ready while the rest convert.** Today, raw cells
  (H5AD, Xenium, MERSCOPE) dropped beside an already-processed image or
  molecule folder are previewed on their own. The aim: open the processed
  parts at once — the image, the molecules — and add the cells to the same
  scene when their conversion finishes, in the browser or on the server,
  instead of making the whole folder wait for its slowest part.
- An image under a saved molecule dataset that has no cells (only cells carry
  overlays today).

## Code

| Path | What |
| --- | --- |
| `app/spatialdata-viewer/` | The two routes |
| `components/spatialdata-viewer.tsx` | Loads the store, mounts the scene and the image layer |
| `components/spatialdata-image-controls.tsx` | Image button and channel panel |
| `components/image-viewer.tsx`, `app/image-viewer/` | An image on its own, by URL |
| `lib/spatialdata/store.ts` | Opening a store by dataset id or URL |
| `lib/spatialdata/localStore.ts` | A zarr store over dropped files; opening a store from them |
| `lib/local/classify-bundle.ts`, `open-bundle.ts`, `save-bundle.ts` | Recognising, opening and saving a dropped folder of processed parts |
| `components/bundle-upload-modal.tsx` | *Upload & Save* for such a folder |
| `lib/storage/fetchRange.ts` | Byte-range reads, for sharded arrays |
| `lib/spatialdata/SpatialDataTableAdapter.ts` | Table → cells, genes, columns, annotations, DE stats |
| `lib/spatialdata/loadSpatialDataMolecules.ts` | The linked molecules of a store that has no cells |
| `lib/spatialdata/image.ts` | OME-Zarr multiscale metadata and placement |
| `lib/spatialdata/ImageLayer.ts` | Tiled image layer in the Three.js scene |
| `lib/spatialdata/imageUrlState.ts` | An image's display settings in the link (`img=`) |
| `lib/spatialdata/resolveUrl.ts` | What a URL is: a store, an image, or a folder holding a store |
| `components/three-scene.tsx` | Draws a dataset's `mapping.json` image under its cells |
| `components/overlay-manager.tsx`, `components/image-overlay-picker.tsx` | The owner's overlay dialog: molecules, and the image browser |
| `lib/datasets/kind.ts` | The three dataset kinds and which viewer each opens in |
| `lib/datasets/register-or-claim.ts` | Registering a URL as a dataset; recognises an OME-Zarr image and a SpatialData store |
| `lib/upload/uploadImageToS3.ts` | Uploading an image folder as an image dataset |
| `lib/ingest/classify-folder.ts` | Finding cells, molecules and image folders in a picked folder |
| `app/api/datasets/[datasetId]/overlay/`, `lib/ingest/overlay.ts` | Overlay API and the `mapping.json` it writes |
| `scripts/spatialdata_points_to_sm.py` | `points/` → per-gene molecule files |
| `scripts/spatialdata_cluster.py` | Leiden clusters + DE stats sidecar |

## Tests

| Where | What | Runs in CI as |
| --- | --- | --- |
| `tests/unit/spatialdata-table-adapter.test.ts` | The table adapter, on both tiny stores (plain, and sharded with nullable strings) | Unit tests |
| `tests/unit/fetch-range.test.ts` | Byte-range reads, including a server that ignores `Range` | Unit tests |
| `tests/unit/ome-zarr-image.test.ts` | Image metadata and placement, including a `sequence` transform | Unit tests |
| `tests/unit/spatialdata-url.test.ts` | What a URL is (store, image, folder by `mapping.json` or by name); image settings in a link | Unit tests |
| `e2e/tests/viewer/spatialdata.spec.ts` | Every combination of cells, molecules and image in a store opened by URL; the sharded store; an image on its own; **Load from S3** with a folder, a store and an image; image settings kept in the link | E2E quick (`@quick`) |
| `e2e/tests/upload/image-overlay.spec.ts` | Adding an image, and a store, to an account by URL; picking the image in the overlay dialog, the channel bounds, who may change it, deleting it | E2E upload (`@upload`) |
| `e2e/tests/upload/combined-upload.spec.ts` | One folder of cells, molecules and an image becoming three linked datasets in a project; the uploaded image opening by id | E2E upload (`@upload`) |
| `tests/unit/classify-folder.test.ts`, `tests/unit/image-dataset.test.ts` | Image-folder detection, dataset kinds, recognising an OME-Zarr image | Unit tests |
| `tests/unit/local-bundle.test.ts` | Recognising the parts of a dropped folder, including a table on its own; reading zarr from local files | Unit tests |
| `e2e/tests/viewer/local-bundle.spec.ts` | Dropping an image, both layouts and their subsets, a table on its own; duplicates; raw cells beside processed parts; normalized cells | E2E quick (`@quick`) |
| `e2e/tests/upload/local-bundle-save.spec.ts` | Saving each layout, a table on its own and an image from its preview, and reopening the saved link overlaid | E2E upload (`@upload`) |

The browser tests serve the committed fixtures from a local stand-in for a
public bucket (`e2e/helpers/fixture-server.ts`), which answers CORS
preflights and byte ranges like S3. One store on disk covers every
combination: `e2e/helpers/spatialdata.ts` trims the consolidated metadata the
viewer lists a store by.

The fixtures are `tests/data/tiny/zarr/` (two stores, 30 cells × 8 genes and
a 64 × 64 image) and `tests/data/tiny/single-molecule/chunked/`. Regenerate
the stores with `python scripts/testdata/make-tiny-zarr.py` (needs
anndata ≥ 0.12 and zarr ≥ 3); see [TESTING.md](../../TESTING.md) for the rest.

Zarr v3 string arrays need zarrita ≥ 0.7. The h5ad-zarr path is pinned to
0.5.1, so the newer one is installed under the alias `zarrita-v3` and used
only here.
