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

## Opening a store

**By public URL** — a store in a public-read bucket:

```
/spatialdata-viewer/from-s3?url=https://<bucket>.s3.<region>.amazonaws.com/<folder>/<name>.zarr
```

**By dataset id** — a store under `datasets/{id}/` in the app's own bucket,
with a `Dataset` row (`COMPLETE`):

```
/spatialdata-viewer/{id}
```

The store is found by its root `zarr.json`, so it can sit at
`datasets/{id}/` or in a subfolder such as `datasets/{id}/data.zarr/`.

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
└── mapping.json          links the cells to molecules/
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

# 3. Link the cells to the molecules
echo "{\"linkColumn\":\"__all__\",\"links\":{\"__all__\":\"$URL/molecules\"}}" > mapping.json
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

## What a store needs

- **Zarr v3** (`zarr.json`). Zarr v2 SpatialData stores are not read.
- **A table** under `tables/` with `obsm/spatial` (or `obsm/X_spatial`).
  With several tables, the one named `table` is used, else the first.
- **Expression that can be sliced by gene.** `X` is CSR in these stores, so
  genes are read from a CSC copy in `layers/` (e.g. `layers/X_csc`): one gene
  is ~17 small chunk reads. A dense or CSC `X` also works. A CSR `X` with no
  CSC copy still works but reads the whole matrix for every gene.
- **Consolidated metadata** in the root `zarr.json`, for stores opened by
  URL: a plain URL can't be listed, so this is how the viewer finds the
  store's elements and columns. Stores opened by dataset id are listed
  through S3 instead.
- **An image** (optional) with axes `c, y, x` and `uint8` or `uint16`
  pixels. Its transforms and those of the shapes/points element may be
  `identity`, `scale` and `translation`; others are ignored with a console
  warning.
- **CORS** on the bucket for the origin the app is served from, allowing
  `GET`.

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
  channel's visibility, colour and brightness. Only the first channel is on
  to start with; its brightness is set from the 99.9th percentile of the
  coarsest level.
- Chunks are whatever the store was written with. The test store uses
  4096 × 4096 chunks of ~19 MB at full resolution, and a zoomed-in view
  needs two to four per visible channel, so that is the slow step on a
  slow connection.

## Not supported yet

- Cell boundaries (`shapes/`).
- Reading transcripts straight from the store's spatially tiled Parquet
  (they are converted to per-gene files instead).
- Labels, multiple images, multiple coordinate systems.
- Uploading a store through the app; stores are synced to S3 by hand.

## Code

| Path | What |
| --- | --- |
| `app/spatialdata-viewer/` | The two routes |
| `components/spatialdata-viewer.tsx` | Loads the store, mounts the scene and the image layer |
| `components/spatialdata-image-controls.tsx` | Image button and channel panel |
| `lib/spatialdata/store.ts` | Opening a store by dataset id or URL |
| `lib/spatialdata/SpatialDataTableAdapter.ts` | Table → cells, genes, columns, annotations, DE stats |
| `lib/spatialdata/image.ts` | OME-Zarr multiscale metadata and placement |
| `lib/spatialdata/ImageLayer.ts` | Tiled image layer in the Three.js scene |
| `scripts/spatialdata_points_to_sm.py` | `points/` → per-gene molecule files |
| `scripts/spatialdata_cluster.py` | Leiden clusters + DE stats sidecar |
| `tests/unit/spatialdata-table-adapter.test.ts` | Adapter tests on a tiny fixture (`scripts/testdata/make-tiny-zarr.py`) |

Zarr v3 string arrays need zarrita ≥ 0.7. The h5ad-zarr path is pinned to
0.5.1, so the newer one is installed under the alias `zarrita-v3` and used
only here.
