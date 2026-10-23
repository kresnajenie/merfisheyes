# Changelog

All notable changes to MERFISHEYES are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
See [docs/RELEASING.md](docs/RELEASING.md) for how a release is cut.

## [Unreleased]

### Added

- **"Upload & process on server" can run `scripts/prepare.py`** (steps 2–3
  of `docs/plans/ONE-PIPELINE.md`), behind `NEXT_PUBLIC_INGEST_PIPELINE=prepare`;
  unset, nothing changes. The worker runs the same script a person runs, in
  its own environment with `scripts/requirements.txt`, and writes one
  SpatialData store per upload: the cells and the molecules beside them as
  one dataset (no separate molecules dataset or project), opened in the
  SpatialData viewer; a molecule table alone becomes a store of points. The
  upload sends what `prepare.py` reads — a Xenium export's seven files, a
  MERSCOPE region's tables, boundaries and transcripts — and the images only
  when "Also upload the images (N GB)" is ticked (it starts unticked). The
  compute tier is chosen from the upload's size instead of found by running
  out of memory. A molecule table whose columns are not recognised asks for
  them, listing the columns in its schema. On the Xenium lung and breast
  exports the worker's store and a local `prepare.py` run's are identical,
  file for file.

- **DEG stats live inside the table.** Every categorical column's per-gene
  mean (log-normalised) and fraction expressing are written into the
  SpatialData table itself — `varm/de_mean__<column>`, `varm/de_pct__<column>`,
  `uns/merfisheyes_de` (`scripts/de_stats.py`) — by `prepare.py` and
  `spatialdata_cluster.py`, so a store needs no `annotations/` folder beside
  it; the viewer's DEG panel reads them (older stores' `annotations/` still
  work).

- **`scripts/prepare.py`: one command to prepare your own data** (step 1 of
  `docs/plans/ONE-PIPELINE.md`). `python scripts/prepare.py INPUT OUT`
  recognises an `.h5ad` / AnnData Zarr, a Xenium export, a MERSCOPE region,
  a Visium output, a molecule table (a lab's own `Gene, X, Y, Z` columns too)
  or a folder of parts, and writes one SpatialData store with `mapping.json`
  and `stats.json`, never changing INPUT. Images are written when the input
  has them (`--no-images` to skip); `--labels` adds a MapMyCells CSV;
  control probes are left out of every table; an `uns` entry the Zarr writer
  refuses is left out with a warning; cells are placed from `obs` coordinate
  columns when `obsm` has none; Xenium no longer keeps a second (parquet)
  copy of the molecules. Packages pinned in `scripts/requirements.txt`
  (`requirements-prepare.txt` now points there); shared rules in
  `scripts/meyes_common.py`. The main README says where it is and how to run
  it. Exit codes 0 / 2 / 3 / 4.

- **Every kind of SpatialData-DB store opens: CODEX, PhenoCycler,
  spatial-CUT&Tag, spatial-ATAC-seq, MALDI-MSI and CosMx beside the RNA
  ones.** Measured on the root metadata of one store per technology
  (`tests/data/sddb/roots`): a table's features come from the index column
  `var` declares (`formula_adduct` in a MALDI table, `gene` in a CosMx
  one), not only `var/_index`; a dense `X` (every image channel a feature,
  as CODEX writes it) is read one feature at a time; a table with no
  `obsm/spatial` (CosMx) draws its cells at the centres of the shapes it
  annotates, read from the Point `<region>_for_vitessce` element when there
  is one (35 MB for 1.46 M cells, against 265 MB of polygons) and matched
  by instance id; an image placed by an `affine` that rotates (a MALDI
  optical image, 91° in one store) is drawn as is, the levels staying in
  the image's own frame and the placement being one matrix on the layer;
  and an image of any number of channels (84 in a CODEX store) compiles,
  the layer showing at most 8 at once through fixed shader slots. The
  validator (`scripts/validate_store.py`) applies the same rules, and the
  profiles say so (`docs/specs/table-profile.md`, `image-profile.md`,
  `store-layout.md`). Fixtures shaped like the real stores, built by the
  `spatialdata` library: `scripts/testdata/make-tiny-modalities.py`.
- **Processed data on "Upload & process on server".** A SpatialData store,
  separate stores, an OME-Zarr image or merfisheyes output folders dropped
  on the server card are uploaded as they are, the same as Upload & Save
  from a preview, without opening them first; it used to say "No
  single-cell data found". Raw exports still go to the worker.
- **Separate SpatialData parts become one store.** A store without a table
  (tiled molecules, an image) and a standalone OME-Zarr image beside a store
  or table join it as its elements, in the preview and when saved; only the
  root `zarr.json` is rewritten to list them.
- **Drop matrix tests**: every combination of a table (CSR `X` + CSC layer,
  or CSC `X`), an image and tiled molecules, in four layouts, dropped through
  the drop zone (`e2e/tests/viewer/drop-matrix.spec.ts`), saved and sent to
  the server card (`e2e/tests/upload/drop-matrix-save.spec.ts`), and timed on
  real stores (`e2e/perf/real-drop.spec.ts`) and real uploads, prepared or
  raw (`e2e/perf/real-upload.spec.ts`); tiny stores whose `X` is CSC
  (`scripts/testdata/make-tiny-csc-x.py`) and whose indexes have names
  (`make-tiny-named-index.py`).

### Fixed

- A store's tiled molecules without a table were dropped: a molecules-only
  store bounced back to the homepage, and one with an image showed the image
  alone. They now open in the molecule viewer, dropped or saved.
- A molecules-only store beside a table store could take the table's place
  (the first by name won), leaving no cells.
- Several folders dragged in at once overwrote one another's root metadata,
  keeping only one of them.
- A store holding only molecules opened with the camera on an empty patch
  of the plane, so none were drawn; the molecule scene now frames them.
- The local `@upload` test stack no longer calls AWS Batch: with a custom S3
  endpoint no job is submitted (it was reaching Batch with a developer's
  staging settings, refused only for the wrong keys).
- **A lab's own molecule table processes on the server.** The worker read
  only Xenium's (`x_location`) and MERSCOPE's (`global_x`) column names and
  failed on anything else (`Gene, X, Y, Z`); it now takes a gene column and
  x / y (/ z) by name in any case, and a failure names the columns it found.
- **A molecule table beside the cells is found by its columns** on the server
  card, not only when "transcript" is in its name.
- **Gene names stored as bytes** (Xenium Onboard Analysis 1.x `transcripts.parquet`)
  crashed the molecule processor ("cannot use a string pattern on a
  bytes-like object"); they are read as text.
- **A failed processing job says why.** When the worker's report did not
  arrive, the dataset only said "container exited with exit code 1"; the
  worker now also leaves its reason beside the dataset, and the app reads it.
- **Error messages stay on screen until closed** instead of vanishing after
  a few seconds.

## [0.13.0] - 2026-10-21

Every dataset of SpatialData-DB on Explore, and a page of its own for every
dataset.

### Added

- **SpatialData-DB on Explore** (#162). Every dataset of the public
  [SpatialData-DB](https://lamin.ai/scverse/spatialdata-db) instance is an
  Explore card, read from the `catalog.parquet` the project publishes with
  each `catalog-<date>` release of `theislab/spatialdata-db-client` (their
  schema vendored in `lib/sddb/catalog-schema.ts`; one row is one card in
  `lib/sddb/catalog-row.ts`). `scripts/sync-sddb-catalog.ts` registers each
  store as an admin-owned dataset named from the catalog, upserts by the
  catalog uid, unpublishes uids that left the release and fails unless what
  it published equals the release's uid set; `--limit=N` syncs a sample
  spread over the technologies. Run by hand, dev first, then production;
  never on a schedule. A **SpatialData-DB** category, a licence chip on
  every card that has one, licence and source in the info modal and a
  "View on SpatialData-DB" link on the page. Catalog rows gain `source`
  (`curated`, `bil`, `community`, `sddb`) and the licence and provenance
  columns (migration `20261021100000_catalog_sddb_source`). Spec:
  `docs/specs/sddb-catalog.md`.

### Changed

- **Every dataset opens its own page** (#162). An Explore card is a plain
  link to the dataset's page (`/explore/[id]`, the BIL page for BIL rows)
  instead of jumping straight to a viewer; the page lists the dataset's
  parts and opens the viewer from there. The page groups entries by kind,
  so a SpatialData store is listed under Single Cell and opens in the
  SpatialData viewer, and an image gets its own row; before, the page
  showed no entry at all for a SpatialData store. Both links come from one
  module, `lib/explore/links.ts`.

## [0.12.0] - 2026-10-20

Molecules streamed by view from a tiled points element, cell polygons,
Visium and Visium HD, one way to prepare data from any source, and the
loading and rendering of every element measured against human time.

### Added

- **Molecules as a multiscale Zarr points element, streamed by view**
  (#150). A SpatialData points element stored as Zarr instead of Parquet
  (`docs/points-zarr-spec.md`): level 0 holds every molecule, each further
  level a uniform quarter; each level in two row orders (by tile, by gene)
  with an index, so the viewer reads only the tiles under the camera, at
  the level whose count fits the budget, by HTTP range requests. A
  **Detail** setting (Auto, More, All), a molecule panel with per-gene
  progress and the bytes fetched, and molecules kept at least 1 px when
  zoomed out. Measured on an 11.24 billion-molecule, ten-section mouse
  brain store: a whole-atlas view of three genes reads 29 MB.
- **Cell polygons from a store's shapes** (#151). A **Shapes** button on the
  rail: outline or filled, with an opacity slider, each polygon in its
  cell's colour, drawn when zoomed in below 1 mm.
- **Visium and Visium HD** (#158, #159). A store's elements are placed
  through its named coordinate systems, so Visium stores, which have no
  `global`, open with the H&E under the spots. Spots and bins are drawn at
  their real size and shape, as the cells, in the cell's colour, with
  nothing extra downloaded (the radius from the shapes file's footer, the
  bin side from the grid's pitch). The curated default table is the one
  opened. A census of every store in SpatialData-DB (1,884) and the
  profile of what they need (`docs/specs/sddb-profile.md`,
  `docs/plans/SDDB.md`).
- **Preparing data: one script per source into one SpatialData store,
  with the tutorial** (#160). `/docs/prepare`, the README and the script
  guide show exactly what you start with and what you get. Scripts:
  `prepare_folder.py` for parts you already have (an AnnData table in
  h5ad or Zarr, a transcripts table, an OME-Zarr image or a per-gene
  molecule folder), `transcripts_to_points.py` (a Parquet or CSV table
  straight into the tiled points element), `image_to_ome_zarr.py` (a TIFF
  as OME-Zarr with a pyramid), `prepare_xenium.py`, `prepare_merscope.py`
  and `prepare_visium.py`. The output is one store: the table with a CSC
  copy of `X` and, with `--cluster`, Leiden clusters in `obs/leiden`; the
  image written by the `spatialdata` library; the molecules as the tiled
  points element; the DEG stats in an `annotations/` sidecar. The expanded
  layout (store, molecule folder and image folder side by side) still
  opens and is what `--expanded` writes. A store the viewer cannot read
  as written is repaired in place: shapes parts into one WKB file, image
  chunks over 1024² rechunked, parquet points into the tiled element.
  Names in `mapping.json` may be relative to its folder.
- **Open in 3D by default** (#148): an Advanced option in the upload
  dialogs, and the owner's saved 2D/3D choice is now applied for
  single-cell datasets.
- **The formats as specifications** (#157): `docs/specs/` profiles for
  the store layout, the table, the shapes, the image and `mapping.json`,
  JSON Schemas, `scripts/validate_store.py` (run in CI on the fixtures),
  `scripts/legacy_to_points.py`, and a Python test job that also checks
  the `spatialdata` library reads what we write.
- **Budgets** (#152): `docs/specs/BUDGETS.md` states the human-time limits
  every change is judged by, with the perf suite measuring open-to-cells,
  first molecules, view settled, bytes per view and frames over 16 ms on
  the tiny and a generated medium store.

### Changed

- **One layer contract and one dataset loader** (#153). Every element
  drawn (image, shapes, points, molecules) is a `Layer` with one lifecycle
  in a registry; a scene renders only when something changed; every page
  and split panel opens a dataset through one loader with one loading
  screen.
- **One store for every reader** (#154). `ElementStore` holds one byte
  cache, shares reads in flight, cancels every read all the way to `fetch`,
  merges neighbouring range reads and bounds requests in flight; the
  consolidated metadata is seeded so opening a node costs no request; the
  points layer keeps a byte budget per view. On the medium store a gene
  settles in 2.7 s where it took 14.9 s, and a view in 13.5 s where it
  took 140 s.
- **Rendering within the frame budget** (#155). Molecule clouds carry
  positions only, with colour, size and shape as uniforms, so a restyle
  uploads nothing; cell buffers are reused and the gene percentile cached;
  hover picks through a grid; the pixel ratio is capped at 2; the shapes
  fill draws under the molecules. A heavy slider tick takes 2.6 s where it
  took 32 s.
- A cell dataset's overlays (a dropped folder's parts, the link registry,
  `mapping.json`, the app's mapping route) are resolved outside the scene
  (#156).
- The molecule overlay's size slider takes the same right-click min/max
  range editing as the cell slider (#149).
- Leiden clusters from `spatialdata_cluster.py` go into the table itself;
  the sidecar carries only the DEG stats (#160).
- The public skin reference store on the bucket was repaired with the new
  scripts: its polygons now load, and its image is in 1024² tiles.

### Fixed

- Dots invisible at the default size on graphics cards that draw nothing
  for sub-pixel points: the point size is clamped to at least 1 px (#149).
- A Visium image was drawn 148× too large (198× on Visium HD): an element's
  transform list holds alternatives, one per system, and was multiplied
  (#158).
- An H&E image showed as one blue channel: channels sharing a chunk are
  filled from one read, and r, g, b images open as RGB (#158).
- Recolouring 1.59 M cells by gene froze for seconds: the 95th percentile
  is taken by selection, not by sorting every cell (#150).
- A store whose index anndata writes as a nullable-string group registered
  with 0 cells and genes (#150).
- The tooltip said "undefined" with no column to colour by; a real-size
  cell keeps its size in every mode (#159).

## [0.11.0] - 2026-10-14

A user guide and a new About page.

### Added

- **A user guide at /docs, in the navbar.** Getting started, data formats
  (with a Python example for writing each one yourself), loading data, using
  the viewer, expression and DEGs (where values come from, raw vs log, and how
  DEGs are ranked), and sharing on Explore — with diagrams, and links that
  land on the exact section.
- **A new About page**: the team and principal investigators with photos,
  the Cleveland and Bintu Labs at UC San Diego, funding, and how to cite the
  bioRxiv preprint.

### Changed

- README, CLAUDE.md, TESTING.md and docs/spatialdata describe v0.10.0's
  features: Load from S3 for stores, images and folders; adding a store opened
  by URL to an account; image settings in the link; a Zarr v3 table on its
  own; images written by the spatialdata library.

## [0.10.0] - 2026-10-12

Data that is already processed opens as it is: drop a folder of cells,
molecules and an image, or paste the URL of one on a public bucket, and they
open overlaid with nothing converted.

### Added

- **Drop an OME-Zarr image to look at it.** An image folder dropped on the
  single-cell Folder card (Preview in browser) opens in the image viewer,
  read from disk as you zoom. Nothing is uploaded.
- **Already-processed cells, molecules and an image open together, as they
  are.** A folder holding a SpatialData store or a merfisheyes cell folder,
  a merfisheyes molecule folder and an OME-Zarr image — any of them, in any
  combination, recognised by content — opens overlaid without converting
  anything, so its size does not matter.
- **Upload & Save keeps the overlays.** Saving such a preview uploads every
  part as it is, each as its own dataset in one project, and the saved link
  opens with the molecules over the cells and the image under them. A
  SpatialData store can now be uploaded this way.
- **A Zarr v3 table opens on its own.** An AnnData written to Zarr v3, or
  the table copied out of a SpatialData store (`tables/<table>/`, or
  `tables/` itself), previews and saves without the store around it. Zarr v2
  AnnData folders open as before.
- **Load from S3 opens stores, images and whole folders.** The homepage box
  takes a SpatialData store, an OME-Zarr image, or the folder holding a store
  with its molecules beside it, which opens with everything overlaid. The
  folder's `mapping.json` names the store (`cells`); a store named after its
  folder is found without it.
- **A SpatialData store opened by URL can be added to your account**, like
  cells, molecules and images opened by URL.
- **Image settings are kept in the link.** Channels on or off, their colours
  and their lower and upper bounds go into the URL as they are changed, so a
  shared link looks the same.
- **The upload guide covers images and processed folders**, on the site and
  in the README.

### Fixed

- An image written by the `spatialdata` library is placed correctly: its
  scale-then-translation is stored as a "sequence" transform, which was
  ignored, drawing the image at one pixel per unit.
- Load from S3 now says why a URL could not be opened; the message was never
  shown.

### Changed

- A single file in a zarr-folder upload may now be up to 5 GB (was 200 MB):
  a sharded array is one file however large it is.

## [0.9.0] - 2026-10-11

Images join cells and molecules: an image is a dataset you can own, draw
under your cells, and upload in the same folder. SpatialData stores open
without conversion, and you can delete what you uploaded.

### Added

- **SpatialData stores open as they are.** `/spatialdata-viewer` reads a
  SpatialData zarr (v3) store straight from S3 — by dataset id or by the public
  URL of the store — without converting it: cells from its table, expression
  per gene on demand, and its image under the cells at the resolution the zoom
  needs. Molecules and Leiden clusters with DEGs come from sidecars written
  beside the store by two scripts. See `docs/spatialdata/README.md`.
- **A store may hold any of cells, molecules and an image.** Cells give the
  single cell scene with the others layered on it; molecules without cells
  give the single molecule scene; an image alone is shown on its own.
- **Tables written by current anndata are read.** anndata ≥ 0.12 shards its
  arrays and stores strings as nullable by default; both are now read (shards
  by byte range, so the bucket's CORS must allow the `Range` header).
- **Images are datasets.** An OME-Zarr image opened by URL at
  `/image-viewer/from-s3?url=…` can be added to your account like any other
  dataset opened by URL, and then appears under Your datasets. Each channel
  has a lower and an upper display bound. The viewer is not linked from the
  homepage yet.
- **An image under any single cell dataset.** A dataset's owner picks one of
  their images in the overlay dialog — e.g. a stitched DAPI mosaic — and it
  is drawn under the cells, alongside the molecule overlay. An image not in
  the account yet can be added by URL from the same dialog.
- **Cells, molecules and an image from one folder.** A folder uploaded for
  server processing may hold a transcripts file and an OME-Zarr image folder
  beside the cells. Each becomes its own dataset in one project, and the
  link you get opens the cells with the molecules over them and the image
  under them.
- **Delete your datasets.** From a dataset's card menu or its page, an owner
  can permanently delete an upload: its files, its row, and any job still
  processing it. A dataset shown on Explore has to be withdrawn first.
  Deleting molecules or an image unlinks them from the cells that overlaid
  them; deleting cells can take their molecule overlay along.
- **What can I upload?** A button on the homepage opens a guide to every
  accepted format — one tab each, with the expected file or folder structure
  and the rules that trip people up. The same guide replaces the empty
  `/docs` page.

### Fixed

- **Each side of a split view has its own molecule overlay.** With two cell
  datasets side by side, the right panel drew the left panel's molecule
  genes, had no legend, and its controls changed nothing. Each panel now
  reads its own overlay, the right one has its own slot in shared links, and
  the sync toggle decides whether they follow each other.
- **Large molecule CSVs upload for server processing.** On a file of a few
  hundred million rows the "Confirm columns" dialog could open with no
  columns in it, leaving its button disabled and the upload impossible to
  start. Server processing detects the columns itself, so it no longer asks;
  the in-browser dialog falls back to reading just the header.
- **A failed upload can be retried.** The row left behind by a failed or
  interrupted upload kept its fingerprint, so uploading the same file again
  was refused as a duplicate while the duplicate check said it was safe.
  Fingerprints are still computed and stored, but no longer block an upload.

## [0.8.0] - 2026-09-29

Explore ranks what people actually open, every card explains itself, and large
uploads no longer outlive their credentials.

### Added

- **Explore orders the grid by views.** It listed catalog rows in insertion
  order, so nothing surfaced what people open most. View counts now live on
  the catalog row and are kept current by the view endpoint, so the list and
  the API agree; admins keep `sortOrder` in front to pin a row, and the
  featured carousel stays hand-curated.
- **An info button on every Explore card.** Only BIL and project-backed cards
  could reach the detail page; a curated single-entry card went straight to
  the viewer, so its description, provenance and stats had nowhere to surface.
  The ⓘ opens a modal, so search, filters and scroll position survive a quick
  look, with a link through to the full page for the gene list.
- **`demo1`–`demo4.merfisheyes.com`** land on one dataset each, the way
  `spiralia.merfisheyes.com` does: the ACE mouse brain with a saved view, the
  spiralia embryo, and the schier and heart sites.

### Fixed

- **Large chunked uploads finish.** The upload modal signed one presigned URL
  per file with a one-hour expiry and uploaded files one at a time, so a big
  pre-chunked folder could outlive its URLs partway through — and the retries
  reused the expired ones. Files now upload eight at a time under a single
  POST policy that refreshes itself before it expires, and completion is
  recorded in batches instead of one request per file.

## [0.7.0] - 2026-09-17

Every embryo reachable in one click from the filmstrip, and a selection that
follows you between them.

### Added

- **Every embryo has its own tile in the filmstrip**, in developmental order,
  and clicking one opens it. The strip grouped 45 embryos into 11 stage tiles,
  so reaching a dataset it was already showing took two clicks — one to open
  the stage, another to pick from a list inside it. Each run of tiles still
  carries its stage as a heading, highlighted when the open embryo is in it, so
  nothing about where you are in the series was lost.
- **Selections carry across embryos.** Switching used to clear everything.
  Selections now follow, pruned against the embryo you are opening so nothing
  stays selected that it does not have — gene selections survive the whole
  series, while cell and domain selections survive within a developmental stage
  and fall away across one, which is what the labels support.
- **Splitting from the filmstrip turns sync on**, since opening a second panel
  from the strip is a request to compare the same selection on another embryo.

### Fixed

- **The filmstrip no longer leaves a bar at its right edge when scrolled.** The
  shared preview canvas was positioned inside the scrolling element, so it slid
  out of view with the content, uncovering the right edge and drawing every
  preview shifted by the scroll offset.
- **A split panel shows the same colours as the main one.** Syncing rebuilt
  each gene's colour from scratch in the receiving panel, so the same gene
  could appear in two different colours across the split.
- **The "made by" credit sits above the filmstrip** rather than on top of it.

## [0.6.0] - 2026-09-15

Nuclei alongside cells in the labelled-molecule viewer, and every featured
atlas reachable on Explore.

### Added

- **Nuclear surfaces in the labelled-molecule viewer**, for all 45 spiralia
  embryos — unlike the cell meshes, which exist for 31. Nuclei are keyed by the
  same cell labels, so they follow the cell menu exactly as the cell surfaces
  do, and they have their own on/off and opacity controls. Built from the
  segmentation volume by `scripts/spiralia/export_nuclei.py`, decimated to
  about 160 KB per embryo.

### Changed

- **The featured row on Explore is a carousel.** It rendered a 4-up grid and
  hid everything past the fourth card behind a "view all" link, so growing from
  3 featured datasets to 8 left most of them invisible. Arrows page a viewport
  at a time, disable at each end, and hide when the row already fits.
- **Labelled-molecule defaults follow a hand-tuned view** rather than first
  guesses: smaller selected molecules, a fainter cell surface, and nuclei shown
  by default. Nuclei are only on because the cell surface is now faint enough
  not to compete with them.
- **Surface toggles are switches**, and the camera panel's section headings now
  match the weight of the slider labels beside them — they were the dimmest
  text in the panel, where the eye needs an anchor to skim.

### Fixed

- **Admin-curated atlases get their detail page.** The link was gated on a
  dataset being a community submission, so a curated atlas that is
  project-backed in exactly the same way only ever offered an inline dropdown.
- **An external-source button names its destination.** It read "View on BIL"
  for every external link regardless of where it pointed.

## [0.5.1] - 2026-09-10

### Added

- **`spiralia.merfisheyes.com`** lands on one spiralia embryo with a set
  camera. Only the root is redirected, so every other route still resolves on
  that host — which is what lets the viewer's in-place dataset switching keep
  working once you are there.

## [0.5.0] - 2026-09-10

A new dataset type for MERFISH data where every molecule carries labels as
well as coordinates, and a browsable atlas of 45 spiralia embryos built on it.

### Added

- **Labelled single-molecule viewer** (`/lm-viewer`). A third dataset type
  alongside single-cell and single-molecule: every molecule carries a gene, an
  RNA domain and a cell, with three filter menus that intersect (OR within a
  menu, AND across menus). Any of the three can drive colour. Selection and
  sizing are applied on the GPU through lookup textures — one texel per
  category — so toggling a value re-uploads under a kilobyte instead of
  rebuilding multi-million-element attribute arrays.
- **Spiralia embryo atlas.** 45 embryos from the 1-cell stage to 24 cells,
  grouped into one project.
- **Stage rail.** An Allen-atlas style filmstrip along the bottom of the
  viewer, one tile per developmental stage, each a live 3D preview rather than
  a screenshot. All tiles share a single WebGL context drawn with per-tile
  scissor rectangles. Hovering opens a card listing that stage's embryos with
  the selected one previewed beside them; opening another embryo swaps it in
  place with no page reload.
- **Cell segmentation meshes** for the 31 embryos that have them, as a
  translucent or wireframe overlay.
- **Split screen for labelled molecules**, syncing selections and the colouring
  column between the two panels. Incoming values are intersected with the
  target dataset's own vocabulary, so gene selections carry across the whole
  series while cell and domain selections carry within a developmental stage
  and drop across one — which is what the data supports.
- **Quantised coordinates** (`q16`) for labelled-molecule datasets, cutting a
  typical embryo from 34.7 MB to 19.5 MB.
- **Dataset previews.** A `preview/` folder beside each dataset holding ~5,000
  sampled points at roughly 32 KB, enough to render a recognisable thumbnail
  without touching the full dataset.
- **Public project reads** (`/api/projects/[id]/public`), served only when
  every dataset in a project is admin-owned.

### Changed

- **Real byte progress** while a viewer downloads coordinates, so the bar
  tracks the download instead of jumping between fixed checkpoints.

### Fixed

- **Scenes redraw only when something changes.** The render loop drew every
  frame unconditionally, so a stationary camera still redrew every molecule 60
  times a second. Opt-in per viewer, currently used by the labelled-molecule
  viewer; the cell and single-molecule viewers are unchanged.
- **Camera defaults could not be saved** by an admin on an admin-owned
  labelled-molecule dataset — the permission check read a session field that
  does not exist, so it was always false.
- **A split panel no longer overwrites the main panel's dataset identity**,
  which had made "save as defaults" target the wrong dataset once a split was
  open.
- **Chunked datasets revalidate their manifest**, so replacing a dataset's
  files in place no longer serves a stale manifest from cache.

## [0.4.1] - 2026-09-04

Single-molecule loading becomes visible and cancellable, and pre-chunked
single-molecule uploads stop timing out on production.

### Fixed

- **Pre-chunked single-molecule uploads no longer 500 on production.** The
  initiate endpoint created one DB row and one presigned URL per file serially
  inside a single transaction (~33s for a 529-gene dataset — past Vercel's
  function timeout). It now batch-inserts the file rows and presigns outside
  the transaction in parallel batches (~2s).
- **Streamed default genes render again.** When every default gene streamed,
  the camera auto-fit ran against still-empty clouds and produced a broken
  (NaN) camera that never recovered, leaving the viewer permanently blank.
- **WebGL failures are no longer silent.** If the browser can't create a WebGL
  context (GPU blocklist, half-applied browser update, software GL), both
  viewers now show a "3D view unavailable" panel with remediation hints
  instead of a blank canvas.

### Changed

- **The single-molecule scene shows on first data.** The camera fits as soon
  as the first gene (or first streamed batch) has points, so molecules appear
  progressively instead of after every default gene finishes downloading.
- **Gene loads show progress and truly cancel.** Whole-file downloads report
  byte progress in the loading toast (e.g. `Drd1: 12.4 / 38.2 MB`); loading
  toasts are blue, and closing one aborts the download and deselects the gene.
  Deselecting a gene or leaving the view also aborts in-flight downloads.
- **Streaming kicks in at 100k molecules** (was 2M), so most genes stream in
  with visible progress.
- **Gene expression no longer fades low-expression points by default**
  (`EXPRESSION_ALPHA_MIN` 0.3 → 1.0; still adjustable in the advanced panel).

## [0.4.0] - 2026-09-01

Single-molecule overlays become first-class on single-cell datasets, the DEG
panel is rebuilt around a plain comparison of means, and signing in to upload
stops getting in the way.

### Added

- **Single-molecule overlays on single-cell datasets.** Link a molecule dataset
  to a cell dataset and its spots render on top of the cells. The overlay picker
  is searchable and paginated, works for both app-uploaded and S3-registered
  molecule datasets, and the selected molecule genes now persist in the viewer
  link (`ov=`) so a shared URL restores them.
- **Route DEG gene clicks to the overlay.** When a molecule overlay is present,
  a Cell / Molecule toggle beside the Expression control chooses whether clicking
  a gene colours the cells by expression or adds it as a molecule layer.
- **Expression-mode explainer.** An info button by the Expression control shows
  how counts-vs-log-normalized is auto-detected and how it affects the means and
  fold change.

### Changed

- **DEG panel shows mean1 / mean2 / linear fold change** instead of
  log2FC / mean_in / mean_pct — a direct comparison of means, with counts-vs-log
  auto-detected (and overridable) and log means un-logged before the ratio.
  Higher-contrast target / reference and mean columns.

### Fixed

- Signing in no longer blocks uploads: the email field is hidden when signed in
  and the account email is preserved, so "Process & Save" / "Process & Upload"
  enable correctly. Sign-in happens in place from the upload modals.
- Retroactive overlays no longer orphan a duplicate re-upload; the finalize step
  redirects to the existing dataset.
- Large single-cell `.h5ad` files with BigInt values no longer crash the loader,
  and the homepage featured datasets and server-fallback modal behave correctly.

## [0.3.0] - 2026-08-25

Large datasets stop being a wall: whole-transcriptome Xenium runs convert with
their full expression matrices, out-of-memory jobs retry themselves on bigger
machines, and big single-molecule genes stream into the viewer instead of
loading whole.

### Added

- **Xenium `cell_feature_matrix.h5` support.** Current Xenium exports ship
  expression only as the 10x HDF5 matrix, which neither pipeline could read —
  datasets loaded cells with zero genes. The server now reads it natively
  (a 717,576-cell x 27,104-gene run converts in 44 minutes), and the browser
  loads standard panels in full; whole-transcriptome matrices too large for a
  tab fall back to cells + gene list instead of crashing.
- **Automatic memory-tier escalation.** Processing jobs that run out of memory
  are retried on larger instances (16 -> 32 -> 64 GB) without user action, and
  admins can retry failures from the dashboard.
- **Progressive single-molecule streaming.** Gene files are written in
  bit-reversal order and streamed into the scene, so a half-loaded gene shows
  the whole slide densifying instead of one corner.
- **Account & Projects overhaul**: dataset cards, a project picker, optimistic
  toggles, and a shared dataset detail page.
- **Benchmark harness** (`bench/`, `scripts/bench/`): measures browser vs
  server processing across a real-dataset databank, with per-machine limits,
  costs, and a self-contained report page; 64-dataset results included.

### Changed

- The ingestion worker reads `.h5ad` selectively and streams the expression
  matrix: a 3-million-cell dataset that previously needed >64 GB now converts
  in a 16 GB instance (peak 13.3 GB).
- Browser and server pipelines produce byte-identical chunked output,
  enforced by a parity harness.
- The single-molecule gene picker paginates (250 genes per page) —
  thousand-gene panels no longer lag the panel open.

### Fixed

- Ingestion failures report the container's real reason instead of the
  generic task status.
- Community catalog entries resolve thumbnails from the live dataset.
- Unnamed single-molecule files upload correctly.

## [0.2.0] - 2026-08-19

Explore overhaul: the catalog is fast to search, gene search actually works
(and covers uploads), and datasets open identically everywhere.

### Added

- **Multi-gene search** on Explore: comma/space-separated terms are AND-matched
  with fuzzy per-term matching; cards show which genes matched; Enter creates
  one exact chip per term.
- **Gene search over uploaded datasets**: new `Dataset.genes` column filled
  automatically when processing completes (all three upload paths), a gene
  search box on "Your datasets", and gene copy-through on community
  submission. Backfill script for pre-existing rows
  (`scripts/backfill-dataset-genes.ts` — must be run once against production).
- **BIL → Projects organizer** (`scripts/ingest-bil-projects.ts`): groups the
  auto-registered BIL datasets into properly titled Projects and renames the
  "combined_output"-style rows (must be run once against production).

### Changed

- **Explore API is card-slim**: list payloads dropped from ~700 KB to
  ~17–41 KB (gene arrays never leave the database; featured/filters fetched
  once, not per keystroke). Search inputs are debounced with stale-response
  cancellation.
- **Dataset cards**: no more clipped titles/descriptions (fixed flexbox
  slicing), taller layout with investigator · age · publication-year byline,
  tooltips for full text, gradient fallback for broken thumbnails.

### Fixed

- **Split screen loads datasets like the full viewer**: one shared
  open-dataset pipeline (reset → owner config → URL-state overlay) for all
  viewer pages and both split panels — per-sample alignment, rotation, and
  saved colors now apply in the right panel, stores reset between datasets,
  and the loading screen clears only after coordinates are aligned (no more
  plot-then-jump). Same dataset in both panels renders identically.
- In-viewer catalog modal fetched 50 items under a "Showing 1–20" counter;
  now fetches exactly one page.

## [0.1.0] - 2026-08-18

First tracked release — a baseline snapshot of the production platform. Earlier
history is in the git log; changes from here on are recorded per release.

### Added

- **Single-cell viewer** (H5AD / Xenium / MERSCOPE): gene, celltype, and numerical
  coloring; combined gene + celltype mode; interactive legends and gene/numerical
  scalebar; a floating quantification plot panel (boxplots, histograms, celltype
  barplots) with secondary "group by", Top-N / axis caps, and CSV export.
- **Single-molecule viewer**: parquet/CSV, one point cloud per gene, S3 lazy
  loading, 2D/3D toggle, automatic control-probe / unassigned filtering.
- **Two ingestion paths**: in-browser (web workers) and **server-side** (AWS Batch
  worker) with email delivery. Server single-molecule uploads include a
  column-confirm/remap step; unassigned cell-ids are detected across MERSCOPE,
  Xenium, and CosMx conventions.
- **Combined single-cell + single-molecule upload**: a cells folder containing a
  transcripts file uploads both, auto-creates a project, and overlays molecules on
  cells.
- **Accounts** (Google sign-in), **Your datasets**, **Projects**, and an
  admin-reviewed **Community** catalog on Explore.
- **Split-screen** comparison of two datasets with synced camera and selection.
- **Admin processing dashboard** with failure classification and owner-notified
  failure emails.
- **Python preprocessing** (`process_spatial_data.py`, `process_single_molecule.py`)
  and BIL HPC / SLURM pipelines for very large datasets.

[Unreleased]: https://github.com/kresnajenie/merfisheyes/compare/v0.13.0...HEAD
[0.13.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.12.0...v0.13.0
[0.12.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.11.0...v0.12.0
[0.11.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.10.0...v0.11.0
[0.10.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/kresnajenie/merfisheyes/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/kresnajenie/merfisheyes/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/kresnajenie/merfisheyes/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/kresnajenie/merfisheyes/releases/tag/v0.1.0
