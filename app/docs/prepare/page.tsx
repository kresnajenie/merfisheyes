import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocSection,
  Prose,
  Steps,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Preparing data — MERFISHEYES Docs",
  description:
    "Cells, molecules and an image you processed yourself, or a Xenium, MERSCOPE or Visium export: one script makes one SpatialData store the viewer opens as it is.",
};

const pre =
  "overflow-x-auto rounded-xl border border-default-200 bg-default-50/40 p-3 text-xs";

export default function PreparePage() {
  return (
    <>
      <DocHeader
        heading="Preparing data"
        lead="If you process your own data, this is the one layout the viewer opens as it is: a SpatialData store holding the cells, the image and the molecules, made by one script. Each section below shows exactly what you start with and exactly what you get."
      />

      <DocSection heading="What you get, whatever you start from" id="layout">
        <Prose>
          <p>
            One folder, with one SpatialData store in it. The store is written
            by the <code>spatialdata</code> library where it can be, so the same
            store opens in spatialdata, napari and Vitessce; the molecules are
            the one part in our own format, a tiled points element the viewer
            streams by view.
          </p>
        </Prose>
        <pre className={pre}>
          {`sample/
├── sample.zarr/                   one SpatialData store (Zarr v3)
│   ├── zarr.json                  spatialdata_attrs + consolidated metadata, written last
│   ├── tables/table/              the cells: AnnData — obs (your columns, + obs/leiden with --cluster),
│   │   ├── X                      obsm/spatial, var; X as written, plus
│   │   └── layers/X_csc           a CSC copy when X is not CSC (one gene = one slice)
│   ├── images/dapi/               the image: OME-Zarr, axes c, y, x, 1024² chunks, a pyramid,
│   │                              placed among the cells by its scale and translation
│   ├── points/molecules/          the molecules, straight from the transcripts table: a tiled
│   │                              points element (docs/points-zarr-spec.md) streamed by view
│   └── shapes/cell_boundaries/    the polygons, when you have them: one GeoParquet file
├── annotations/                   (--cluster) de/leiden.bin.gz: the DEG stats, which SpatialData
│                                  has no place for; the clusters themselves are in the table
├── sources/                       the image folder and the transcripts table the store was made from
└── mapping.json                   { "cells": "sample.zarr" }`}
        </pre>
        <Callout>
          <b>The expanded layout works too.</b> With <code>--expanded</code> the
          image and the molecule folder stay beside the store and{" "}
          <code>mapping.json</code> names them: <code>sample.zarr/</code>,{" "}
          <code>molecules/</code>, <code>dapi.zarr/</code> side by side. A
          folder in that form, or any subset of it, opens exactly the same way.
        </Callout>
      </DocSection>

      <DocSection heading="From parts you already have" id="parts">
        <Prose>
          <p>
            You start with a table of cells, a table of transcripts and a
            microscope image, in one coordinate system (microns from the same
            corner is the usual one):
          </p>
        </Prose>
        <pre className={pre}>
          {`you have/
├── cells.h5ad                 an AnnData table: obsm["spatial"] (x, y), X, obs columns
│                              (or cells.zarr from adata.write_zarr, Zarr v2 or v3)
├── transcripts.parquet        one row per molecule: gene, x, y (, z)   (or .csv)
└── dapi.tif                   one plane, a (c, y, x) stack or RGB; uint8 or uint16`}
        </pre>
        <Steps
          steps={[
            {
              title: "Get the scripts",
              children: (
                <>
                  <pre className={pre}>
                    {`git clone https://github.com/kresnajenie/merfisheyes.git && cd merfisheyes
python -m venv .venv && source .venv/bin/activate     # Python 3.12 or newer
pip install -r scripts/requirements-prepare.txt spatialdata`}
                  </pre>
                  <p>
                    One environment for every script on this page. Clustering
                    (scanpy, igraph, leidenalg) is in the file too; leave those
                    lines out if you never use <code>--cluster</code>.
                  </p>
                </>
              ),
            },
            {
              title: "The cells and the molecules: as they are",
              children: (
                <>
                  <pre className={pre}>
                    {`mkdir sample && cp cells.h5ad transcripts.parquet sample/`}
                  </pre>
                  <p>
                    Nothing to convert. The transcripts table is read straight
                    into the store&apos;s tiled points element. Xenium&apos;s
                    columns (<code>feature_name</code>, <code>x_location</code>,{" "}
                    <code>y_location</code>) and MERSCOPE&apos;s (
                    <code>gene</code>, <code>global_x</code>,{" "}
                    <code>global_y</code>) are recognised; any other names are
                    given to the script in the next step with{" "}
                    <code>--molecule-columns GENE,X,Y</code>. Molecules are
                    written flat on the cells&apos; plane.
                  </p>
                </>
              ),
            },
            {
              title: "The image: a TIFF becomes OME-Zarr",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/image_to_ome_zarr.py dapi.tif sample/dapi.zarr --um-per-pixel 0.108 --origin 0 0 --channels DAPI`}
                  </pre>
                  <p>
                    The pixel size and the origin (where the centre of pixel 0,
                    0 sits among the cells, in their units) are what place it
                    under them. Now the folder is:
                  </p>
                  <pre className={pre}>
                    {`sample/
├── cells.h5ad
├── transcripts.parquet
└── dapi.zarr/                 OME-Zarr: zarr.json + 0/, 1/, … (the pyramid)`}
                  </pre>
                </>
              ),
            },
            {
              title: "Run the script",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_folder.py sample/              # add --cluster for Leiden clusters + DEGs`}
                  </pre>
                  <p>
                    It finds the three parts by what is in them (names do not
                    matter), wraps the table into the store, adds the CSC copy,
                    writes the image into the store with the{" "}
                    <code>spatialdata</code> library and the transcripts table
                    straight into the tiled points element (the script that does
                    that, <code>transcripts_to_points.py</code>, also runs on
                    its own), moves the image folder and the table into{" "}
                    <code>sources/</code>, writes <code>mapping.json</code> and
                    runs the validator, one line per check:
                  </p>
                  <pre className={pre}>
                    {`OK   root: consolidated metadata names all 220 nodes
OK   tables/table: obsm/spatial [500, 2]
OK   tables/table: X is CSR with a CSC copy in layers/X_csc
OK   points/molecules: metadata valid (points:version 0.1, 8 levels)
OK   points/molecules: levels 0: 20000, 1: 5041, 2: 1219, 3: 278, …
OK   images/dapi: uint16 [1, 800, 1000], chunks [1, 800, 1000], 1 levels
OK   coordinates: target 'global' (global)
OK   mapping.json: names the store only`}
                  </pre>
                  <p>
                    The result is the layout at the top of this page. A{" "}
                    <code>FAIL</code> names what the viewer cannot read. With{" "}
                    <code>--expanded</code> the parts stay beside the store and
                    the script also checks that each sits where the cells are: a
                    part in other units or another frame gets a{" "}
                    <code>WARN</code> naming both extents.
                  </p>
                </>
              ),
            },
            {
              title: "Open it",
              children: (
                <p>
                  Drop <code>sample/</code> on the single cell <b>Folder</b>{" "}
                  card in <b>Preview in browser</b>: the cells open with the
                  image under them and the molecules streamed over them, gene by
                  gene. <b>Upload &amp; Save</b> keeps it in your account.
                </p>
              ),
            },
            {
              title: "Or put it on a bucket",
              children: (
                <>
                  <pre className={pre}>
                    {`aws s3 sync sample/ s3://<bucket>/<path>/sample/ --exclude "sources/*"`}
                  </pre>
                  <p>
                    The bucket must allow public reads and this site in its CORS
                    settings, with the <code>Range</code> header allowed. Then
                    paste the folder&apos;s URL into <b>Load from S3</b>. In the
                    expanded layout, <code>mapping.json</code> names the parts
                    by name and the viewer resolves them against the folder;{" "}
                    <code>--url</code> writes absolute links instead.
                  </p>
                </>
              ),
            },
          ]}
        />
        <Callout kind="note">
          <b>anndata 0.12 writes Zarr v2 by default</b> (0.13 writes v3). A v2
          table opens from a drop, but a store needs v3. The script converts it,
          reading the table whole; for a large table set{" "}
          <code>anndata.settings.zarr_write_format = 3</code> before{" "}
          <code>write_zarr</code> and the script moves it into the store without
          rewriting it.
        </Callout>
        <Callout>
          <b>Clusters and DEGs.</b> With <code>--cluster</code>, when the table
          has no categorical column, Leiden clusters are written into the table
          as <code>obs/leiden</code>, so every reader of the store sees them;
          only the per-cluster differential-expression stats go beside the store
          in <code>annotations/</code>, because SpatialData has no place for
          them. <code>scripts/spatialdata_cluster.py</code> does the same for a
          store you already have.
        </Callout>
      </DocSection>

      <DocSection heading="Every combination opens" id="combinations">
        <CompareTable
          headers={["In the folder", "What you see"]}
          rows={[
            [
              "A store with a table, an image and points",
              "Molecules over the cells, the image under both",
            ],
            ["A store with a table and points", "Molecules over the cells"],
            ["A store with a table and an image", "The image under the cells"],
            [
              "A store with an image and points, no table",
              "The molecules as the scene, the image under them",
            ],
            ["A table alone (a store, or a bare AnnData zarr)", "The cells"],
            [
              "A molecule folder alone",
              "The molecules, in the molecule viewer",
            ],
            ["An OME-Zarr image alone", "The image, in the image viewer"],
            [
              "Expanded: a store, a molecule folder, an image folder side by side",
              "The same as the wrapped forms above",
            ],
          ]}
        />
      </DocSection>

      <DocSection
        heading="From a platform export (Xenium, MERSCOPE, Visium)"
        id="platforms"
      >
        <Prose>
          <p>
            A raw export needs no preparation: switch the homepage to{" "}
            <b>Upload &amp; process on server</b> and drop the folder. To
            process it yourself into the layout above, one script per platform
            reads the export and ends in the same checks. Our own readers,{" "}
            <code>scripts/process_spatial_data.py</code> and{" "}
            <code>scripts/process_single_molecule.py</code>, still produce
            merfisheyes&apos; own folders, which open the same way.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Xenium" id="xenium">
        <Prose>
          <p>You start with the Xenium Onboard Analysis output folder:</p>
        </Prose>
        <pre className={pre}>
          {`Xenium_outs/
├── experiment.xenium              the run's specs (pixel size, version)
├── cell_feature_matrix.h5         counts
├── cells.parquet (or cells.csv.gz) cell metadata
├── cell_boundaries.parquet        polygons
├── nucleus_boundaries.parquet     polygons
├── transcripts.parquet            the molecules
├── morphology_focus.ome.tif       the image (or morphology_focus/ in newer versions)
└── cells.zarr.zip                 (1.3+) label ids`}
        </pre>
        <Steps
          steps={[
            {
              title: "Install the reader",
              children: (
                <pre className={pre}>
                  {`pip install -r scripts/requirements-prepare.txt spatialdata-io    # Python 3.12 or newer`}
                </pre>
              ),
            },
            {
              title: "Run it",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_xenium.py Xenium_outs/ sample/ --workers 8   # add --cluster when the export has no cell types`}
                  </pre>
                  <p>You get:</p>
                  <pre className={pre}>
                    {`sample/
├── Xenium.zarr/
│   ├── tables/table/              counts (+ layers/X_csc), cell metadata, obsm/spatial in microns
│   ├── shapes/cell_boundaries/    one GeoParquet each
│   ├── shapes/nucleus_boundaries/
│   ├── images/morphology_focus/   1024² chunks, 6 levels
│   ├── points/molecules/          transcripts.parquet, tiled, flat on the cells' plane
│   └── points/transcripts/        the reader's own parquet copy; not read by the viewer
├── annotations/                   (--cluster) the DEG stats
└── mapping.json                   { "cells": "Xenium.zarr" }`}
                  </pre>
                  <p>
                    Coordinates are microns, as the export has them; the image
                    is placed through the store&apos;s <code>global</code>{" "}
                    system (pixels are 0.2125 µm). The raster cell labels are
                    left out unless <code>--labels</code>: the viewer does not
                    draw them and they are the largest part of the store; the
                    table then annotates the polygons, which carry the same cell
                    ids.
                  </p>
                  <p>
                    Measured on a 130,870-cell, 42 M-transcript mouse brain
                    (Xenium 1.0): 73 s and 9 GB of memory at the peak, a 3.4 GB
                    store; <code>--cluster</code> another 47 s for 33 Leiden
                    clusters. The report&apos;s <code>WARN</code> lines on{" "}
                    <code>points/transcripts</code> (not read) and on the
                    polygon files (one row group) are expected.
                  </p>
                </>
              ),
            },
            {
              title: "Open it",
              children: (
                <p>
                  As above: drop <code>sample/</code> on the single cell{" "}
                  <b>Folder</b> card, <b>Upload &amp; Save</b>, or sync it to a
                  bucket. Zoomed in below 1 mm the cell polygons draw over the
                  image; the molecules stream by view.
                </p>
              ),
            },
          ]}
        />
      </DocSection>

      <DocSection heading="MERSCOPE" id="merscope">
        <Prose>
          <p>
            You start with a region folder as the instrument (or vpt) writes it:
          </p>
        </Prose>
        <pre className={pre}>
          {`region_R1/
├── cell_metadata.csv              one row per cell: EntityID, center_x, center_y, …
├── cell_by_gene.csv               counts, one row per cell, one column per gene
├── detected_transcripts.csv       gene, global_x, global_y, global_z, …   (or .parquet)
├── cell_boundaries.parquet        (optional) vpt polygons: EntityID, ZIndex, Geometry
└── cell_categories.csv            (optional) cell types / clusters per cell`}
        </pre>
        <Prose>
          <p>
            The store is built with the <code>spatialdata</code> library itself,
            not spatialdata-io, whose MERSCOPE reader stops without the{" "}
            <code>images/</code> folder that most shared exports lack.
          </p>
        </Prose>
        <Steps
          steps={[
            {
              title: "Install the library",
              children: (
                <pre className={pre}>
                  {`pip install -r scripts/requirements-prepare.txt spatialdata    # Python 3.12 or newer`}
                </pre>
              ),
            },
            {
              title: "Run it",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_merscope.py region_R1/ sample/ --workers 8   # --z 3 picks the polygons' plane; --cluster when there are no cell types`}
                  </pre>
                  <p>You get:</p>
                  <pre className={pre}>
                    {`sample/
├── region_R1.zarr/
│   ├── tables/table/              counts without blanks, CSC; the metadata and categories as columns
│   ├── shapes/cell_boundaries/    each cell's polygon on one z-plane (--z), one GeoParquet
│   └── points/molecules/          detected_transcripts, tiled, flat on the cells' plane
├── annotations/                   (--cluster) the DEG stats
└── mapping.json                   { "cells": "region_R1.zarr" }`}
                  </pre>
                  <p>
                    Everything is in microns, as the files have it. The counts
                    are aligned to the metadata&apos;s cell ids, so the two CSVs
                    may list cells in different orders. The boundaries file
                    holds every cell once per z-plane; one plane is kept (the
                    largest part of a MultiPolygon), and a plane that is not in
                    the file is refused with the planes that are. The mosaic
                    images are not written yet: a MERSCOPE mosaic is tens of
                    gigabytes, more than the image script holds in memory.
                  </p>
                  <p>
                    Measured: a 137,627-cell, 960-gene region with a 6 GB
                    transcripts CSV (70 M molecules) took about 3½ min and 16 GB
                    of memory, clustering included; a 703,879-cell region with
                    boundaries and categories about 4 min and 22 GB. The counts
                    CSV is read whole, which is where the memory goes.
                  </p>
                </>
              ),
            },
            {
              title: "Open it",
              children: (
                <p>
                  As above: drop <code>sample/</code> on the single cell{" "}
                  <b>Folder</b> card, <b>Upload &amp; Save</b>, or sync it to a
                  bucket.
                </p>
              ),
            },
          ]}
        />
      </DocSection>

      <DocSection heading="Visium and Visium HD" id="visium">
        <Prose>
          <p>You start with the Space Ranger output folder:</p>
        </Prose>
        <pre className={pre}>
          {`outs/
├── filtered_feature_bc_matrix.h5  counts
├── spatial/
│   ├── tissue_positions.csv       (or tissue_positions_list.csv) spot centres in full-res pixels
│   ├── scalefactors_json.json     spot diameter, hires / lowres scale
│   ├── tissue_hires_image.png
│   └── tissue_lowres_image.png
├── binned_outputs/                Visium HD only: square_002um/, square_008um/, square_016um/
└── (the microscope image Space Ranger was given, if you kept it)`}
        </pre>
        <Prose>
          <p>
            The spots, or the bins, are the cells: they are drawn at their real
            size, there are no molecules, and the H&amp;E lies under them.
          </p>
        </Prose>
        <Steps
          steps={[
            {
              title: "Install the reader",
              children: (
                <pre className={pre}>
                  {`pip install -r scripts/requirements-prepare.txt spatialdata-io    # Python 3.12 or newer`}
                </pre>
              ),
            },
            {
              title: "Run it",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_visium.py outs/ sample/ --fullres image.tif   # the microscope image, if you have it
python scripts/prepare_visium.py outs/ sample/ --bins 8 16           # Visium HD: which bin sizes to read (8 and 16 µm by default)`}
                  </pre>
                  <p>You get:</p>
                  <pre className={pre}>
                    {`sample/
├── sample.zarr/
│   ├── tables/table/              counts (+ layers/X_csc), obsm/spatial in full-res pixels
│   ├── shapes/sample/             the spots, with their radius
│   ├── images/sample_hires_image/
│   ├── images/sample_lowres_image/
│   └── images/sample_full_image/  with --fullres: 1024² chunks, a pyramid
└── mapping.json                   { "cells": "sample.zarr" }

# Visium HD: tables/square_008um, square_016um and shapes/sample_square_008um, …;
# the 8 µm table is the default (metadata_general["Default Table"])`}
                  </pre>
                  <p>
                    Coordinates are full-resolution image pixels, as Space
                    Ranger has them. A Visium store has no shared coordinate
                    system: the spots and each image declare their own, and the
                    viewer targets the spots&apos; and draws the finest image
                    that reaches it, the full-resolution one when you gave it,
                    else the hires one. The 2 µm bins are too many for a browser
                    and are not read unless you ask.
                  </p>
                </>
              ),
            },
            {
              title: "Open it",
              children: (
                <p>
                  As above: drop <code>sample/</code> on the single cell{" "}
                  <b>Folder</b> card, <b>Upload &amp; Save</b>, or sync it to a
                  bucket.
                </p>
              ),
            },
          ]}
        />
      </DocSection>

      <DocSection heading="A store you already have" id="repairs">
        <Prose>
          <p>
            Put the store in a folder and run <code>prepare_folder.py</code> on
            it. Where the viewer cannot read the store as written, the script
            repairs it in place: a shapes element stored as a directory of
            parquet parts, or in the geoarrow encoding, becomes one WKB file
            (the original stays beside it as <code>shapes.parquet.parts</code>);
            image levels in chunks over 1024² are rewritten in 1024² chunks; and{" "}
            <code>--molecules-from-points</code> turns the store&apos;s own
            parquet points element into the tiled one the viewer streams. The
            CSC copy of <code>X</code> is added as for any table.
          </p>
        </Prose>
      </DocSection>
    </>
  );
}
