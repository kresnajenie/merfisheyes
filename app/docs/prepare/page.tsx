import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocSection,
  Flow,
  Prose,
  Steps,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Preparing data — MERFISHEYES Docs",
  description:
    "Cells, molecules and an image you processed yourself: one folder, one script, then drop it, save it or put it on a bucket.",
};

const pre =
  "overflow-x-auto rounded-xl border border-default-200 bg-default-50/40 p-3 text-xs";

export default function PreparePage() {
  return (
    <>
      <DocHeader
        heading="Preparing data"
        lead="If you process your own data, this is the one layout the viewer opens as it is: cells, molecules and an image in one folder, made ready by one script."
      />

      <DocSection heading="The layout" id="layout">
        <Flow
          steps={[
            {
              title: "Cells",
              detail:
                "An AnnData table (h5ad or Zarr) → a SpatialData store, Zarr v3, with a CSC copy of X",
            },
            {
              title: "+ Molecules",
              detail:
                "A merfisheyes molecule folder: manifest.json.gz + genes/*.bin.gz",
            },
            {
              title: "+ Image",
              detail: "An OME-Zarr image: axes c, y, x; 1024² chunks",
            },
            {
              title: "mapping.json",
              detail: "Names the three, relative to the folder",
            },
          ]}
        />
        <pre className={pre}>
          {`sample/
├── sample.zarr/        the store the script writes: tables/table (+ layers/X_csc)
├── molecules/          manifest.json.gz + genes/        (process_single_molecule.py)
├── dapi.zarr/          OME-Zarr image                   (any tool writing OME-Zarr 0.5)
├── annotations/        (optional) Leiden clusters + DEGs (--cluster)
└── mapping.json        { "cells": "sample.zarr", "links": { "__all__": "molecules" }, "image": "dapi.zarr" }`}
        </pre>
        <Prose>
          <p>
            Any part may be missing. Folder names don&apos;t matter: each part
            is recognised by what is inside it.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="What each part must be" id="parts">
        <CompareTable
          headers={["Part", "You have", "The viewer needs", "Who makes it so"]}
          rows={[
            [
              "Cells",
              'An AnnData table: .h5ad, or a Zarr folder from write_zarr (v2 or v3), with obsm["spatial"] in raw units and any obs columns',
              "A SpatialData store in Zarr v3; X readable one gene at a time (CSC), consolidated metadata",
              "The script: wraps the table, converts v2 and h5ad to v3, adds layers/X_csc when X is not CSC",
            ],
            [
              "Molecules",
              "A transcripts file (parquet / CSV) with gene, x, y (and z) columns",
              "manifest.json.gz + genes/<gene>.bin.gz, coordinates in the cells' units",
              "scripts/process_single_molecule.py (step 2)",
            ],
            [
              "Image",
              "A TIFF (one plane, a channel stack, or RGB), or an OME-Zarr folder",
              "OME-Zarr 0.5 (Zarr v3), axes c, y, x, uint8 or uint16, chunks ≤ 1024², scale + translation in the cells' units",
              "scripts/image_to_ome_zarr.py (step 2), or any OME-Zarr writer",
            ],
            [
              "Clusters & DEGs",
              "Perhaps nothing",
              "A categorical obs column, or an annotations/ sidecar",
              "The script with --cluster, when the table has no categorical column",
            ],
          ]}
        />
        <Callout kind="note">
          <b>anndata 0.12 writes Zarr v2 by default</b> (0.13 writes v3). A v2
          table opens from a drop, but a store and a bucket need v3. The script
          converts it, reading the table whole; for a large table set{" "}
          <code>anndata.settings.zarr_write_format = 3</code> before{" "}
          <code>write_zarr</code> and the script moves it into the store without
          rewriting it.
        </Callout>
      </DocSection>

      <DocSection heading="Run the script" id="run">
        <Steps
          steps={[
            {
              title: "Get the scripts",
              children: (
                <>
                  <pre className={pre}>
                    {`git clone https://github.com/kresnajenie/merfisheyes.git && cd merfisheyes
python -m venv .venv && source .venv/bin/activate     # Python 3.11 or newer
pip install -r scripts/requirements-prepare.txt`}
                  </pre>
                  <p>
                    One environment for every script on this page. Clustering
                    (scanpy, igraph, leidenalg) is in it too; leave those lines
                    out of the file if you never use <code>--cluster</code>.
                  </p>
                </>
              ),
            },
            {
              title: "Make the three parts",
              children: (
                <>
                  <pre className={pre}>
                    {`# cells: an AnnData table with obsm["spatial"] (x, y in microns), written by scanpy, anndata or Seurat's converters
#        cells.h5ad, or cells.zarr from adata.write_zarr("cells.zarr")

# molecules: a transcripts table → one file per gene
python scripts/process_single_molecule.py transcripts.parquet molecules/ --dataset-type xenium
#   --dataset-type merscope for gene/global_x/global_y/global_z columns;
#   any other names: --gene-col gene --x-col x --y-col y [--z-col z]

# image: a TIFF → OME-Zarr, placed among the cells by its pixel size and origin
python scripts/image_to_ome_zarr.py dapi.tif dapi.zarr --um-per-pixel 0.108 --origin 0 0 --channels DAPI`}
                  </pre>
                  <p>
                    Put them in one folder, side by side. Two candidates for one
                    part (two tables, say) make the script ask you to choose
                    with <code>--cells</code>, <code>--molecules</code> or{" "}
                    <code>--image</code>.
                  </p>
                </>
              ),
            },
            {
              title: "Run it",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_folder.py sample/              # add --cluster for Leiden clusters + DEGs`}
                  </pre>
                  <p>
                    It writes the store, <code>mapping.json</code> and, with{" "}
                    <code>--cluster</code>, <code>annotations/</code>, then runs
                    the validator and prints one line per check:
                  </p>
                  <pre className={pre}>
                    {`OK   root: consolidated metadata names all 26 nodes
OK   tables/table: X is CSR with a CSC copy in layers/X_csc
WARN tables/table: no instance_key: polygons cannot be matched to cells
OK   images/dapi.zarr: uint16 [1, 8192, 7941], chunks [1, 1024, 1024], 4 levels
OK   images/dapi.zarr: spans x 0–860, y 0–887, covering 100% of the cells' extent
WARN molecules/molecules: spans x 0–7941, y 0–8192, 9.2× the cells' extent: not the same units?`}
                  </pre>
                  <p>
                    A <code>FAIL</code> names what the viewer cannot read. A{" "}
                    <code>WARN</code> on a part&apos;s extent means it does not
                    sit where the cells are: the last line above is molecules
                    still in pixels while the cells are in microns. Fix the
                    source and run again; a second run changes nothing it has
                    already done.
                  </p>
                </>
              ),
            },
            {
              title: "Open it",
              children: (
                <p>
                  Drop the folder on the single cell <b>Folder</b> card in{" "}
                  <b>Preview in browser</b>: it opens at once with the molecules
                  over the cells and the image under them.{" "}
                  <b>Upload &amp; Save</b> keeps it in your account, every part
                  as it is.
                </p>
              ),
            },
            {
              title: "Or put it on a bucket",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_folder.py sample/ --url https://<bucket>.s3.<region>.amazonaws.com/<path>/sample
aws s3 sync sample/ s3://<bucket>/<path>/sample/`}
                  </pre>
                  <p>
                    <code>--url</code> is optional: names in{" "}
                    <code>mapping.json</code> are resolved against the folder
                    they are fetched from. The bucket must allow public reads
                    and this site in its CORS settings, with the{" "}
                    <code>Range</code> header allowed. Then paste the
                    folder&apos;s URL into <b>Load from S3</b>.
                  </p>
                </>
              ),
            },
          ]}
        />
        <Callout>
          <b>Same units, same frame.</b> The cells&apos;{" "}
          <code>obsm[&quot;spatial&quot;]</code>, the molecules&apos; x and y
          and the image&apos;s pixel size and origin must all be in one
          coordinate system — microns from the same corner is the usual one. The
          viewer draws them exactly where the numbers say; nothing is aligned
          for you.
        </Callout>
        <Callout>
          <b>A store you already have is repaired</b> where the viewer cannot
          read it as written: a shapes element stored as a directory of parquet
          parts, or in the geoarrow encoding, becomes one WKB file (the original
          stays beside it as <code>shapes.parquet.parts</code>); image levels in
          chunks over 1024² are rewritten in 1024² chunks; and{" "}
          <code>--molecules-from-points</code> turns the store&apos;s own points
          element into the molecule folder, since the viewer does not read
          parquet points.
        </Callout>
      </DocSection>

      <DocSection heading="Every combination opens" id="combinations">
        <CompareTable
          headers={["In the folder", "What you see"]}
          rows={[
            [
              "Cells + molecules + image",
              "Molecules over the cells, the image under both",
            ],
            ["Cells + molecules", "Molecules over the cells"],
            ["Cells + image", "The image under the cells"],
            [
              "Molecules + image",
              "The molecules as the scene, the image under them",
            ],
            ["Cells only", "The cells"],
            ["Molecules only", "The molecules, in the molecule viewer"],
            ["Image only", "The image, in the image viewer"],
            [
              "A SpatialData store holding several elements, no table",
              "Its molecules (a linked folder) with its image under them",
            ],
          ]}
        />
        <Callout>
          <b>Clusters and DEGs</b> stay a sidecar beside the store (
          <code>annotations/</code>): SpatialData has no place for
          differential-expression stats, and the viewer reads them from there
          for any store. <code>scripts/spatialdata_cluster.py</code> writes one
          for a table that already has cell types too.
        </Callout>
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
            reads the export with <code>spatialdata-io</code> and ends in the
            same checks as above. Xenium is below; MERSCOPE and Visium follow.
            Our own readers, <code>scripts/process_spatial_data.py</code> and{" "}
            <code>scripts/process_single_molecule.py</code>, still produce
            merfisheyes&apos; own folders, which open the same way.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Xenium" id="xenium">
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
              title: "Run it on the output folder",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_xenium.py Xenium_outs/ sample/ --workers 8   # add --cluster when the export has no cell types`}
                  </pre>
                  <p>
                    The folder with <code>experiment.xenium</code> in it. The
                    script reads the cell and nucleus boundaries, the
                    transcripts, the morphology image and the cell table, and
                    writes <code>sample/</code>:
                  </p>
                  <pre className={pre}>
                    {`sample/
├── Xenium.zarr/       tables/table (+ layers/X_csc), shapes/cell_boundaries + nucleus_boundaries,
│                      images/morphology_focus (1024² chunks, 6 levels), points/transcripts
├── molecules/         the transcripts as one file per gene, flat on the cells' plane (--keep-z keeps z)
├── annotations/       (--cluster) Leiden clusters + DEGs
└── mapping.json       { "cells": "Xenium.zarr", "links": { "__all__": "molecules" } }`}
                  </pre>
                  <p>
                    Coordinates are microns, as the export has them: the cells
                    and the molecules as written, the image placed through the
                    store&apos;s <code>global</code> system (pixels are 0.2125
                    µm). The raster cell labels are left out unless{" "}
                    <code>--labels</code>: the viewer does not draw them and
                    they are the largest part of the store; the table then
                    annotates the polygons, which carry the same cell ids.
                  </p>
                  <p>
                    Measured on a 130,870-cell, 7 M-transcript mouse brain
                    (Xenium 1.0): 70 s, 9 GB of memory at the peak, a 2.9 GB
                    store and 290 MB of molecule files; <code>--cluster</code>{" "}
                    another 47 s for 33 Leiden clusters. The report&apos;s{" "}
                    <code>WARN</code> lines on <code>points/transcripts</code>{" "}
                    (not read; the molecule folder is) and on the polygon files
                    (one row group) are expected.
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
                  image; the molecules stream per gene.
                </p>
              ),
            },
          ]}
        />
      </DocSection>

      <DocSection heading="MERSCOPE" id="merscope">
        <Prose>
          <p>
            A region folder as the instrument (or vpt) writes it:{" "}
            <code>cell_metadata.csv</code>, <code>cell_by_gene.csv</code>,{" "}
            <code>detected_transcripts.csv</code> (or <code>.parquet</code>)
            and, when you have them, <code>cell_boundaries.parquet</code> and{" "}
            <code>cell_categories.csv</code>. The store is built with the{" "}
            <code>spatialdata</code> library itself, not spatialdata-io, whose
            MERSCOPE reader stops without the <code>images/</code> folder that
            most shared exports lack.
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
              title: "Run it on the region folder",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_merscope.py region_R1/ sample/ --workers 8   # --z 3 picks the polygons' plane; --cluster when there are no cell types`}
                  </pre>
                  <pre className={pre}>
                    {`sample/
├── region_R1.zarr/    tables/table (counts without blanks, CSC; the metadata and categories as columns),
│                      shapes/cell_boundaries (each cell's polygon on one z-plane, --z)
├── molecules/         detected_transcripts as one file per gene, flat on the cells' plane (--keep-z keeps z)
├── annotations/       (--cluster) Leiden clusters + DEGs
└── mapping.json       { "cells": "region_R1.zarr", "links": { "__all__": "molecules" } }`}
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
                    transcripts CSV took 2 min 22 s and 15 GB of memory (65 MB
                    store, 530 MB of molecule files); a 703,879-cell region with
                    boundaries and categories took 3 min 19 s and 22 GB (550 MB
                    store with a 290 MB polygon file, 1.1 GB of molecules). The
                    counts CSV is read whole, which is where the memory goes.
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
          <p>
            A Space Ranger output folder: the counts h5, <code>spatial/</code>{" "}
            with the tissue positions, the scale factors and the hires and
            lowres images, and for Visium HD <code>binned_outputs/</code>. The
            spots, or the bins, are the cells: they are drawn at their real
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
              title: "Run it on the output folder",
              children: (
                <>
                  <pre className={pre}>
                    {`python scripts/prepare_visium.py outs/ sample/ --fullres image.tif   # the microscope image Space Ranger was given, if you have it
python scripts/prepare_visium.py outs/ sample/ --bins 8 16           # Visium HD: which bin sizes to read (8 and 16 µm by default)`}
                  </pre>
                  <pre className={pre}>
                    {`sample/
├── sample.zarr/       tables/table (+ layers/X_csc), shapes/sample (the spots, with their radius),
│                      images/sample_hires_image, sample_lowres_image (+ sample_full_image with --fullres)
└── mapping.json       { "cells": "sample.zarr" }`}
                  </pre>
                  <p>
                    Coordinates are full-resolution image pixels, as Space
                    Ranger has them. A Visium store has no shared coordinate
                    system: the spots and each image declare their own, and the
                    viewer targets the spots&apos; and draws the finest image
                    that reaches it, the full-resolution one when you gave it,
                    else the hires one. For Visium HD the 8 µm table is the one
                    the viewer opens; the 2 µm bins are too many for a browser
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
    </>
  );
}
