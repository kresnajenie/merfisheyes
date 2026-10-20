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
            process it yourself into the layout above, the per-platform scripts
            built on <code>spatialdata-io</code> follow in the next chapters;
            until then <code>scripts/process_spatial_data.py</code> and{" "}
            <code>scripts/process_single_molecule.py</code> produce
            merfisheyes&apos; own folders, which open the same way.
          </p>
        </Prose>
      </DocSection>
    </>
  );
}
