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
              "A transcripts file (parquet / CSV)",
              "manifest.json.gz + genes/<gene>.bin.gz",
              "scripts/process_single_molecule.py",
            ],
            [
              "Image",
              "A stitched plane, or an OME-Zarr folder",
              "OME-Zarr 0.5 (Zarr v3), axes c, y, x, uint8 or uint16, chunks ≤ 1024², scale + translation in the cells' units",
              "Any OME-Zarr writer; a 30-line example is in docs/spatialdata/README.md",
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
          <b>anndata writes Zarr v2 by default</b> (0.12). A v2 table opens from
          a drop, but a store and a bucket need v3. The script converts it; or
          set <code>anndata.settings.zarr_write_format = 3</code> before{" "}
          <code>write_zarr</code> and the script moves the table into the store
          without rewriting it.
        </Callout>
      </DocSection>

      <DocSection heading="Run the script" id="run">
        <Steps
          steps={[
            {
              title: "Put the parts in one folder",
              children: (
                <p>
                  The table (or a SpatialData store), the molecule folder and
                  the image folder, side by side. Two candidates for one part
                  (two tables, say) make the script ask you to choose with{" "}
                  <code>--cells</code>, <code>--molecules</code> or{" "}
                  <code>--image</code>.
                </p>
              ),
            },
            {
              title: "Run it",
              children: (
                <>
                  <pre className={pre}>
                    {`pip install -r tests/python/requirements.txt   # once: anndata, zarr, pyarrow, jsonschema
python scripts/prepare_folder.py sample/              # add --cluster for Leiden clusters (needs scanpy)`}
                  </pre>
                  <p>
                    It writes the store, <code>mapping.json</code> and, with{" "}
                    <code>--cluster</code>, <code>annotations/</code>, then runs
                    the validator and prints one line per check. A{" "}
                    <code>FAIL</code> names what the viewer could not read.
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
