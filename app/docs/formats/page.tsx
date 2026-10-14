import type { Metadata } from "next";

import { Suspense } from "react";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocSection,
  Prose,
} from "@/components/docs/doc-ui";
import { FormatGuideFromLink } from "@/components/docs/format-guide-from-link";
import { FormatGuide } from "@/components/format-guide";

export const metadata: Metadata = {
  title: "Data formats — MERFISHEYES Docs",
  description:
    "Which files MERFISHEYES reads for cells, molecules and images, and how each is expected to look.",
};

export default function FormatsPage() {
  return (
    <>
      <DocHeader
        heading="Data formats"
        lead="MERFISHEYES reads the raw outputs of the common spatial platforms, so you don't need to pre-format anything. Drop a whole folder or a single file on the homepage; the format is worked out from what is inside."
      />

      <DocSection heading="At a glance" id="glance">
        <CompareTable
          headers={["What to drop", "Required", "Also read if present"]}
          labelHeader="Format"
          rows={[
            [
              "H5AD (AnnData)",
              "A single .h5ad file",
              "Coordinates in obsm['X_spatial']",
              "Embeddings such as obsm['X_umap'], and any cell type or metadata columns in obs",
            ],
            [
              "Xenium (10x)",
              "The Xenium output folder",
              "cells.csv(.gz) — the cell centroids",
              "An expression matrix (cell_feature_matrix.h5, the cell_feature_matrix/ folder, or cell_by_gene.csv); a transcripts file is read as molecules",
            ],
            [
              "MERSCOPE (Vizgen)",
              "The MERSCOPE output folder",
              "cell_metadata.csv — per-cell coordinates",
              "cell_by_gene.csv (expression matrix)",
            ],
            [
              "Molecules",
              "One .parquet or .csv of per-molecule rows",
              "A gene column, x and y (z for 3D data)",
              "A cell id column",
            ],
            [
              "Image",
              "An OME-Zarr image folder",
              "Zarr v3, axes c, y, x; uint8 or uint16",
              "Several resolution levels",
            ],
            [
              "Already processed",
              "A SpatialData store, an AnnData table in Zarr, or merfisheyes' own output folders",
              "See the tabs below",
              "Opens as it is — nothing is converted",
            ],
          ]}
        />
      </DocSection>

      <DocSection heading="Format by format" id="formats">
        <Prose>
          <p>
            Each tab shows the folder or file as it should look, which parts are
            required, and anything to watch out for. The same guide opens from{" "}
            <b>What can I upload?</b> on the homepage.
          </p>
        </Prose>
        <Suspense fallback={<FormatGuide />}>
          <FormatGuideFromLink />
        </Suspense>
      </DocSection>

      <DocSection heading="Molecule columns" id="molecule-columns">
        <Prose>
          <p>
            Column names are detected per platform, and you can remap them when
            you upload. Parquet is strongly preferred for millions of molecules:
            it is columnar and needs far less memory than CSV.
          </p>
        </Prose>
        <CompareTable
          headers={["Gene", "X", "Y", "Z", "Cell ID"]}
          labelHeader="Platform"
          rows={[
            [
              "Xenium",
              <code key="g">feature_name</code>,
              <code key="x">x_location</code>,
              <code key="y">y_location</code>,
              <code key="z">z_location</code>,
              <span key="c">
                <code>cell_id</code> (<code>UNASSIGNED</code> = no cell)
              </span>,
            ],
            [
              "MERSCOPE",
              <code key="g">gene</code>,
              <code key="x">global_x</code>,
              <code key="y">global_y</code>,
              <code key="z">global_z</code>,
              <span key="c">
                <code>cell_id</code> (<code>-1</code> = no cell)
              </span>,
            ],
            ["Custom", "You pick each column from a dropdown", "", "", "", ""],
          ]}
        />
      </DocSection>

      <DocSection heading="Cell types you assigned yourself" id="annotations">
        <Prose>
          <p>
            Cell types don&apos;t have to be in the file. Run{" "}
            <a
              href="https://knowledge.brain-map.org/mapmycells/process"
              rel="noreferrer"
              target="_blank"
            >
              MapMyCells
            </a>{" "}
            (or any tool) yourself, then choose{" "}
            <b>Upload &amp; process on server</b> and add the resulting
            one-row-per-cell CSV in the{" "}
            <b>MapMyCells celltype CSV (optional)</b> box on the homepage. Its
            columns become columns to colour by, with their own palettes and
            differential expression. There is no cell typing on our side, so the
            labels stay yours.
          </p>
        </Prose>
      </DocSection>

      <Callout>
        An image is drawn under cells only when the cells kept their original
        coordinates (usually microns). Folders written by older versions of the
        Python scripts normalized them; they still load, but an image can&apos;t
        be placed under them.
      </Callout>
    </>
  );
}
