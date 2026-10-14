import type { Metadata } from "next";

import {
  CompareTable,
  DocHeader,
  DocSection,
  Flow,
  LinkCards,
  Prose,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Docs — MERFISHEYES",
  description:
    "How to bring spatial transcriptomics data to MERFISHEYES, look at it, and share it.",
};

export default function DocsPage() {
  return (
    <>
      <DocHeader
        heading="Getting started"
        lead="MERFISHEYES shows spatial transcriptomics data in 3D in your browser — cells, the molecules inside them, and the image underneath — and lets you share it."
      />

      <DocSection heading="From your files to a shared link" id="overview">
        <Flow
          steps={[
            {
              title: "Bring your data",
              detail:
                "Raw platform output (H5AD, Xenium, MERSCOPE, a transcripts file), an image, or data that is already processed.",
              href: "/docs/formats",
            },
            {
              title: "Load it",
              detail:
                "Preview it in your browser, upload it to be processed on our servers, or open it from a public bucket.",
              href: "/docs/loading",
            },
            {
              title: "Explore it",
              detail:
                "Colour by gene or cell type, compare groups in plots, overlay molecules and an image, view two datasets side by side.",
              href: "/docs/viewer",
            },
            {
              title: "Share it",
              detail:
                "Keep it in your account, group it into projects, and submit it to the public Explore page.",
              href: "/docs/sharing",
            },
          ]}
        />
      </DocSection>

      <DocSection heading="Three kinds of data" id="kinds">
        <Prose>
          <p>
            Everything you load is one of three kinds. They are separate
            datasets, but they are shown together when they come from the same
            sample: molecules over the cells, the image under both.
          </p>
        </Prose>
        <CompareTable
          headers={["Cells", "Molecules", "Image"]}
          rows={[
            [
              "What it is",
              "One point per cell, with a gene expression matrix and cell metadata (cell types, clusters).",
              "One point per detected transcript, labelled with its gene.",
              "A microscopy image, such as a stitched DAPI mosaic.",
            ],
            [
              "Typical files",
              ".h5ad, a Xenium or MERSCOPE folder, a SpatialData store",
              "transcripts .parquet or .csv",
              "An OME-Zarr folder (.zarr)",
            ],
            [
              "On the homepage",
              "Single cell",
              "Single molecule",
              "Single cell → Folder",
            ],
          ]}
        />
      </DocSection>

      <DocSection heading="Where to go next" id="next">
        <CompareTable
          headers={["Go to"]}
          labelHeader="I have…"
          rows={[
            [
              "Raw platform output",
              "Data formats — exactly which files each platform needs",
            ],
            ["An image", "Data formats → Image"],
            [
              "Data that is already processed",
              "Loading data — it opens as it is, whatever its size",
            ],
            [
              "A dataset too big for your browser",
              "Loading data → Upload & process on server",
            ],
            ["Data on a public S3 bucket", "Loading data → Load from S3"],
            ["Putting your dataset on Explore", "Sharing on Explore"],
          ]}
        />
        <LinkCards
          links={[
            {
              title: "Data formats",
              detail: "What you can bring, file by file.",
              href: "/docs/formats",
            },
            {
              title: "Loading data",
              detail: "Which way to load it, and what happens to it.",
              href: "/docs/loading",
            },
            {
              title: "Using the viewer",
              detail: "Colouring, plots, overlays, images and links.",
              href: "/docs/viewer",
            },
            {
              title: "Sharing on Explore",
              detail: "Your account, projects, and getting reviewed.",
              href: "/docs/sharing",
            },
          ]}
        />
      </DocSection>
    </>
  );
}
