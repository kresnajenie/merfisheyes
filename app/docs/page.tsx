import type { Metadata } from "next";

import {
  CompareTable,
  DocHeader,
  DocSection,
  Flow,
  LinkCards,
  LinkList,
  Prose,
} from "@/components/docs/doc-ui";

/** A tab of the format guide. */
const format = (key: string) => `/docs/formats?format=${key}#formats`;

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
              <LinkList
                key="cells"
                links={[
                  { label: ".h5ad", href: format("h5ad") },
                  { label: "Xenium", href: format("xenium") },
                  { label: "MERSCOPE", href: format("merscope") },
                  { label: "SpatialData", href: format("spatialdata") },
                  { label: "Zarr", href: format("zarr") },
                ]}
              />,
              <LinkList
                key="molecules"
                links={[
                  { label: "Parquet / CSV", href: format("molecules") },
                  {
                    label: "column names",
                    href: "/docs/formats#molecule-columns",
                  },
                ]}
              />,
              <LinkList
                key="image"
                links={[{ label: "OME-Zarr (.zarr)", href: format("image") }]}
              />,
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
          labelHeader="I have… / I want to…"
          rows={[
            [
              "Raw cells from a platform",
              <LinkList
                key="raw"
                links={[
                  { label: "H5AD", href: format("h5ad") },
                  { label: "Xenium", href: format("xenium") },
                  { label: "MERSCOPE", href: format("merscope") },
                  { label: "Preview it", href: "/docs/loading#preview" },
                  {
                    label: "or process it on our servers",
                    href: "/docs/loading#server",
                  },
                ]}
              />,
            ],
            [
              "A transcripts file",
              <LinkList
                key="molecules"
                links={[
                  { label: "Parquet / CSV", href: format("molecules") },
                  {
                    label: "column names",
                    href: "/docs/formats#molecule-columns",
                  },
                  {
                    label: "the single molecule viewer",
                    href: "/docs/viewer#single-molecule",
                  },
                ]}
              />,
            ],
            [
              "An image",
              <LinkList
                key="image"
                links={[
                  { label: "Image format", href: format("image") },
                  {
                    label: "drawing it under cells",
                    href: "/docs/viewer#overlays",
                  },
                ]}
              />,
            ],
            [
              "Cells, molecules and an image of one sample",
              <LinkList
                key="together"
                links={[
                  {
                    label: "Upload them in one go",
                    href: "/docs/loading#combined",
                  },
                  {
                    label: "overlays in the viewer",
                    href: "/docs/viewer#overlays",
                  },
                ]}
              />,
            ],
            [
              "Data that is already processed",
              <LinkList
                key="processed"
                links={[
                  {
                    label: "How it opens as it is",
                    href: "/docs/loading#processed",
                  },
                  { label: "processed folders", href: format("processed") },
                  { label: "SpatialData", href: format("spatialdata") },
                  { label: "Zarr table", href: format("zarr") },
                ]}
              />,
            ],
            [
              "A dataset too big for my browser",
              <LinkList
                key="big"
                links={[
                  {
                    label: "Upload & process on server",
                    href: "/docs/loading#server",
                  },
                  {
                    label: "or the Python scripts",
                    href: "/docs/loading#scripts",
                  },
                ]}
              />,
            ],
            [
              "Data on a public S3 bucket",
              <LinkList
                key="s3"
                links={[{ label: "Load from S3", href: "/docs/loading#s3" }]}
              />,
            ],
            [
              "Not sure how to load it",
              <LinkList
                key="choose"
                links={[
                  {
                    label: "Which way should I load my data?",
                    href: "/docs/loading#choose",
                  },
                  {
                    label: "the three ways side by side",
                    href: "/docs/loading#compare",
                  },
                ]}
              />,
            ],
            [
              "My own cell type labels",
              <LinkList
                key="labels"
                links={[
                  {
                    label: "Cell types you assigned yourself",
                    href: "/docs/formats#annotations",
                  },
                ]}
              />,
            ],
            [
              "Find my way around the viewer",
              <LinkList
                key="viewer"
                links={[
                  { label: "What is where", href: "/docs/viewer#layout" },
                  {
                    label: "colouring and plots",
                    href: "/docs/viewer#single-cell",
                  },
                  {
                    label: "two datasets side by side",
                    href: "/docs/viewer#split",
                  },
                ]}
              />,
            ],
            [
              "Understand the expression values and DEGs",
              <LinkList
                key="expression"
                links={[
                  {
                    label: "Where expression comes from",
                    href: "/docs/expression#source",
                  },
                  {
                    label: "raw counts or log",
                    href: "/docs/expression#raw-log",
                  },
                  {
                    label: "how DEGs are calculated",
                    href: "/docs/expression#deg",
                  },
                ]}
              />,
            ],
            [
              "Share what I am looking at",
              <LinkList
                key="link"
                links={[
                  { label: "Links keep your view", href: "/docs/viewer#links" },
                ]}
              />,
            ],
            [
              "Put my dataset on Explore",
              <LinkList
                key="explore"
                links={[
                  { label: "Step by step", href: "/docs/sharing#steps" },
                  {
                    label: "what happens after I submit",
                    href: "/docs/sharing#review",
                  },
                  { label: "projects", href: "/docs/sharing#projects" },
                ]}
              />,
            ],
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
              title: "Expression & DEGs",
              detail:
                "Where values come from, raw vs log, and how DEGs are ranked.",
              href: "/docs/expression",
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
