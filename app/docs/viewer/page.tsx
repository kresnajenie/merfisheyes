import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocLink,
  DocSection,
  Prose,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Using the viewer — MERFISHEYES Docs",
  description:
    "Colouring by gene or cell type, plots, molecule and image overlays, split screen, and links that keep your view.",
};

/**
 * The controls in the screenshot below, numbered where they sit. Positions
 * are pixels of the 1440 × 900 capture, measured from the page itself.
 */
const CONTROLS: { x: number; y: number; name: string; what: string }[] = [
  {
    x: 72,
    y: 68,
    name: "Image",
    what: "Show or hide the image, and for each channel whether it is on, its colour, and its lower and upper bound.",
  },
  {
    x: 72,
    y: 112,
    name: "Celltype",
    what: "Colour by a cell type, cluster or any other column; tick types to highlight them.",
  },
  {
    x: 72,
    y: 176,
    name: "Gene",
    what: "Colour the cells by one gene's expression. With molecules shown, it also picks molecule genes.",
  },
  {
    x: 72,
    y: 240,
    name: "DEG",
    what: "The genes that set one cell type apart from the rest, or from another type.",
  },
  {
    x: 72,
    y: 304,
    name: "Plot",
    what: "Boxplots, histograms and cell counts in a panel you can move and resize.",
  },
  {
    x: 72,
    y: 368,
    name: "Split screen",
    what: "Open a second dataset beside this one.",
  },
  { x: 30, y: 440, name: "Point size", what: "How big the dots are drawn." },
  {
    x: 103,
    y: 433,
    name: "Cells / molecules",
    what: "Whether the size slider and the gene picker act on the cells or the molecules. Shown when molecules are overlaid.",
  },
  {
    x: 72,
    y: 600,
    name: "SC · SM",
    what: "Show or hide the cells (SC) and the molecules (SM).",
  },
  {
    x: 1176,
    y: 72,
    name: "Camera controls",
    what: "2D or 3D, rotation, flips, reset the view, and a screenshot or thumbnail.",
  },
  {
    x: 1240,
    y: 72,
    name: "Hide UI",
    what: "Hide every control for a clean screenshot. Press H to toggle.",
  },
  {
    x: 1304,
    y: 72,
    name: "Advanced settings",
    what: "Further display options.",
  },
  {
    x: 1368,
    y: 72,
    name: "Copy link to this view",
    what: "A link that opens exactly what you see, image settings included.",
  },
  {
    x: 1268,
    y: 124,
    name: "Selected genes",
    what: "The molecule genes on screen. The eye hides one (⌘-click shows only it), × removes it, Clear All removes all.",
  },
  {
    x: 111,
    y: 806,
    name: "Scale bar",
    what: "The length shown, in the data's units.",
  },
];

/** The viewer as it looks, with its controls numbered. */
function ViewerMap() {
  return (
    <figure className="flex flex-col gap-4">
      <div className="relative w-full overflow-hidden rounded-xl border border-default-200">
        {/* eslint-disable-next-line @next/next/no-img-element */}
        <img
          alt="The SpatialData viewer zoomed into a skin section: the four-channel fluorescence stain underneath, cells coloured by cluster, and the molecules of DRD1, TH and A2ML1 on top"
          className="block h-auto w-full"
          height={900}
          src="/docs/viewer-overview.webp"
          width={1440}
        />
        {CONTROLS.map((c, i) => (
          <span
            key={c.name}
            aria-hidden
            className="absolute flex h-3.5 w-3.5 -translate-x-1/2 -translate-y-1/2 items-center justify-center rounded-full bg-primary text-[8px] font-bold text-white shadow ring-1 ring-black sm:h-6 sm:w-6 sm:text-xs sm:ring-2"
            style={{
              left: `${(c.x / 1440) * 100}%`,
              top: `${(c.y / 900) * 100}%`,
            }}
          >
            {i + 1}
          </span>
        ))}
      </div>
      <ol className="grid gap-x-6 gap-y-2 text-sm sm:grid-cols-2">
        {CONTROLS.map((c, i) => (
          <li key={c.name} className="flex gap-2">
            <span className="flex h-5 w-5 shrink-0 items-center justify-center rounded-full bg-primary text-[10px] font-bold text-white">
              {i + 1}
            </span>
            <span>
              <b>{c.name}</b>{" "}
              <span className="text-default-500">— {c.what}</span>
            </span>
          </li>
        ))}
      </ol>
      <figcaption className="text-xs text-default-500">
        A skin section (112,551 cells, a SpatialData store) at 50 µm scale: its
        four-channel image underneath, the cells, and three molecule genes on
        top.{" "}
        <a
          className="text-primary hover:underline"
          href="https://www.merfisheyes.com/spatialdata-viewer/from-s3?url=https%3A%2F%2Fmerfisheyes-bil.s3.us-west-2.amazonaws.com%2Fspatialdata%2Fskin_adapt_v2&img=eyJjIjpbWzEsIiMzZDdkZmYiLDAsMzA4M10sWzEsIiMwMGUwNjAiLDAsNTY4M10sWzEsIiNmZjQwZDAiLDAsMjMyNV0sWzEsIiNmZmQwMDAiLDAsOTUwNV1dfQ&ov=eyJnZW5lcyI6W1siRFJEMSIsIiMwMEZGRkYiLDEsdHJ1ZV0sWyJUSCIsIiNGRjAwRkYiLDEsdHJ1ZV0sWyJBMk1MMSIsIiNGRkZGMDAiLDEsdHJ1ZV1dfQ"
          rel="noreferrer"
          target="_blank"
        >
          Open this dataset
        </a>
        .
      </figcaption>
    </figure>
  );
}

export default function ViewerPage() {
  return (
    <>
      <DocHeader
        heading="Using the viewer"
        lead="Cells, molecules and an image in one 3D scene. Everything you change — the gene, the cell types, the image settings — is kept in the page's link."
      />

      <DocSection heading="Finding your way around" id="layout">
        <ViewerMap />
        <CompareTable
          headers={["How"]}
          labelHeader="To…"
          rows={[
            ["Rotate", "Left-drag"],
            ["Pan", "Right- or middle-drag"],
            ["Zoom", "Mouse wheel"],
            [
              "Read a cell",
              "Hover it: its cell type (and expression, when a gene is selected)",
            ],
            ["Select a cell type", "Double-click one of its cells"],
            ["Hide the controls for a screenshot", "Press H"],
          ]}
        />
      </DocSection>

      <DocSection heading="Single cell viewer" id="single-cell">
        <Prose>
          <ul className="flex list-disc flex-col gap-2 pl-5">
            <li>
              <b>Colour by gene</b> — a blue-white-red gradient, scaled to the
              95th percentile. Drag the numbers on the scale bar to change its
              range, and switch it between <b>Raw</b> and <b>Log</b> (see{" "}
              <DocLink href="/docs/expression#raw-log">
                raw counts or log values
              </DocLink>
              ).
            </li>
            <li>
              <b>Colour by a column</b> — cell types, clusters or any metadata.
              Columns with up to 100 values get one colour each and can be
              filtered; columns with more are shown as a gradient.
            </li>
            <li>
              <b>Both at once</b> — select a gene and some cell types to see the
              gene only in those cells; the rest are greyed out.
            </li>
            <li>
              <b>DEG</b> — the genes that set one cell type apart from the rest,
              or from another cell type, ranked by fold change. See{" "}
              <DocLink href="/docs/expression#deg">
                how DEGs are calculated
              </DocLink>
              .
            </li>
            <li>
              <b>Plot panel</b> — boxplots and histograms of a gene per cell
              type, cell counts per type, and histograms of numeric columns.
              Group by a second column (e.g. treatment within each cell type),
              cap the axes, and download the numbers as CSV.
            </li>
          </ul>
        </Prose>
      </DocSection>

      <DocSection heading="Single molecule viewer" id="single-molecule">
        <Prose>
          <ul className="flex list-disc flex-col gap-2 pl-5">
            <li>A few genes are shown when it opens.</li>
            <li>
              <b>Select genes</b> by searching and ticking them; each gene is
              its own colour, loaded when you pick it.
            </li>
            <li>
              <b>2D / 3D</b> — a top-down view, or full rotation.
            </li>
            <li>
              <b>Point size</b> for all genes, and for each gene on its own.
            </li>
          </ul>
        </Prose>
      </DocSection>

      <DocSection
        heading="Molecules and an image over your cells"
        id="overlays"
      >
        <Prose>
          <p>
            A single cell dataset can carry its molecules and an image. They are
            linked when you upload them together, and an owner can link them
            later with the <b>Overlays</b> button: pick one of your molecule
            datasets and one of your images.
          </p>
          <p>
            The <b>Image</b> button shows or hides the image and sets, for each
            channel, whether it is on, its colour, and its display range: values
            at or below the lower bound draw black, at or above the upper bound
            full colour.
          </p>
        </Prose>
        <Callout>
          An image can be placed only under cells that kept their original
          coordinates (usually microns).
        </Callout>
      </DocSection>

      <DocSection heading="Two datasets side by side" id="split">
        <Prose>
          <p>
            <b>Split screen</b> opens a second dataset next to the first —
            single cell or single molecule, in any combination — with the camera
            and the selection kept in step between the two.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Links keep your view" id="links">
        <Prose>
          <p>
            The page&apos;s address follows what you do: the selected gene,
            column and cell types, the molecule genes, and the image&apos;s
            channels, colours and bounds. <b>Copy link to this view</b> (top
            right) and anyone opening it sees the same thing.
          </p>
        </Prose>
      </DocSection>
    </>
  );
}
