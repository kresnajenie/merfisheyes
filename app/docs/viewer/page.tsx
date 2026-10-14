import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocSection,
  Prose,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Using the viewer — MERFISHEYES Docs",
  description:
    "Colouring by gene or cell type, plots, molecule and image overlays, split screen, and links that keep your view.",
};

/** Where things are on screen, as a labelled sketch of the viewer. */
function ViewerMap() {
  const label =
    "rounded-md border border-default-300 bg-background/90 px-2 py-1 text-[11px] leading-tight text-default-700";

  return (
    <figure className="flex flex-col gap-2">
      <div className="relative aspect-[16/10] w-full overflow-hidden rounded-xl border border-default-200 bg-black">
        {/* Cells, sketched as dots */}
        <svg
          aria-hidden
          className="absolute inset-0 h-full w-full opacity-70"
          preserveAspectRatio="none"
          viewBox="0 0 160 100"
        >
          {Array.from({ length: 220 }, (_, i) => {
            // A fixed scatter, so the sketch is the same on every render
            const rand = (n: number) => {
              const v = Math.sin(n * 12.9898) * 43758.5453;

              return v - Math.floor(v);
            };
            const hues = [
              "#5EA2EF",
              "#17c964",
              "#f5a524",
              "#FF72E1",
              "#00b7fa",
            ];

            return (
              <circle
                key={i}
                cx={45 + rand(i + 1) * 80}
                cy={20 + rand(i + 501) * 65}
                fill={hues[Math.floor(rand(i + 901) * hues.length)]}
                r={0.8}
              />
            );
          })}
        </svg>

        <div className="absolute left-2 top-2 flex flex-col gap-1 sm:left-3 sm:top-3">
          <span className={label}>Image — show/hide, channels, bounds</span>
          <span className={label}>Celltype — colour by a column</span>
          <span className={label}>Gene — colour by expression</span>
          <span className={label}>DEG — top genes of a group</span>
          <span className={label}>Plot — boxplots, histograms</span>
          <span className={label}>Split screen</span>
          <span className={label}>Point size · SC / SM layers</span>
        </div>

        <div className="absolute right-2 top-2 flex max-w-[45%] flex-col items-end gap-1 sm:right-3 sm:top-3">
          <span className={label}>
            Camera · Hide UI · Settings · Copy link to this view
          </span>
          <span className={label}>Legends — gene, scale bar, cell types</span>
        </div>

        <div className="absolute bottom-2 left-2 sm:bottom-3 sm:left-3">
          <span className={label}>Scale bar</span>
        </div>
        <div className="absolute bottom-2 right-2 sm:bottom-3 sm:right-3">
          <span className={label}>Overlays — molecules and image</span>
        </div>
      </div>
      <figcaption className="text-xs text-default-500">
        The single cell viewer. Hover a button in the app for its name.
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
              range.
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
              or from another cell type, ranked by fold change.
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
