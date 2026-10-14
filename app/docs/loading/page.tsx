import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  Decision,
  DocHeader,
  DocSection,
  Flow,
  Prose,
  Steps,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Loading data — MERFISHEYES Docs",
  description:
    "Preview in your browser, upload to be processed on our servers, or open data from a public bucket.",
};

export default function LoadingPage() {
  return (
    <>
      <DocHeader
        heading="Loading data"
        lead="On the homepage, pick single cell or single molecule, then choose how to load: in your browser, on our servers, or straight from a public bucket."
      />

      <DocSection heading="Which way should I load my data?" id="choose">
        <Decision
          otherwise={{
            answer: "Upload & process on server",
            href: "#server",
          }}
          questions={[
            {
              question: "Is it already on a public S3 bucket?",
              answer: "Load from S3",
              href: "#s3",
            },
            {
              question:
                "Is it already processed — merfisheyes output folders, a SpatialData store, an AnnData table in Zarr, an OME-Zarr image?",
              answer: "Preview in browser — opens as it is",
              href: "#processed",
            },
            {
              question:
                "Does it fit in your browser — roughly under 500,000 cells, or tens of millions of molecules?",
              answer: "Preview in browser",
              href: "#preview",
            },
          ]}
        />
      </DocSection>

      <DocSection heading="The three ways, side by side" id="compare">
        <CompareTable
          headers={[
            "Preview in browser",
            "Upload & process on server",
            "Load from S3",
          ]}
          rows={[
            [
              "Where it is processed",
              "In your browser; nothing leaves your computer",
              "On our servers, by a background worker",
              "Not at all: read from the bucket as you look",
            ],
            [
              "Sign in",
              "Not to look; yes to save it",
              "Yes",
              "Not to look; yes to add it to your account",
            ],
            [
              "Size",
              "What your browser's memory allows (already-processed data: any size)",
              "No limit from your browser",
              "Any size",
            ],
            [
              "You end up with",
              "A view now; Upload & Save gives a link that lasts",
              "An email with a link when processing is done",
              "A link that reopens it",
            ],
            [
              "Best for",
              "A quick look; data that is already processed",
              "Large raw datasets",
              "Data you already host publicly",
            ],
          ]}
        />
      </DocSection>

      <DocSection heading="Preview in browser" id="preview">
        <Prose>
          <p>
            The file is parsed entirely in your browser and opens straight away.
            Nothing leaves your computer. It suits quick looks and datasets that
            fit in your browser&apos;s memory — roughly under 500,000 cells or
            tens of millions of molecules.
          </p>
          <p>
            To keep it, press <b>Upload &amp; Save</b> in the viewer (you need
            to be signed in). It is then in your account with a link that lasts.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Upload & process on server" id="server">
        <Prose>
          <p>
            For large datasets, upload the raw files and let a background worker
            process them. You need to be signed in, so we know where to email
            you.
          </p>
        </Prose>
        <Steps
          steps={[
            {
              title: "Drop your file or folder",
              children: (
                <p>
                  Switch the homepage to <b>Upload &amp; process on server</b>{" "}
                  and drop it as you would for a preview, then confirm the
                  dataset name. Single cell data can also take a cell type CSV
                  you made yourself, e.g. with MapMyCells.
                </p>
              ),
            },
            {
              title: "Confirm the columns (molecules only)",
              children: (
                <p>
                  A preview of the first rows, with a dropdown for each of gene,
                  x, y, z and cell id, filled in from the detected platform.
                  This is what makes sure your columns are read correctly.
                </p>
              ),
            },
            {
              title: "Wait for the upload",
              children: (
                <p>
                  A progress bar follows the upload. You can leave once it says{" "}
                  <b>Uploaded — we&apos;ll email you</b>.
                </p>
              ),
            },
            {
              title: "Get the email",
              children: (
                <p>
                  The worker processes the data and emails you a link to the
                  viewer. The dataset is already in your account.
                </p>
              ),
            },
          ]}
        />
        <Callout>
          <b>Cells, molecules and an image in one go.</b> Drop a single cell
          folder that also holds a transcripts file (a typical Xenium export)
          and/or an OME-Zarr image folder. Each becomes its own dataset, grouped
          into one project, and the cells open with the molecules over them and
          the image under them.
        </Callout>
      </DocSection>

      <DocSection
        heading="Already-processed data opens as it is"
        id="processed"
      >
        <Prose>
          <p>
            If your data is already in a format the viewer reads, nothing is
            converted. Drop the folder on the single cell <b>Folder</b> card in{" "}
            <b>Preview in browser</b> and it opens at once, whatever its size:
            each part is read from disk a piece at a time, as you look at it.
            Cells, molecules and an image in the same folder are overlaid.
            Folder names don&apos;t matter — each part is recognised by what is
            inside it.
          </p>
        </Prose>
        <Flow
          steps={[
            {
              title: "Cells",
              detail:
                "A SpatialData store, an AnnData table in Zarr (v2 or v3), or a merfisheyes cell folder",
            },
            {
              title: "+ Molecules",
              detail:
                "A merfisheyes molecule folder (manifest.json.gz + genes/)",
            },
            {
              title: "+ Image",
              detail:
                "An OME-Zarr image, or the one inside the SpatialData store",
            },
            {
              title: "One overlaid view",
              detail: "Molecules over the cells, the image under both",
            },
          ]}
        />
        <Prose>
          <p>
            Any part may be missing. <b>Upload &amp; Save</b> from this preview
            uploads every part exactly as it is — nothing is processed — each as
            its own dataset in one project, and the saved link opens with the
            same overlays. Raw cells (H5AD, Xenium, MERSCOPE) beside processed
            parts are shown on their own here; use{" "}
            <b>Upload &amp; process on server</b> to combine those.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Load from S3" id="s3">
        <Prose>
          <p>
            Data already on a public bucket opens straight from it. Click{" "}
            <b>Load from S3</b> on the homepage and paste one of:
          </p>
          <ul className="list-disc pl-5">
            <li>a merfisheyes cell or molecule folder;</li>
            <li>
              a SpatialData store (the <code>.zarr</code> folder);
            </li>
            <li>an OME-Zarr image;</li>
            <li>
              the folder holding a store with its molecules beside it — it opens
              with everything overlaid.
            </li>
          </ul>
          <p>
            Nothing is copied. A bucket&apos;s folders can&apos;t be listed from
            a browser, so a folder names its store in its{" "}
            <code>mapping.json</code> (
            <code>&quot;cells&quot;: &quot;name.zarr&quot;</code>
            ); without that, a store named after the folder (
            <code>sample/sample.zarr</code>) is tried.
          </p>
          <p>
            The bucket must allow public reads and allow this site in its CORS
            settings (including the <code>Range</code> header, for stores
            written by recent anndata). A dataset opened this way that nobody
            owns yet shows <b>Add to my account</b>.
          </p>
        </Prose>
      </DocSection>

      <DocSection
        heading="Very large datasets on your own machine"
        id="scripts"
      >
        <Prose>
          <p>
            You can also process data yourself with the Python scripts in the{" "}
            <a
              href="https://github.com/kresnajenie/merfisheyes/tree/main/scripts"
              rel="noreferrer"
              target="_blank"
            >
              GitHub repository
            </a>{" "}
            and drop the output folder on the homepage, where it opens as it is:
          </p>
          <pre className="overflow-x-auto rounded-xl border border-default-200 bg-default-50/40 p-3 text-xs">
            {`python scripts/process_spatial_data.py    data.h5ad    output/   # cells
python scripts/process_single_molecule.py data.parquet output/   # molecules`}
          </pre>
        </Prose>
      </DocSection>
    </>
  );
}
