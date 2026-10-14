import type { Metadata } from "next";

import {
  Callout,
  CompareTable,
  DocHeader,
  DocSection,
  Flow,
  Prose,
  Status,
  Steps,
} from "@/components/docs/doc-ui";

export const metadata: Metadata = {
  title: "Sharing on Explore — MERFISHEYES Docs",
  description:
    "Your account, projects, and how to submit a dataset to the public Explore page.",
};

/** Submit → review → published or not, as the app labels each state. */
function ReviewDiagram() {
  const box =
    "flex flex-col items-start gap-2 rounded-xl border border-default-200 bg-default-50/40 p-4";

  return (
    <figure className="flex flex-col gap-2">
      <div className="grid items-center gap-3 md:grid-cols-[1fr_auto_1fr_auto_1.2fr]">
        <div className={box}>
          <span className="text-sm font-semibold">
            You press Submit to Explore
          </span>
          <span className="text-xs text-default-500">
            On the dataset&apos;s page, or in the ⋯ menu on its card.
          </span>
        </div>
        <span aria-hidden className="text-center text-default-400">
          <span className="md:hidden">↓</span>
          <span className="hidden md:inline">→</span>
        </span>
        <div className={box}>
          <Status color="warning">Pending review</Status>
          <span className="text-xs text-default-500">
            An admin looks at it. Nothing is public yet.
          </span>
        </div>
        <span aria-hidden className="text-center text-default-400">
          <span className="md:hidden">↓</span>
          <span className="hidden md:inline">→</span>
        </span>
        <div className="flex flex-col gap-3">
          <div className={box}>
            <Status color="success">Published on Explore</Status>
            <span className="text-xs text-default-500">
              Listed on Explore under <b>Community</b> (and in <b>All</b>).
            </span>
          </div>
          <div className={box}>
            <Status color="danger">Not approved</Status>
            <span className="text-xs text-default-500">
              Hover the badge for the reviewer&apos;s note. Change what it asks
              for and submit again.
            </span>
          </div>
        </div>
      </div>
      <figcaption className="text-xs text-default-500">
        The badge on your dataset&apos;s card shows where it is.{" "}
        <b>Withdraw from Explore</b> takes it off again at any point.
      </figcaption>
    </figure>
  );
}

export default function SharingPage() {
  return (
    <>
      <DocHeader
        heading="Sharing on Explore"
        lead="Explore lists datasets anyone can open. Next to our curated collection there is a Community section: datasets people submitted from their own accounts and an admin approved."
      />

      <DocSection heading="The path to Explore" id="overview">
        <Flow
          steps={[
            {
              title: "Sign in",
              detail: "With Google, or with a code sent to your email.",
            },
            {
              title: "Have the dataset in your account",
              detail: "Upload it, save a preview, or add one you host.",
            },
            {
              title: "Describe it",
              detail: "Title, species, tissue, authors, a thumbnail…",
            },
            {
              title: "Submit and get reviewed",
              detail: "An admin approves it, and it appears on Explore.",
            },
          ]}
        />
      </DocSection>

      <DocSection heading="Step by step" id="steps">
        <Steps
          steps={[
            {
              title: "Sign in",
              children: (
                <p>
                  Use <b>Sign In</b> at the top right of any page, with Google
                  or with a one-time code sent to your email. Your datasets are
                  then under your name → <b>Your datasets</b>.
                </p>
              ),
            },
            {
              title: "Get the dataset into your account",
              children: (
                <>
                  <p>Any of these makes you its owner:</p>
                  <CompareTable
                    headers={["How"]}
                    labelHeader="You have"
                    rows={[
                      [
                        "Raw data",
                        <>
                          Upload it with <b>Upload &amp; process on server</b>.
                        </>,
                      ],
                      [
                        "A preview open in your browser",
                        <>
                          Press <b>Upload &amp; Save</b> in the viewer.
                        </>,
                      ],
                      [
                        "Data on your own public bucket",
                        <>
                          Open it with <b>Load from S3</b>, then{" "}
                          <b>Add to my account</b> on the banner. Nothing is
                          copied.
                        </>,
                      ],
                    ]}
                  />
                </>
              ),
            },
            {
              title: "Add its details",
              children: (
                <>
                  <p>
                    In <b>Your datasets</b>, open the dataset and press{" "}
                    <b>Edit details</b>. What you fill in is what people see on
                    Explore and what they can filter by:
                  </p>
                  <ul className="list-disc pl-5">
                    <li>
                      <b>Title</b> and <b>Description</b>
                    </li>
                    <li>
                      <b>Species</b>, <b>Tissue</b>, <b>Platform</b>,{" "}
                      <b>Disease</b>, <b>Institute</b> and <b>Tags</b>
                    </li>
                    <li>
                      <b>Source Link</b> and <b>Publication Link</b>
                    </li>
                    <li>
                      Publication details: authors, investigator,
                      co-investigators, institution, age, sex, genotype,
                      technique, funding, year, license and citation
                    </li>
                  </ul>
                  <p>
                    For a thumbnail, open the dataset in the viewer, frame a
                    view you like, and use{" "}
                    <b>Camera controls → Set thumbnail from screenshot</b>.
                  </p>
                </>
              ),
            },
            {
              title: "Submit to Explore",
              children: (
                <p>
                  On the dataset&apos;s page press <b>Submit to Explore</b>, or
                  use the ⋯ menu on its card. Any dataset you own that has
                  finished processing can be submitted.
                </p>
              ),
            },
            {
              title: "Wait for the review",
              children: (
                <p>
                  An admin reviews every submission before it is public. The
                  badge on your card follows it, as shown below.
                </p>
              ),
            },
          ]}
        />
      </DocSection>

      <DocSection heading="What happens after you submit" id="review">
        <ReviewDiagram />
        <Prose>
          <p>
            Once a dataset is published, editing its details and pressing{" "}
            <b>Submit to Explore</b> again updates its Explore listing straight
            away — it does not go back into review. A submission that is pending
            or not approved goes back into the queue.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Projects" id="projects">
        <Prose>
          <p>
            A project groups datasets that belong together — the cells,
            molecules and image of one experiment, or the sections of one study.
            In <b>Your datasets</b>, use <b>Add to project</b> on a dataset (or
            create a <b>New project</b> there). A project has its own page,
            details and thumbnail, and can be submitted to Explore as a whole,
            where it gets its own page listing its datasets.
          </p>
          <p>
            Uploading cells, molecules and an image together makes a project for
            you.
          </p>
        </Prose>
      </DocSection>

      <DocSection heading="Taking it down" id="remove">
        <Prose>
          <p>
            <b>Withdraw from Explore</b> removes a dataset or project from
            Explore; it stays in your account. To delete a dataset you uploaded,
            open it and press <b>Delete</b> — withdraw it from Explore first if
            it is there.
          </p>
        </Prose>
        <Callout kind="note">
          Questions, or a dataset you&apos;d like in the curated collection?
          Open an issue on{" "}
          <a
            href="https://github.com/kresnajenie/merfisheyes/issues"
            rel="noreferrer"
            target="_blank"
          >
            GitHub
          </a>
          .
        </Callout>
      </DocSection>
    </>
  );
}
