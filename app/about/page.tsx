import type { Metadata } from "next";

import NextLink from "next/link";

import { title } from "@/components/primitives";
import { siteConfig } from "@/config/site";
import { publicCatalogWhere } from "@/lib/explore/public-filter";
import { prisma } from "@/lib/prisma";
import packageJson from "@/package.json";

import { AboutRays } from "./about-rays";

export const metadata: Metadata = {
  title: "About — MERFISHEYES",
  description: "Who makes MERFISHEYES, how it is funded, and how to cite it.",
};

// The numbers are re-counted at most once an hour
export const revalidate = 3600;

const LAB_URL = "https://b.bintulab.com";
const PEOPLE_URL = "https://b.bintulab.com/people/";
const PREPRINT_URL =
  "https://www.biorxiv.org/content/10.64898/2026.09.20.753029v2";
const EMAIL = "ijenie@ucsd.edu";

const PIS = [
  { name: "Don W. Cleveland", role: "Principal Investigator" },
  { name: "Bogdan Bintu", role: "Principal Investigator" },
];

const TEAM = [
  { name: "Kresna Jenie", role: "Founder" },
  { name: "Evan Mishkin", role: "Lead Developer" },
  { name: "Elvis Smith", role: "Lead Data Analyst" },
];

const CITATION =
  "Jenie I, Mishkin E, Smith E, Kenney M, Cao-Berg I, Maimon R, Ropelewski AJ, Cleveland DW, Bintu B. MERFISHEYES: Web-based Visualization and Sharing Platform for Single-cell and Single-molecule Spatial Transcriptomics. bioRxiv (2026). doi:10.64898/2026.09.20.753029";

const BIBTEX = `@article{jenie2026merfisheyes,
  title   = {MERFISHEYES: Web-based Visualization and Sharing Platform
             for Single-cell and Single-molecule Spatial Transcriptomics},
  author  = {Jenie, Ignatius and Mishkin, Evan and Smith, Elvis and
             Kenney, Mariah and Cao-Berg, Ivan and Maimon, Roy and
             Ropelewski, Alexander J. and Cleveland, Don W. and Bintu, Bogdan},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.09.20.753029},
}`;

/** What Explore holds right now. Null when the database can't be reached. */
async function exploreNumbers() {
  try {
    const distinct = (field: "species" | "tissue" | "platform") =>
      prisma.catalogDataset
        .findMany({
          where: { ...publicCatalogWhere, [field]: { not: null } },
          select: { [field]: true },
          distinct: [field],
        })
        .then((rows) => rows.length);
    const [datasets, species, tissues, platforms] = await Promise.all([
      prisma.catalogDataset.count({ where: publicCatalogWhere }),
      distinct("species"),
      distinct("tissue"),
      distinct("platform"),
    ]);

    // An empty catalog (a fresh database) has nothing worth showing
    if (datasets === 0) return null;

    return [
      { value: datasets, label: "datasets on Explore" },
      { value: species, label: "species" },
      { value: tissues, label: "tissues" },
      { value: platforms, label: "platforms" },
    ];
  } catch {
    return null;
  }
}

const card = "rounded-2xl border border-default-200 bg-default-50/40";

function Section({
  heading,
  children,
}: {
  heading: string;
  children: React.ReactNode;
}) {
  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-xl font-semibold tracking-tight">{heading}</h2>
      {children}
    </section>
  );
}

function Person({ name, role }: { name: string; role: string }) {
  // First and last name, so "Don W. Cleveland" is DC
  const words = name.split(" ");
  const initials = `${words[0][0]}${words[words.length - 1][0]}`;

  return (
    <a
      className={`${card} flex items-center gap-4 p-4 transition-colors hover:border-primary/50`}
      href={PEOPLE_URL}
      rel="noreferrer"
      target="_blank"
    >
      <span className="flex h-12 w-12 shrink-0 items-center justify-center rounded-full bg-primary/15 text-sm font-semibold text-primary">
        {initials}
      </span>
      <span className="flex flex-col">
        <span className="font-semibold">{name}</span>
        <span className="text-sm text-default-500">{role}</span>
      </span>
    </a>
  );
}

const link = "text-primary hover:underline";

export default async function AboutPage() {
  const numbers = await exploreNumbers();

  return (
    <div>
      <div className="pointer-events-none fixed inset-0 z-0 h-full w-full opacity-40">
        <AboutRays />
      </div>

      <div className="relative z-10 mx-auto flex max-w-4xl flex-col gap-14 px-4 py-12 md:px-8 md:py-16">
        <header className="flex flex-col items-center gap-5 text-center">
          <h1 className={title({ size: "md" })}>
            <span className={title({ color: "green", size: "md" })}>About</span>{" "}
            MERFISHEYES
          </h1>
          <p className="max-w-2xl text-lg text-default-600">
            MERFISHEYES shows spatial transcriptomics data in 3D, right in the
            browser: cells, the molecules inside them, and the tissue image
            underneath. There is nothing to install, it is free, and anything
            you load can be shared as a link.
          </p>
          <div className="flex flex-wrap justify-center gap-3">
            <NextLink
              className="rounded-full bg-primary px-5 py-2 text-sm font-medium text-white hover:opacity-90"
              href="/explore"
            >
              Explore datasets
            </NextLink>
            <NextLink
              className="rounded-full border border-default-300 px-5 py-2 text-sm font-medium hover:border-primary/50"
              href="/docs"
            >
              Read the docs
            </NextLink>
          </div>
        </header>

        {numbers && (
          <dl className="grid grid-cols-2 gap-3 sm:grid-cols-4">
            {numbers.map((n) => (
              <div
                key={n.label}
                className={`${card} flex flex-col items-center gap-1 px-3 py-5 text-center`}
              >
                <dt className="order-2 text-xs text-default-500">{n.label}</dt>
                <dd className="order-1 text-3xl font-semibold tabular-nums">
                  {n.value.toLocaleString("en-US")}
                </dd>
              </div>
            ))}
          </dl>
        )}

        <figure className="flex flex-col gap-2">
          {/* eslint-disable-next-line @next/next/no-img-element */}
          <img
            alt="A skin section in MERFISHEYES: its four-channel stain underneath, the cells, and three genes' molecules on top"
            className="w-full rounded-2xl border border-default-200"
            height={900}
            src="/docs/viewer-overview.webp"
            width={1440}
          />
          <figcaption className="text-center text-xs text-default-500">
            A skin section at 50 µm: the tissue image, the cells and the
            molecules of three genes, together in one view.
          </figcaption>
        </figure>

        <Section heading="Why it exists">
          <div className="flex flex-col gap-3 text-default-700">
            <p>
              Spatial transcriptomics experiments produce millions of cells and
              hundreds of millions of molecules, and looking at them usually
              means a workstation, a Python environment and a lot of memory.
              Sharing them is harder still.
            </p>
            <p>
              MERFISHEYES opens them in any browser, reads only what the view
              needs, and turns what you are looking at into a link — for a
              collaborator, a reviewer, or anyone on{" "}
              <NextLink className={link} href="/explore">
                Explore
              </NextLink>
              .
            </p>
          </div>
        </Section>

        <Section heading="The team">
          <p className="text-default-600">
            MERFISHEYES is made at the{" "}
            <a className={link} href={LAB_URL} rel="noreferrer" target="_blank">
              Bintu Lab
            </a>{" "}
            at UC San Diego.
          </p>
          <div className="grid gap-3 sm:grid-cols-2">
            {PIS.map((p) => (
              <Person key={p.name} {...p} />
            ))}
          </div>
          <div className="grid gap-3 sm:grid-cols-3">
            {TEAM.map((p) => (
              <Person key={p.name} {...p} />
            ))}
          </div>
        </Section>

        <section
          className={`${card} flex flex-col items-center gap-2 border-primary/30 bg-primary/5 px-6 py-10 text-center`}
        >
          <p className="text-sm uppercase tracking-widest text-default-500">
            Funding
          </p>
          <p className="text-2xl font-semibold md:text-3xl">
            Funded by Don W. Cleveland
          </p>
          <p className="text-default-600">UC San Diego</p>
        </section>

        <Section heading="How to cite">
          <p className="text-default-700">
            If MERFISHEYES helped your work, please cite the preprint:
          </p>
          <blockquote className={`${card} px-5 py-4 text-sm leading-relaxed`}>
            {CITATION}.{" "}
            <a
              className={link}
              href={PREPRINT_URL}
              rel="noreferrer"
              target="_blank"
            >
              bioRxiv
            </a>
          </blockquote>
          <pre
            className={`${card} overflow-x-auto px-5 py-4 text-xs leading-relaxed`}
          >
            <code>{BIBTEX}</code>
          </pre>
        </Section>

        <Section heading="Acknowledgements">
          <ul className="flex list-disc flex-col gap-2 pl-5 text-default-700">
            <li>The Brain Image Library, whose datasets are on Explore.</li>
            <li>
              MERFISHEYES is open source under the MIT license, and built on
              open-source work: Next.js, Three.js, AnnData and Zarr readers, and
              the scverse SpatialData formats.
            </li>
          </ul>
        </Section>

        <Section heading="Get in touch">
          <div className="grid gap-3 sm:grid-cols-2">
            {[
              {
                label: "Email",
                detail: EMAIL,
                href: `mailto:${EMAIL}`,
              },
              {
                label: "Discord",
                detail: "Questions and feedback",
                href: siteConfig.links.discord,
              },
              {
                label: "GitHub",
                detail: "Code, issues and feature requests",
                href: siteConfig.links.github,
              },
              {
                label: "Sponsor",
                detail: "Support MERFISHEYES on Patreon",
                href: siteConfig.links.sponsor,
              },
            ].map((c) => (
              <a
                key={c.label}
                className={`${card} flex flex-col gap-1 p-4 transition-colors hover:border-primary/50`}
                href={c.href}
                rel="noreferrer"
                target="_blank"
              >
                <span className="font-semibold">{c.label}</span>
                <span className="text-sm text-default-500">{c.detail}</span>
              </a>
            ))}
          </div>
        </Section>

        <p className="text-center text-sm text-default-500">
          Version {packageJson.version} ·{" "}
          <a
            className={link}
            href={`${siteConfig.links.github}/releases`}
            rel="noreferrer"
            target="_blank"
          >
            What&apos;s new
          </a>
        </p>
      </div>
    </div>
  );
}
