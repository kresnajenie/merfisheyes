import type { Metadata } from "next";

import { siteConfig } from "@/config/site";
import packageJson from "@/package.json";

import { AboutRays } from "./about-rays";

export const metadata: Metadata = {
  title: "About — MERFISHEYES",
  description: "Who makes MERFISHEYES, how it is funded, and how to cite it.",
};

const LAB_URL = "https://b.bintulab.com";
const CLEVELAND_LAB_URL =
  "https://cmm.ucsd.edu/research/labs/don-cleveland/index.html";
const PEOPLE_URL = "https://b.bintulab.com/people/";
const PREPRINT_URL =
  "https://www.biorxiv.org/content/10.64898/2026.09.20.753029v2";
const EMAIL = "ijenie@ucsd.edu";

interface Member {
  name: string;
  /** Where the card leads: a team member's LinkedIn, a PI's lab page. */
  href?: string;
  /** Other names the person goes by. */
  aka?: string[];
  role: string;
  photo: string;
}

const TEAM: Member[] = [
  {
    name: "Ignatius Kresnathan Sjahnir Jenie",
    href: "https://www.linkedin.com/in/ignatius-jenie-1023521b3/",
    aka: ["KJ", "Jenie", "Kresna"],
    role: "Founder",
    photo: "/about/jenie.webp",
  },
  {
    name: "Evan Mishkin",
    href: "https://www.linkedin.com/in/evan-mishkin-1a76a4190/",
    role: "Lead Developer",
    photo: "/about/mishkin.webp",
  },
  {
    name: "Elvis Smith",
    href: "https://www.linkedin.com/in/hire-elvis-smith/",
    role: "Lead Data Analyst",
    photo: "/about/smith.webp",
  },
];

const PIS: Member[] = [
  {
    name: "Don W. Cleveland",
    href: CLEVELAND_LAB_URL,
    role: "Principal Investigator",
    photo: "/about/cleveland.webp",
  },
  {
    name: "Bogdan Bintu",
    href: LAB_URL,
    role: "Principal Investigator",
    photo: "/about/bintu.webp",
  },
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

const card = "rounded-2xl border border-default-200 bg-default-50/40";
const link = "text-primary hover:underline";

function Section({
  heading,
  children,
}: {
  heading: string;
  children: React.ReactNode;
}) {
  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-2xl font-semibold tracking-tight">{heading}</h2>
      {children}
    </section>
  );
}

/** A square photo, the name large, the role small beneath it. */
function Person({ name, aka, href, role, photo }: Member) {
  const body = (
    <>
      {/* eslint-disable-next-line @next/next/no-img-element */}
      <img
        alt={name}
        className="aspect-square w-full rounded-xl object-cover transition-opacity group-hover:opacity-90"
        height={600}
        loading="lazy"
        src={photo}
        width={600}
      />
      <span className="flex flex-col gap-1">
        <span className="text-xl font-medium tracking-tight md:text-2xl">
          {name}
        </span>
        {aka && (
          <span className="text-sm text-default-500">
            a.k.a. {aka.map((a) => `"${a}"`).join(", ")}
          </span>
        )}
        <span className="font-mono text-xs uppercase tracking-wider text-primary">
          {role}
        </span>
      </span>
    </>
  );

  return href ? (
    <a
      className="group flex flex-col gap-3"
      href={href}
      rel="noreferrer"
      target="_blank"
    >
      {body}
    </a>
  ) : (
    <div className="flex flex-col gap-3">{body}</div>
  );
}

export default function AboutPage() {
  return (
    <div>
      <div className="pointer-events-none fixed inset-0 z-0 h-full w-full opacity-40">
        <AboutRays />
      </div>

      <div className="relative z-10 mx-auto flex max-w-5xl flex-col gap-16 px-4 py-12 md:px-8 md:py-16">
        <section className="flex flex-col gap-10">
          <div className="grid gap-4 md:grid-cols-2 md:gap-10">
            <h1 className="text-4xl font-semibold tracking-tight md:text-5xl">
              Team
            </h1>
            <p className="text-lg text-default-600">
              MERFISHEYES is made at the{" "}
              <a
                className={link}
                href={CLEVELAND_LAB_URL}
                rel="noreferrer"
                target="_blank"
              >
                Cleveland Lab
              </a>{" "}
              and the{" "}
              <a
                className={link}
                href={LAB_URL}
                rel="noreferrer"
                target="_blank"
              >
                Bintu Lab
              </a>{" "}
              at UC San Diego: a browser-based viewer for spatial
              transcriptomics, with cells, molecules and the tissue image in one
              shareable view.
            </p>
          </div>

          <div className="grid grid-cols-2 gap-x-6 gap-y-10 md:grid-cols-3">
            {TEAM.map((p) => (
              <Person key={p.name} {...p} />
            ))}
          </div>

          <div className="flex flex-col gap-6">
            <h2 className="text-sm font-medium uppercase tracking-widest text-default-500">
              Principal Investigators
            </h2>
            <div className="grid grid-cols-2 gap-x-6 gap-y-10 md:grid-cols-3">
              {PIS.map((p) => (
                <Person key={p.name} {...p} />
              ))}
            </div>
          </div>
        </section>

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
              { label: "Email", detail: EMAIL, href: `mailto:${EMAIL}` },
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

        <div className="flex flex-col items-center gap-2 text-center text-xs text-default-500">
          <p className="text-sm">
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
          <p>
            Photo of Don W. Cleveland by Flynn Larsen,{" "}
            <a
              className={link}
              href="https://commons.wikimedia.org/wiki/File:Cleveland_Don.jpg"
              rel="noreferrer"
              target="_blank"
            >
              CC BY-SA 4.0
            </a>
            , cropped. Other photos from the{" "}
            <a
              className={link}
              href={PEOPLE_URL}
              rel="noreferrer"
              target="_blank"
            >
              Bintu Lab
            </a>
            .
          </p>
        </div>
      </div>
    </div>
  );
}
