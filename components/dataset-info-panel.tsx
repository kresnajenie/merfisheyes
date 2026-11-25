"use client";

import type { RefObject } from "react";
import type {
  MetadataFact,
  StoreMetadata,
  ValidationStatus,
} from "@/lib/spatialdata/metadata";

import { Button } from "@heroui/button";
import { Tooltip } from "@heroui/tooltip";
import { useRef, useState } from "react";

import { glassButton, glassPanel } from "@/components/primitives";
import { useStoreMetadata } from "@/lib/contexts/StoreMetadataContext";
import { useClickOutside } from "@/lib/hooks/useClickOutside";

const BTN = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

const VALIDATION: Record<
  ValidationStatus,
  { label: string; className: string; hint: string }
> = {
  pass: {
    label: "Validated",
    className: "bg-emerald-500/20 text-emerald-300",
    hint: "Every ontology term was checked and found",
  },
  fail: {
    label: "Validation failed",
    className: "bg-danger-500/20 text-danger-300",
    hint: "A term was checked and not found",
  },
  not_validated: {
    label: "Not validated",
    className: "bg-warning-500/20 text-warning-300",
    hint: "The terms have not been checked against their ontologies yet",
  },
};

function display(value: unknown): string {
  if (typeof value === "number") return value.toLocaleString("en-US");
  if (typeof value === "boolean") return value ? "yes" : "no";
  if (typeof value === "string") return value;

  return JSON.stringify(value, null, 1);
}

/** `mh1:d0ff37d6…dfb025`: enough to compare by eye, all of it on hover. */
const shortHash = (h: string) =>
  h.length > 24 ? `${h.slice(0, 12)}…${h.slice(-6)}` : h;

function FactValue({ fact }: { fact: MetadataFact }) {
  const text = display(fact.value);

  return (
    <span className="min-w-0 break-words">
      {fact.href && !fact.id ? (
        <a
          className="text-primary-300 hover:underline"
          href={fact.href}
          rel="noopener noreferrer"
          target="_blank"
        >
          {text}
        </a>
      ) : (
        text
      )}
      {fact.id && (
        <>
          <span className="text-default-400"> · </span>
          {fact.href ? (
            <a
              className="font-mono text-[11px] text-primary-300 hover:underline"
              href={fact.href}
              rel="noopener noreferrer"
              target="_blank"
              title="Look this term up"
            >
              {fact.id}
            </a>
          ) : (
            <span className="font-mono text-[11px] text-default-400">
              {fact.id}
            </span>
          )}
        </>
      )}
    </span>
  );
}

function Panel({
  metadata,
  onClose,
  buttonRef,
}: {
  metadata: StoreMetadata;
  onClose: () => void;
  buttonRef: RefObject<HTMLElement | null>;
}) {
  const panel = useRef<HTMLDivElement>(null);
  const { title, description, groups, sections, matrixHash } = metadata;
  const { validation, provenance, blocks } = metadata;

  useClickOutside([panel, buttonRef], onClose);

  return (
    <div
      ref={panel}
      className={`absolute top-14 right-0 z-[var(--z-panel)] max-h-[calc(100dvh-6rem)] w-[360px] max-w-[calc(100vw-2rem)] overflow-y-auto ${glassPanel()}`}
      data-testid="dataset-info-panel"
    >
      <div className="flex flex-col gap-3 p-4 text-xs">
        <div className="flex items-start justify-between gap-2">
          <span className="text-sm font-medium">About this dataset</span>
          <button
            aria-label="Close dataset info"
            className="text-default-400 hover:text-foreground"
            type="button"
            onClick={onClose}
          >
            ✕
          </button>
        </div>

        {title && <p className="text-sm font-semibold leading-snug">{title}</p>}
        {description && (
          <p className="whitespace-pre-line text-default-500">{description}</p>
        )}

        <div className="flex flex-wrap gap-1.5">
          {validation && (
            <span
              className={`rounded-full px-2 py-0.5 font-medium ${VALIDATION[validation.status].className}`}
              data-testid="dataset-info-validation"
              title={
                VALIDATION[validation.status].hint +
                (validation.at ? ` (${validation.at})` : "")
              }
            >
              {VALIDATION[validation.status].label}
            </span>
          )}
          {provenance && (
            <span
              className="rounded-full bg-default-100 px-2 py-0.5 text-default-600"
              title={[
                provenance.canonicalKey,
                provenance.input && `input: ${provenance.input}`,
                provenance.writtenAt,
              ]
                .filter(Boolean)
                .join(" · ")}
            >
              {provenance.writer === "merfisheyes"
                ? `MERFISHEYES pipeline${provenance.version ? ` ${provenance.version}` : ""}`
                : "Converted by SpatialData-DB"}
            </span>
          )}
        </div>

        {validation && validation.errors.length > 0 && (
          <ul className="list-disc pl-4 text-danger-300">
            {validation.errors.map((e) => (
              <li key={e}>{e}</li>
            ))}
          </ul>
        )}

        {groups.map((g) => (
          <section key={g.title}>
            <h3 className="mb-1 text-[11px] font-semibold uppercase tracking-wide text-default-400">
              {g.title}
            </h3>
            <dl className="grid grid-cols-[7rem_1fr] gap-x-2 gap-y-1">
              {g.facts.map((f) => (
                <div key={f.key} className="contents">
                  <dt
                    className="text-default-500"
                    title={
                      f.from === "general"
                        ? "As SpatialData-DB labels it"
                        : f.source
                          ? `From ${f.source}`
                          : undefined
                    }
                  >
                    {f.label}
                  </dt>
                  <dd>
                    <FactValue fact={f} />
                  </dd>
                </div>
              ))}
              {g.title === "Data" && matrixHash && (
                <div className="contents">
                  <dt
                    className="text-default-500"
                    title="Equal hashes are the same counts, however stored"
                  >
                    Matrix hash
                  </dt>
                  <dd className="font-mono text-[11px]" title={matrixHash}>
                    {shortHash(matrixHash)}
                  </dd>
                </div>
              )}
            </dl>
          </section>
        ))}

        {sections.length > 0 && (
          <details>
            <summary className="cursor-pointer text-[11px] font-semibold uppercase tracking-wide text-default-400">
              Sections ({sections.length})
            </summary>
            <table className="mt-1 w-full text-left tabular-nums">
              <thead className="text-default-400">
                <tr>
                  <th className="font-normal">Section</th>
                  <th className="font-normal">Sample</th>
                  <th className="text-right font-normal">Cells</th>
                  <th className="text-right font-normal">Transcripts</th>
                </tr>
              </thead>
              <tbody>
                {sections.map((s) => (
                  <tr
                    key={s.id}
                    title={
                      [
                        s.donorId && `donor ${s.donorId}`,
                        s.timepoint,
                        s.sex,
                        s.matrixHash,
                      ]
                        .filter(Boolean)
                        .join(" · ") || undefined
                    }
                  >
                    <td className="max-w-[7rem] truncate">{s.id}</td>
                    <td className="max-w-[7rem] truncate">
                      {s.sampleId ?? ""}
                    </td>
                    <td className="text-right">
                      {s.nObs !== null ? display(s.nObs) : ""}
                    </td>
                    <td className="text-right">
                      {s.nTranscripts !== null ? display(s.nTranscripts) : ""}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </details>
        )}

        {validation && Object.keys(validation.sources).length > 0 && (
          <details>
            <summary className="cursor-pointer text-[11px] font-semibold uppercase tracking-wide text-default-400">
              Checked against
            </summary>
            <dl className="mt-1 grid grid-cols-[7rem_1fr] gap-x-2 gap-y-0.5">
              {Object.entries(validation.sources).map(([k, v]) => (
                <div key={k} className="contents">
                  <dt className="text-default-500">{k}</dt>
                  <dd>{v}</dd>
                </div>
              ))}
            </dl>
          </details>
        )}

        <details data-testid="dataset-info-all">
          <summary className="cursor-pointer text-[11px] font-semibold uppercase tracking-wide text-default-400">
            All fields
          </summary>
          <div className="mt-1 flex flex-col gap-2">
            {blocks.map((b) => (
              <div key={b.name}>
                <div className="font-mono text-[11px] text-default-400">
                  {b.name}
                </div>
                <dl className="grid grid-cols-[minmax(0,9rem)_1fr] gap-x-2 gap-y-0.5">
                  {b.entries.map(([k, v]) => (
                    <div key={k} className="contents">
                      <dt className="break-words text-default-500">{k}</dt>
                      <dd className="min-w-0 whitespace-pre-wrap break-words font-mono text-[11px]">
                        {v === null || v === undefined ? "—" : display(v)}
                      </dd>
                    </div>
                  ))}
                </dl>
              </div>
            ))}
          </div>
        </details>
      </div>
    </div>
  );
}

/**
 * "About this dataset": the info button in a viewer's top-right row and the
 * panel it opens, showing what the store says about itself in its root
 * attributes (docs/specs/metadata.md §4) — the same for a store opened by
 * id, by URL or from disk. Nothing to show (no store, or a store with no
 * metadata), no button.
 */
export function DatasetInfoButton() {
  const metadata = useStoreMetadata();
  const [open, setOpen] = useState(false);
  const button = useRef<HTMLSpanElement>(null);

  if (!metadata) return null;

  return (
    <>
      <Tooltip content="About this dataset" placement="bottom">
        <span ref={button} className="inline-flex">
          <Button
            aria-label="About this dataset"
            className={`${BTN} ${glassButton()}`}
            color="default"
            data-testid="dataset-info-button"
            variant="light"
            onPress={() => setOpen((o) => !o)}
          >
            <svg
              className="w-5 h-5"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M11.25 11.25l.041-.02a.75.75 0 011.063.852l-.708 2.836a.75.75 0 001.063.853l.041-.021M21 12a9 9 0 11-18 0 9 9 0 0118 0zm-9-3.75h.008v.008H12V8.25z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </Button>
        </span>
      </Tooltip>
      {open && (
        <Panel
          buttonRef={button}
          metadata={metadata}
          onClose={() => setOpen(false)}
        />
      )}
    </>
  );
}
