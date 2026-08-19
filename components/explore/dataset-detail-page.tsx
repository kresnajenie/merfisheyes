"use client";

import type { CatalogDatasetItem, CatalogDatasetEntry } from "./types";

import { Button } from "@heroui/button";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { useRouter } from "next/navigation";
import { useState } from "react";

// Metadata keys to display with human-readable labels
const METADATA_LABELS: Record<string, string> = {
  authors: "Authors",
  investigator: "Investigator",
  institution: "Institution",
  coInvestigators: "Co-Investigators",
  funding: "Funding",
  publicationYear: "Year",
  license: "License",
  age: "Age",
  sex: "Sex",
  genotype: "Genotype",
  technique: "Technique",
  citation: "Citation",
};

const COLOR_STYLES = {
  blue: {
    border: "border-blue-500/30 hover:border-blue-400",
    bg: "hover:bg-blue-500/10",
    accent: "bg-blue-500",
    text: "text-blue-400",
    arrow: "text-blue-400/50 group-hover:text-blue-400",
  },
  purple: {
    border: "border-purple-500/30 hover:border-purple-400",
    bg: "hover:bg-purple-500/10",
    accent: "bg-purple-500",
    text: "text-purple-400",
    arrow: "text-purple-400/50 group-hover:text-purple-400",
  },
};

function getEntryHref(entry: CatalogDatasetEntry): string {
  const base =
    entry.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";

  if (entry.s3BaseUrl) {
    const url = entry.s3BaseUrl.replace(/\/+$/, "");

    return `${base}/from-s3?url=${encodeURIComponent(url)}`;
  }
  if (entry.datasetId) {
    return `${base}/${entry.datasetId}`;
  }

  return "#";
}

/**
 * Thumbnail for an entry. App datasets resolve LIVE through the thumbnail
 * proxy — `/api/datasets/{id}/thumbnail` always serves the dataset's current
 * image — so an owner changing their thumbnail is reflected here without
 * re-submitting. The entry's stored thumbnailUrl remains the source for
 * BIL/S3 entries (static images, no datasetId) and as a fallback.
 */
function entryThumbnailSrc(entry: CatalogDatasetEntry): string | null {
  if (entry.datasetId) return `/api/datasets/${entry.datasetId}/thumbnail`;

  return entry.thumbnailUrl;
}

function EntryCard({
  entry,
  color,
  large,
  onSelect,
  actions,
}: {
  entry: CatalogDatasetEntry;
  color: "blue" | "purple";
  large?: boolean;
  onSelect?: (entry: CatalogDatasetEntry) => void;
  /** Overlay controls (e.g. an owner's Remove button), top-right of the card. */
  actions?: React.ReactNode;
}) {
  const s = COLOR_STYLES[color];
  const className = `group rounded-xl border ${s.border} ${s.bg} overflow-hidden transition-all block w-full text-left`;
  const thumbSrc = entryThumbnailSrc(entry);
  const body = (
    <>
      {thumbSrc ? (
        <div
          className={`relative w-full ${large ? "aspect-[4/3]" : "aspect-[16/10]"} bg-default-100`}
        >
          <img
            alt={entry.label}
            className="w-full h-full object-cover"
            src={thumbSrc}
            onError={(e) => {
              // Dataset has no thumbnail (proxy 404s) — collapse to the
              // accent strip instead of a broken-image icon.
              e.currentTarget.parentElement!.style.display = "none";
            }}
          />
          <div
            className={`absolute top-2 left-2 w-2 h-2 rounded-full ${s.accent}`}
          />
        </div>
      ) : (
        <div className={`w-full ${large ? "h-2" : "h-1"} ${s.accent}`} />
      )}
      <div
        className={`flex items-center justify-between gap-2 ${large ? "px-4 py-3" : "px-3 py-2.5"}`}
      >
        <span
          className={`${large ? "text-base font-medium" : "text-sm"} truncate`}
        >
          {entry.label}
        </span>
        <svg
          className={`w-3.5 h-3.5 shrink-0 transition-colors ${s.arrow}`}
          fill="none"
          stroke="currentColor"
          strokeWidth={2}
          viewBox="0 0 24 24"
        >
          <path
            d="M13 7l5 5-5 5M6 12h12"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      </div>
    </>
  );

  const card = onSelect ? (
    <button className={className} type="button" onClick={() => onSelect(entry)}>
      {body}
    </button>
  ) : (
    <a className={className} href={getEntryHref(entry)}>
      {body}
    </a>
  );

  if (!actions) return card;

  // Overlay is a sibling of the link/button so clicking it doesn't navigate.
  return (
    <div className="relative">
      {card}
      <div className="absolute top-2 right-2 z-10">{actions}</div>
    </div>
  );
}

interface Props {
  dataset: CatalogDatasetItem;
  /** When provided, entry clicks call this instead of navigating to the viewer. */
  onSelectEntry?: (entry: CatalogDatasetEntry) => void;
  /** When provided, the back button calls this instead of router.push("/explore"). */
  onBack?: () => void;
  /** Overrides the back-button label (default "← Back to Explore"). */
  backLabel?: string;
  /** Owner controls (Edit / Submit / …) rendered top-right of the header. */
  headerActions?: React.ReactNode;
  /** Per-entry overlay controls (e.g. an owner's Remove button). */
  entryActions?: (entry: CatalogDatasetEntry) => React.ReactNode;
}

export function DatasetDetailPage({
  dataset,
  onSelectEntry,
  onBack,
  backLabel = "← Back to Explore",
  headerActions,
  entryActions,
}: Props) {
  const router = useRouter();
  const [showAllGenes, setShowAllGenes] = useState(false);

  const metadata = (dataset.metadata ?? {}) as Record<string, unknown>;
  const genes = dataset.genes ?? [];
  const scEntries = dataset.entries.filter(
    (e) => e.datasetType === "single_cell",
  );
  const smEntries = dataset.entries.filter(
    (e) => e.datasetType === "single_molecule",
  );

  const displayGenes = showAllGenes ? genes : genes.slice(0, 50);

  return (
    <div className="space-y-4">
      {/* Back button */}
      <Button
        className="text-default-500"
        size="sm"
        variant="light"
        onPress={() => (onBack ? onBack() : router.push("/explore"))}
      >
        {backLabel}
      </Button>

      {/* Header */}
      <div className="space-y-2">
        {headerActions && (
          <div className="flex flex-wrap justify-end gap-2">
            {headerActions}
          </div>
        )}
        <div className="flex items-start gap-2 flex-wrap">
          {dataset.bilCode && (
            <Chip color="secondary" size="sm" variant="flat">
              {dataset.bilCode}
            </Chip>
          )}
          {dataset.species && (
            <Chip color="primary" size="sm" variant="flat">
              {dataset.species}
            </Chip>
          )}
          {dataset.tissue && (
            <Chip size="sm" variant="flat">
              {dataset.tissue}
            </Chip>
          )}
          {dataset.platform && (
            <Chip size="sm" variant="flat">
              {dataset.platform}
            </Chip>
          )}
        </div>
        <h1 className="text-2xl font-bold">{dataset.title}</h1>
        {metadata.investigator ? (
          <p className="text-default-500 text-lg">
            {String(metadata.investigator)}
            {metadata.institution ? (
              <span className="text-default-400">
                {" "}
                — {String(metadata.institution)}
              </span>
            ) : null}
          </p>
        ) : null}
        {dataset.description && (
          <p className="text-default-400 leading-relaxed">
            {dataset.description}
          </p>
        )}
      </div>

      {/* Entries — full-width horizontal sections: SC rows first, SM below.
          Each section flows its cards in responsive rows so projects with
          many datasets of one type use the whole page width instead of a
          single narrow column. */}
      {scEntries.length > 0 && (
        <div className="space-y-2">
          <h2 className="text-sm font-medium text-default-500 uppercase tracking-wider">
            Single Cell{scEntries.length > 1 ? ` (${scEntries.length})` : ""}
          </h2>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-3">
            {scEntries.map((entry) => (
              <EntryCard
                key={entry.id}
                large
                actions={entryActions?.(entry)}
                color="blue"
                entry={entry}
                onSelect={onSelectEntry}
              />
            ))}
          </div>
        </div>
      )}

      {smEntries.length > 0 && (
        <div className="space-y-2">
          <h2 className="text-sm font-medium text-default-500 uppercase tracking-wider">
            Single Molecule ({smEntries.length})
          </h2>
          <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-4 xl:grid-cols-5 gap-2">
            {smEntries.map((entry) => (
              <EntryCard
                key={entry.id}
                actions={entryActions?.(entry)}
                color="purple"
                entry={entry}
                onSelect={onSelectEntry}
              />
            ))}
          </div>
        </div>
      )}

      {/* Metadata */}
      {Object.keys(metadata).length > 0 && (
        <Card>
          <CardBody className="space-y-3">
            <h2 className="text-lg font-semibold">Details</h2>
            <div className="grid grid-cols-1 md:grid-cols-2 gap-x-8 gap-y-2">
              {Object.entries(metadata).map(([key, value]) => {
                if (!value) return null;
                const label = METADATA_LABELS[key] || key;
                const display = Array.isArray(value)
                  ? value.join(", ")
                  : String(value);

                return (
                  <div key={key} className="flex gap-2 py-1">
                    <span className="text-default-400 text-sm min-w-[120px] shrink-0">
                      {label}
                    </span>
                    <span className="text-sm">{display}</span>
                  </div>
                );
              })}
            </div>
          </CardBody>
        </Card>
      )}

      {/* Genes */}
      {genes.length > 0 && (
        <Card>
          <CardBody className="space-y-3">
            <h2 className="text-lg font-semibold">Genes ({genes.length})</h2>
            <div className="flex flex-wrap gap-1.5">
              {displayGenes.map((gene) => (
                <Chip key={gene} size="sm" variant="flat">
                  {gene}
                </Chip>
              ))}
            </div>
            {genes.length > 50 && (
              <Button
                size="sm"
                variant="light"
                onPress={() => setShowAllGenes(!showAllGenes)}
              >
                {showAllGenes ? "Show less" : `Show all ${genes.length} genes`}
              </Button>
            )}
          </CardBody>
        </Card>
      )}

      {/* External links */}
      {(dataset.externalLink || dataset.publicationLink) && (
        <div className="flex gap-3">
          {dataset.externalLink && (
            <Button
              as="a"
              color="primary"
              href={dataset.externalLink}
              rel="noopener noreferrer"
              size="sm"
              target="_blank"
              variant="flat"
            >
              View on BIL
            </Button>
          )}
          {dataset.publicationLink && (
            <Button
              as="a"
              href={dataset.publicationLink}
              rel="noopener noreferrer"
              size="sm"
              target="_blank"
              variant="flat"
            >
              Publication
            </Button>
          )}
        </div>
      )}
    </div>
  );
}
