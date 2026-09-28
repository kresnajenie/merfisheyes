"use client";

import type { CatalogDatasetItem } from "./types";

import { Modal, ModalContent, ModalHeader, ModalBody } from "@heroui/modal";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";
import NextLink from "next/link";

import { metadataEntries } from "./metadata-labels";

/**
 * `datasetType` is a free-form column, not an enum — the catalog already holds
 * values beyond the two the viewer routes on. Label the known pair and show
 * anything else verbatim rather than dropping or mislabelling it.
 */
function entryTypeLabel(datasetType: string): string {
  if (datasetType === "single_cell") return "SC";
  if (datasetType === "single_molecule") return "SM";

  return datasetType;
}

/**
 * Quick-look information for one catalog row — the dataset's (or project's)
 * description, provenance metadata, stats and entries, without leaving the
 * grid. Opened by the card's ⓘ button.
 *
 * Renders entirely from the card payload (CARD_SELECT), so it costs no extra
 * request. The one thing it can't show is the full gene list, which list
 * payloads deliberately omit — "View full details" goes to /explore/[id] for
 * that.
 */
export function DatasetInfoModal({
  dataset,
  isOpen,
  onClose,
}: {
  dataset: CatalogDatasetItem;
  isOpen: boolean;
  onClose: () => void;
}) {
  const metadata = (dataset.metadata ?? {}) as Record<string, unknown>;
  const details = metadataEntries(metadata);
  const entries = dataset.entries ?? [];

  const formatCount = (n: number | null) => {
    if (n == null) return null;
    if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
    if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;

    return String(n);
  };

  const facts = [
    dataset.species && { label: "Species", value: dataset.species },
    dataset.tissue && { label: "Tissue", value: dataset.tissue },
    dataset.disease && { label: "Disease", value: dataset.disease },
    dataset.platform && { label: "Platform", value: dataset.platform },
    dataset.institute && { label: "Institute", value: dataset.institute },
    dataset.numCells != null && {
      label: "Cells / molecules",
      value: formatCount(dataset.numCells)!,
    },
    dataset.numGenes != null && {
      label: "Genes",
      value: formatCount(dataset.numGenes)!,
    },
  ].filter(Boolean) as { label: string; value: string }[];

  return (
    <Modal isOpen={isOpen} scrollBehavior="inside" size="2xl" onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-2 pr-10">
          <div className="flex flex-wrap items-center gap-1.5">
            {dataset.bilCode && (
              <Chip color="secondary" size="sm" variant="flat">
                {dataset.bilCode}
              </Chip>
            )}
            {dataset.isCommunity && (
              <Chip size="sm" variant="flat">
                Community
              </Chip>
            )}
            {entries.length > 1 && (
              <Chip size="sm" variant="flat">
                {entries.length} entries
              </Chip>
            )}
          </div>
          <span className="text-lg leading-snug [overflow-wrap:anywhere]">
            {dataset.title}
          </span>
          {metadata.investigator ? (
            <span className="text-sm font-normal text-default-500">
              {String(metadata.investigator)}
              {metadata.institution ? (
                <span className="text-default-400">
                  {" "}
                  — {String(metadata.institution)}
                </span>
              ) : null}
            </span>
          ) : null}
        </ModalHeader>

        <ModalBody className="gap-5 pb-6">
          {dataset.description && (
            <p className="text-sm leading-relaxed text-default-500">
              {dataset.description}
            </p>
          )}

          {facts.length > 0 && (
            <Section title="Overview">
              <dl className="grid grid-cols-1 gap-x-8 gap-y-1.5 sm:grid-cols-2">
                {facts.map((f) => (
                  <Row key={f.label} label={f.label} value={f.value} />
                ))}
              </dl>
            </Section>
          )}

          {details.length > 0 && (
            <Section title="Details">
              <dl className="grid grid-cols-1 gap-x-8 gap-y-1.5 sm:grid-cols-2">
                {details.map((d) => (
                  <Row key={d.key} label={d.label} value={d.value} />
                ))}
              </dl>
            </Section>
          )}

          {entries.length > 0 && (
            <Section title={entries.length === 1 ? "Dataset" : "Datasets"}>
              <div className="flex flex-col gap-1">
                {entries.map((entry) => (
                  <div
                    key={entry.id}
                    className="flex items-center gap-2 text-sm"
                  >
                    <Chip
                      color={
                        entry.datasetType === "single_cell"
                          ? "primary"
                          : entry.datasetType === "single_molecule"
                            ? "secondary"
                            : "default"
                      }
                      size="sm"
                      variant="flat"
                    >
                      {entryTypeLabel(entry.datasetType)}
                    </Chip>
                    <span className="truncate" title={entry.label}>
                      {entry.label}
                    </span>
                  </div>
                ))}
              </div>
            </Section>
          )}

          {dataset.tags.length > 0 && (
            <Section title="Tags">
              <div className="flex flex-wrap gap-1">
                {dataset.tags.map((tag) => (
                  <Chip key={tag} size="sm" variant="dot">
                    {tag}
                  </Chip>
                ))}
              </div>
            </Section>
          )}

          <div className="flex flex-wrap gap-2 pt-1">
            <Button
              as={NextLink}
              color="primary"
              href={`/explore/${dataset.id}`}
              size="sm"
              variant="flat"
            >
              View full details
            </Button>
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
            {dataset.externalLink && (
              <Button
                as="a"
                href={dataset.externalLink}
                rel="noopener noreferrer"
                size="sm"
                target="_blank"
                variant="flat"
              >
                Source
              </Button>
            )}
          </div>
        </ModalBody>
      </ModalContent>
    </Modal>
  );
}

function Section({
  title,
  children,
}: {
  title: string;
  children: React.ReactNode;
}) {
  return (
    <section className="space-y-2">
      <h3 className="text-xs font-medium uppercase tracking-wider text-default-400">
        {title}
      </h3>
      {children}
    </section>
  );
}

function Row({ label, value }: { label: string; value: string }) {
  return (
    <div className="flex gap-2 py-0.5">
      <dt className="w-[130px] shrink-0 text-sm text-default-400">{label}</dt>
      <dd className="text-sm [overflow-wrap:anywhere]">{value}</dd>
    </div>
  );
}
