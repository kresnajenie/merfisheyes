"use client";

import { Button } from "@heroui/button";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { useRouter } from "next/navigation";
import { useState } from "react";

import type { CatalogDatasetItem, CatalogDatasetEntry } from "./types";

function navigateToEntry(
  entry: CatalogDatasetEntry,
  router: ReturnType<typeof useRouter>,
) {
  const base =
    entry.datasetType === "single_molecule" ? "/sm-viewer" : "/viewer";

  if (entry.s3BaseUrl) {
    const url = entry.s3BaseUrl.replace(/\/+$/, "");

    router.push(`${base}/from-s3?url=${encodeURIComponent(url)}`);
  } else if (entry.datasetId) {
    router.push(`${base}/${entry.datasetId}`);
  }
}

// Metadata keys to display with human-readable labels
const METADATA_LABELS: Record<string, string> = {
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
};

interface Props {
  dataset: CatalogDatasetItem;
}

export function DatasetDetailPage({ dataset }: Props) {
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
    <div className="space-y-6">
      {/* Back button */}
      <Button
        className="text-default-500"
        size="sm"
        variant="light"
        onPress={() => router.push("/explore")}
      >
        ← Back to Explore
      </Button>

      {/* Header */}
      <div className="space-y-3">
        <div className="flex items-start gap-3 flex-wrap">
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

      {/* Thumbnail */}
      {dataset.thumbnailUrl && (
        <div className="rounded-xl overflow-hidden border border-default-200">
          <img
            alt={dataset.title}
            className="w-full max-h-96 object-cover"
            src={dataset.thumbnailUrl}
          />
        </div>
      )}

      {/* Entries */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        {/* Single Cell entries */}
        {scEntries.length > 0 && (
          <Card>
            <CardBody className="space-y-3">
              <h2 className="text-lg font-semibold flex items-center gap-2">
                <Chip color="primary" size="sm" variant="flat">
                  SC
                </Chip>
                Single Cell ({scEntries.length})
              </h2>
              <div className="space-y-2">
                {scEntries.map((entry) => (
                  <Button
                    key={entry.id}
                    className="w-full justify-start"
                    size="sm"
                    variant="flat"
                    onPress={() => navigateToEntry(entry, router)}
                  >
                    {entry.label}
                  </Button>
                ))}
              </div>
            </CardBody>
          </Card>
        )}

        {/* Single Molecule entries */}
        {smEntries.length > 0 && (
          <Card>
            <CardBody className="space-y-3">
              <h2 className="text-lg font-semibold flex items-center gap-2">
                <Chip color="success" size="sm" variant="flat">
                  SM
                </Chip>
                Single Molecule ({smEntries.length})
              </h2>
              <div className="space-y-2 max-h-80 overflow-y-auto">
                {smEntries.map((entry) => (
                  <Button
                    key={entry.id}
                    className="w-full justify-start"
                    size="sm"
                    variant="flat"
                    onPress={() => navigateToEntry(entry, router)}
                  >
                    {entry.label}
                  </Button>
                ))}
              </div>
            </CardBody>
          </Card>
        )}
      </div>

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
            <h2 className="text-lg font-semibold">
              Genes ({genes.length})
            </h2>
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
                {showAllGenes
                  ? "Show less"
                  : `Show all ${genes.length} genes`}
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
