"use client";

import type { DatasetRow } from "./types";
import type { CatalogDatasetItem } from "@/components/explore/types";

import { Chip } from "@heroui/chip";
import { Button } from "@heroui/button";
import {
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
} from "@heroui/dropdown";
import { useRouter } from "next/navigation";

import { SubmissionBadge } from "./submission-badge";

import { ExploreDatasetCard } from "@/components/explore/explore-dataset-card";

interface AccountDatasetCardProps {
  dataset: DatasetRow;
  projectNames: string[];
  /** Active gene search text — makes the card show its matched-gene chips. */
  geneHighlight?: string;
  onEdit: (dataset: DatasetRow) => void;
  onSubmit: (dataset: DatasetRow) => void;
  onWithdraw: (dataset: DatasetRow) => void;
  onAddToProject: (dataset: DatasetRow) => void;
}

/** Map an owned dataset onto the shared Explore card shape (single entry). */
function toItem(d: DatasetRow): CatalogDatasetItem {
  return {
    id: d.id,
    title: d.title ?? "Untitled dataset",
    description: d.description,
    species: d.species,
    disease: d.disease,
    institute: d.institute,
    tissue: d.tissue,
    platform: d.platform,
    tags: d.tags ?? [],
    matchedGenes: d.matchedGenes,
    thumbnailUrl: d.thumbnailUrl,
    bilCode: null,
    metadata: (d.metadata ?? {}) as Record<string, unknown>,
    externalLink: d.externalLink,
    publicationLink: d.publicationLink,
    entries: [
      {
        id: d.id,
        label: d.title ?? "Dataset",
        // Legacy rows store raw formats ("h5ad", "merscope") here; anything
        // not single_molecule renders as single cell (badge + viewer route).
        datasetType:
          d.datasetType === "single_molecule"
            ? "single_molecule"
            : "single_cell",
        s3BaseUrl: null,
        datasetId: d.id,
        thumbnailUrl: d.thumbnailUrl,
        sortOrder: 0,
      },
    ],
    isPublished: false,
    isFeatured: false,
    isBil: false,
    isInternal: false,
    isCommunity: false,
    sortOrder: 0,
    numCells: d.numCells,
    numGenes: d.numGenes,
  };
}

export function AccountDatasetCard({
  dataset,
  projectNames,
  geneHighlight,
  onEdit,
  onSubmit,
  onWithdraw,
  onAddToProject,
}: AccountDatasetCardProps) {
  const router = useRouter();
  const sub = dataset.submission;
  const canSubmit = dataset.status === "COMPLETE";

  const overlay = (
    <>
      {sub && <SubmissionBadge submission={sub} />}
      <Dropdown placement="bottom-end">
        <DropdownTrigger>
          <Button
            isIconOnly
            aria-label="Dataset actions"
            className="min-w-7 w-7 h-7 bg-default-800/70 text-white backdrop-blur"
            size="sm"
            variant="flat"
          >
            ⋯
          </Button>
        </DropdownTrigger>
        <DropdownMenu aria-label="Dataset actions">
          <DropdownItem key="edit" onPress={() => onEdit(dataset)}>
            Edit details
          </DropdownItem>
          <DropdownItem key="add" onPress={() => onAddToProject(dataset)}>
            Add to project
          </DropdownItem>
          {sub ? (
            <DropdownItem
              key="withdraw"
              className="text-danger"
              color="danger"
              onPress={() => onWithdraw(dataset)}
            >
              Withdraw from Explore
            </DropdownItem>
          ) : canSubmit ? (
            <DropdownItem key="submit" onPress={() => onSubmit(dataset)}>
              Submit to Explore
            </DropdownItem>
          ) : null}
        </DropdownMenu>
      </Dropdown>
    </>
  );

  // Project-membership badges, overlaid on the header (bottom-left) so they
  // never collide with the card's stats text.
  const projectBadges =
    projectNames.length > 0 ? (
      <>
        {projectNames.slice(0, 2).map((name) => (
          <Chip
            key={name}
            className="max-w-[130px] bg-default-800/70 text-white backdrop-blur"
            classNames={{ content: "truncate" }}
            size="sm"
            variant="flat"
          >
            {name}
          </Chip>
        ))}
        {projectNames.length > 2 && (
          <Chip
            className="bg-default-800/70 text-white backdrop-blur"
            size="sm"
            variant="flat"
          >
            +{projectNames.length - 2}
          </Chip>
        )}
      </>
    ) : null;

  return (
    <ExploreDatasetCard
      dataset={toItem(dataset)}
      geneHighlight={geneHighlight}
      headerBadges={projectBadges}
      ownerOverlay={overlay}
      onCardClick={() => router.push(`/account/datasets/${dataset.id}`)}
    />
  );
}
