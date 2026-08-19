"use client";

import type { ProjectRow } from "./types";
import type { CatalogDatasetItem } from "@/components/explore/types";

import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";
import {
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
} from "@heroui/dropdown";
import { useRouter } from "next/navigation";

import { SubmissionBadge } from "./submission-badge";

import { ExploreDatasetCard } from "@/components/explore/explore-dataset-card";

interface AccountProjectCardProps {
  project: ProjectRow;
  onEdit: (project: ProjectRow) => void;
  onSubmit: (project: ProjectRow) => void;
  onWithdraw: (project: ProjectRow) => void;
  onDelete: (project: ProjectRow) => void;
}

/** Map a project onto the shared Explore card shape (member datasets = entries). */
function toItem(p: ProjectRow): CatalogDatasetItem {
  const thumb =
    p.thumbnailUrl ??
    p.datasets.find((d) => d.thumbnailUrl)?.thumbnailUrl ??
    null;

  return {
    id: p.id,
    title: p.title,
    description: p.description,
    species: p.species,
    disease: p.disease,
    institute: p.institute,
    tissue: p.tissue,
    platform: p.platform,
    tags: p.tags ?? [],
    genes: [],
    thumbnailUrl: thumb,
    bilCode: null,
    metadata: (p.metadata ?? {}) as Record<string, unknown>,
    externalLink: p.externalLink,
    publicationLink: p.publicationLink,
    entries: p.datasets.map((d, i) => ({
      id: d.id,
      label: d.title ?? "Dataset",
      // Legacy rows store raw formats ("h5ad", "merscope") here; anything
      // not single_molecule renders as single cell.
      datasetType:
        d.datasetType === "single_molecule" ? "single_molecule" : "single_cell",
      s3BaseUrl: null,
      datasetId: d.id,
      thumbnailUrl: d.thumbnailUrl,
      sortOrder: i,
    })),
    isPublished: false,
    isFeatured: false,
    isBil: false,
    isInternal: false,
    isCommunity: false,
    sortOrder: 0,
    numCells: null,
    numGenes: null,
  };
}

export function AccountProjectCard({
  project,
  onEdit,
  onSubmit,
  onWithdraw,
  onDelete,
}: AccountProjectCardProps) {
  const router = useRouter();
  const sub = project.submission;
  const canSubmit = project.datasetCount > 0;

  const overlay = (
    <>
      {sub && <SubmissionBadge submission={sub} />}
      <Dropdown placement="bottom-end">
        <DropdownTrigger>
          <Button
            isIconOnly
            aria-label="Project actions"
            className="min-w-7 w-7 h-7 bg-default-800/70 text-white backdrop-blur"
            size="sm"
            variant="flat"
          >
            ⋯
          </Button>
        </DropdownTrigger>
        <DropdownMenu aria-label="Project actions">
          <DropdownItem key="edit" onPress={() => onEdit(project)}>
            Edit details
          </DropdownItem>
          {canSubmit && !sub ? (
            <DropdownItem key="submit" onPress={() => onSubmit(project)}>
              Submit to Explore
            </DropdownItem>
          ) : null}
          {canSubmit && sub ? (
            <DropdownItem key="resubmit" onPress={() => onSubmit(project)}>
              {sub.reviewStatus === "REJECTED"
                ? "Resubmit to Explore"
                : "Update submission"}
            </DropdownItem>
          ) : null}
          {sub ? (
            <DropdownItem key="withdraw" onPress={() => onWithdraw(project)}>
              Withdraw from Explore
            </DropdownItem>
          ) : null}
          <DropdownItem
            key="delete"
            className="text-danger"
            color="danger"
            onPress={() => onDelete(project)}
          >
            Delete project
          </DropdownItem>
        </DropdownMenu>
      </Dropdown>
    </>
  );

  // The owner overlay suppresses the card's built-in entry-count badge, so
  // surface the dataset count on the header's bottom-left instead.
  const countBadge = (
    <Chip
      className="bg-default-800/70 text-white backdrop-blur"
      size="sm"
      variant="flat"
    >
      {project.datasetCount} dataset{project.datasetCount === 1 ? "" : "s"}
    </Chip>
  );

  return (
    <ExploreDatasetCard
      dataset={toItem(project)}
      headerBadges={countBadge}
      ownerOverlay={overlay}
      onCardClick={() => router.push(`/account/projects/${project.id}`)}
    />
  );
}
