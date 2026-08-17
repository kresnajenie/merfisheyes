"use client";

import type { DatasetRow } from "./types";

import { Card, CardBody, CardFooter } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { Button } from "@heroui/button";
import {
  Dropdown,
  DropdownTrigger,
  DropdownMenu,
  DropdownItem,
} from "@heroui/dropdown";
import NextLink from "next/link";

import { SubmissionBadge } from "./submission-badge";
import { ThumbnailPlaceholder } from "./thumbnail-placeholder";

const statusColor: Record<
  string,
  "success" | "warning" | "danger" | "default"
> = {
  COMPLETE: "success",
  PROCESSING: "warning",
  QUEUED: "warning",
  FAILED: "danger",
};

function formatCount(n: number | null): string | null {
  if (n == null) return null;
  if (n >= 1_000_000) return `${(n / 1_000_000).toFixed(1)}M`;
  if (n >= 1_000) return `${(n / 1_000).toFixed(0)}K`;

  return String(n);
}

interface AccountDatasetCardProps {
  dataset: DatasetRow;
  projectNames: string[];
  onEdit: (dataset: DatasetRow) => void;
  onSubmit: (dataset: DatasetRow) => void;
  onWithdraw: (dataset: DatasetRow) => void;
  onAddToProject: (dataset: DatasetRow) => void;
}

export function AccountDatasetCard({
  dataset,
  projectNames,
  onEdit,
  onSubmit,
  onWithdraw,
  onAddToProject,
}: AccountDatasetCardProps) {
  const isMolecule = dataset.datasetType === "single_molecule";
  const typeLabel = isMolecule ? "SM" : "SC";
  const isComplete = dataset.status === "COMPLETE";
  const cells = formatCount(dataset.numCells);
  const genes = formatCount(dataset.numGenes);
  const sub = dataset.submission;
  const canSubmit = isComplete;

  return (
    <Card
      isBlurred
      className="border-none bg-background/60 dark:bg-default-100/50"
      shadow="sm"
    >
      <CardBody className="p-0 overflow-hidden">
        {/* Header: real thumbnail when present, else a tinted placeholder so
            every card reads as intentional (matching the Explore card look). */}
        <div className="relative h-36 w-full bg-default-200/30">
          {dataset.thumbnailUrl ? (
            // eslint-disable-next-line @next/next/no-img-element
            <img
              alt={dataset.title ?? "dataset"}
              className="w-full h-full object-cover"
              src={dataset.thumbnailUrl}
            />
          ) : (
            <ThumbnailPlaceholder
              className="h-full"
              kind={isMolecule ? "single_molecule" : "single_cell"}
            />
          )}
          <div className="absolute top-2 left-2 flex gap-1">
            <Chip
              color={isMolecule ? "secondary" : "primary"}
              size="sm"
              variant="solid"
            >
              {typeLabel}
            </Chip>
            {dataset.adminOwned && (
              <Chip color="secondary" size="sm" variant="flat">
                admin
              </Chip>
            )}
          </div>
        </div>

        <div className="p-4 flex flex-col gap-2">
          <h3 className="font-semibold text-foreground text-sm line-clamp-2">
            {dataset.title || "Untitled dataset"}
          </h3>

          {dataset.description && (
            <p className="text-xs text-default-500 line-clamp-2">
              {dataset.description}
            </p>
          )}

          <div className="flex flex-wrap gap-1">
            {dataset.species && (
              <Chip size="sm" variant="flat">
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

          <div className="flex gap-3 text-xs text-default-500">
            {cells && <span>{cells} cells/molecules</span>}
            {genes && <span>{genes} genes</span>}
            <span>{dataset.viewCount.toLocaleString()} views</span>
          </div>

          <div className="flex flex-wrap items-center gap-1.5">
            {!isComplete && (
              <Chip
                color={statusColor[dataset.status] ?? "default"}
                size="sm"
                variant="flat"
              >
                {dataset.status.toLowerCase()}
              </Chip>
            )}
            <SubmissionBadge submission={sub} />
          </div>

          {projectNames.length > 0 && (
            <div className="flex flex-wrap items-center gap-1">
              {projectNames.map((name) => (
                <Chip
                  key={name}
                  className="max-w-full"
                  color="default"
                  size="sm"
                  startContent={
                    <span className="text-xs leading-none pl-0.5">📁</span>
                  }
                  variant="flat"
                >
                  <span className="truncate">{name}</span>
                </Chip>
              ))}
            </div>
          )}
        </div>
      </CardBody>

      <CardFooter className="gap-2 pt-0 flex-wrap">
        {dataset.viewerUrl ? (
          <Button
            as={NextLink}
            color="primary"
            href={dataset.viewerUrl}
            size="sm"
            variant="flat"
          >
            Open
          </Button>
        ) : null}
        <Button size="sm" variant="flat" onPress={() => onEdit(dataset)}>
          Edit
        </Button>
        <Button
          size="sm"
          startContent={<span className="text-sm leading-none">＋</span>}
          variant="flat"
          onPress={() => onAddToProject(dataset)}
        >
          Project
        </Button>
        {canSubmit || sub ? (
          <Dropdown>
            <DropdownTrigger>
              <Button isIconOnly size="sm" variant="light">
                <span className="text-lg leading-none">⋯</span>
              </Button>
            </DropdownTrigger>
            <DropdownMenu aria-label="Dataset actions">
              {canSubmit && !sub ? (
                <DropdownItem key="submit" onPress={() => onSubmit(dataset)}>
                  Submit to Explore
                </DropdownItem>
              ) : null}
              {canSubmit && sub ? (
                <DropdownItem key="resubmit" onPress={() => onSubmit(dataset)}>
                  {sub.reviewStatus === "REJECTED"
                    ? "Resubmit to Explore"
                    : "Update submission"}
                </DropdownItem>
              ) : null}
              {sub ? (
                <DropdownItem
                  key="withdraw"
                  className="text-danger"
                  color="danger"
                  onPress={() => onWithdraw(dataset)}
                >
                  Withdraw from Explore
                </DropdownItem>
              ) : null}
            </DropdownMenu>
          </Dropdown>
        ) : null}
      </CardFooter>
    </Card>
  );
}
