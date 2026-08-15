"use client";

import type { ProjectRow } from "./types";

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

interface AccountProjectCardProps {
  project: ProjectRow;
  onEdit: (project: ProjectRow) => void;
  onSubmit: (project: ProjectRow) => void;
  onWithdraw: (project: ProjectRow) => void;
  onDelete: (project: ProjectRow) => void;
}

export function AccountProjectCard({
  project,
  onEdit,
  onSubmit,
  onWithdraw,
  onDelete,
}: AccountProjectCardProps) {
  const sub = project.submission;
  const canSubmit = project.datasetCount > 0;
  const href = `/account/projects/${project.id}`;
  // Fall back to a member dataset's thumbnail when the project has none set.
  const headerThumb =
    project.thumbnailUrl ??
    project.datasets.find((d) => d.thumbnailUrl)?.thumbnailUrl ??
    null;

  return (
    <Card
      isBlurred
      className="border-none bg-background/60 dark:bg-default-100/50"
      shadow="sm"
    >
      <CardBody className="p-0 overflow-hidden">
        <div className="relative h-36 w-full bg-default-200/30">
          {headerThumb ? (
            // eslint-disable-next-line @next/next/no-img-element
            <img
              alt={project.title}
              className="w-full h-full object-cover"
              src={headerThumb}
            />
          ) : (
            <ThumbnailPlaceholder className="h-full" kind="project" />
          )}
          <div className="absolute top-2 left-2">
            <Chip className="bg-default-800/70 text-white" size="sm">
              Project · {project.datasetCount}
            </Chip>
          </div>
        </div>

        <div className="p-4 flex flex-col gap-2">
          <h3 className="font-semibold text-foreground text-sm line-clamp-2">
            {project.title}
          </h3>

          {project.description && (
            <p className="text-xs text-default-500 line-clamp-2">
              {project.description}
            </p>
          )}

          <div className="flex flex-wrap gap-1">
            {project.species && (
              <Chip size="sm" variant="flat">
                {project.species}
              </Chip>
            )}
            {project.tissue && (
              <Chip size="sm" variant="flat">
                {project.tissue}
              </Chip>
            )}
          </div>

          {project.thumbnailUrl && (
            <p className="text-xs text-default-500">
              {project.datasetCount} dataset
              {project.datasetCount === 1 ? "" : "s"}
            </p>
          )}

          <div className="flex flex-wrap items-center gap-1.5">
            <SubmissionBadge submission={sub} />
          </div>
        </div>
      </CardBody>

      <CardFooter className="gap-2 pt-0 flex-wrap">
        <Button
          as={NextLink}
          color="primary"
          href={href}
          size="sm"
          variant="flat"
        >
          Open
        </Button>
        <Button size="sm" variant="flat" onPress={() => onEdit(project)}>
          Edit
        </Button>
        <Dropdown>
          <DropdownTrigger>
            <Button isIconOnly size="sm" variant="light">
              <span className="text-lg leading-none">⋯</span>
            </Button>
          </DropdownTrigger>
          <DropdownMenu aria-label="Project actions">
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
      </CardFooter>
    </Card>
  );
}
