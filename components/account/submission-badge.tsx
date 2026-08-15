"use client";

import type { SubmissionStatus } from "./types";

import { Chip } from "@heroui/chip";
import { Tooltip } from "@heroui/tooltip";
import NextLink from "next/link";

/**
 * Compact badge showing where a dataset/project stands in the Explore
 * publishing pipeline. Renders nothing when the item was never submitted.
 */
export function SubmissionBadge({
  submission,
}: {
  submission: SubmissionStatus | null;
}) {
  if (!submission) return null;

  if (submission.reviewStatus === "APPROVED" && submission.isPublished) {
    return (
      <Chip
        as={NextLink}
        className="cursor-pointer"
        color="success"
        href={`/explore/${submission.catalogId}`}
        size="sm"
        variant="flat"
      >
        Published on Explore
      </Chip>
    );
  }

  if (submission.reviewStatus === "REJECTED") {
    const chip = (
      <Chip color="danger" size="sm" variant="flat">
        Not approved
      </Chip>
    );

    return submission.reviewNote ? (
      <Tooltip content={submission.reviewNote}>{chip}</Tooltip>
    ) : (
      chip
    );
  }

  return (
    <Chip color="warning" size="sm" variant="flat">
      Pending review
    </Chip>
  );
}
