"use client";

import type { MetadataDraft } from "./metadata-edit-modal";
import type { SubmissionStatus, MetadataBag } from "./types";
import type { CatalogDatasetItem } from "@/components/explore/types";

import { useState } from "react";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";
import { useRouter } from "next/navigation";
import { toast } from "react-toastify";

import { MetadataEditModal } from "./metadata-edit-modal";
import { SubmissionBadge } from "./submission-badge";

import { DatasetDetailPage } from "@/components/explore/dataset-detail-page";
import { OverlayManager } from "@/components/overlay-manager";

interface DatasetData {
  id: string;
  title: string | null;
  status: string;
  datasetType: string | null;
  s3BaseUrl: string | null;
  ingestSource: string | null;
  thumbnailUrl: string | null;
  numCells: number | null;
  numGenes: number | null;
  description: string | null;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  tags: string[];
  externalLink: string | null;
  publicationLink: string | null;
  metadata: MetadataBag | null;
}

interface Props {
  initialDataset: DatasetData;
  submission: SubmissionStatus | null;
  viewerUrl: string | null;
}

export function AccountDatasetDetailClient({
  initialDataset,
  submission,
  viewerUrl,
}: Props) {
  const router = useRouter();
  const [dataset, setDataset] = useState(initialDataset);
  const [editOpen, setEditOpen] = useState(false);
  const [submitting, setSubmitting] = useState(false);

  // Map the owned dataset onto the shared Explore detail shape. A single
  // synthetic entry drives the "Open in viewer" card; app-uploaded datasets
  // load by id, S3-registered ones by URL.
  const item: CatalogDatasetItem = {
    id: dataset.id,
    title: dataset.title ?? "Untitled dataset",
    description: dataset.description,
    species: dataset.species,
    disease: dataset.disease,
    institute: dataset.institute,
    tissue: dataset.tissue,
    platform: dataset.platform,
    tags: dataset.tags ?? [],
    genes: [],
    thumbnailUrl: dataset.thumbnailUrl,
    bilCode: null,
    metadata: (dataset.metadata ?? {}) as Record<string, unknown>,
    externalLink: dataset.externalLink,
    publicationLink: dataset.publicationLink,
    entries: viewerUrl
      ? [
          {
            id: dataset.id,
            label: dataset.title ?? "Open viewer",
            datasetType: dataset.datasetType ?? "single_cell",
            s3BaseUrl:
              dataset.ingestSource === "s3_registered"
                ? dataset.s3BaseUrl
                : null,
            datasetId:
              dataset.ingestSource === "s3_registered" ? null : dataset.id,
            thumbnailUrl: dataset.thumbnailUrl,
            sortOrder: 0,
          },
        ]
      : [],
    isPublished: false,
    isFeatured: false,
    isBil: false,
    isInternal: false,
    isCommunity: false,
    sortOrder: 0,
    numCells: dataset.numCells,
    numGenes: dataset.numGenes,
  };

  const saveMetadata = async (draft: MetadataDraft): Promise<boolean> => {
    const res = await fetch(`/api/ingest/mine/${dataset.id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        title: draft.title,
        description: draft.description,
        species: draft.species,
        disease: draft.disease,
        tissue: draft.tissue,
        platform: draft.platform,
        institute: draft.institute,
        tags: draft.tags,
        externalLink: draft.externalLink,
        publicationLink: draft.publicationLink,
        metadata: draft.metadata,
      }),
    });

    if (!res.ok) return false;

    // Reflect the edit locally (title/metadata/etc.) without a full reload.
    setDataset((prev) => ({
      ...prev,
      title: draft.title,
      description: draft.description,
      species: draft.species,
      disease: draft.disease,
      tissue: draft.tissue,
      platform: draft.platform,
      institute: draft.institute,
      tags: draft.tags,
      externalLink: draft.externalLink,
      publicationLink: draft.publicationLink,
      metadata: draft.metadata,
    }));
    router.refresh();

    return true;
  };

  const submit = async () => {
    setSubmitting(true);
    try {
      const res = await fetch("/api/community/submit", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ type: "dataset", id: dataset.id }),
      });

      if (res.ok) {
        const j = await res.json();

        toast.success(
          j.reviewStatus === "APPROVED"
            ? "Submission updated — live on Explore."
            : "Submitted for review. An admin will approve it soon.",
        );
        router.refresh();
      } else {
        const j = await res.json().catch(() => ({}));

        toast.error(j.error ?? "Couldn't submit.");
      }
    } finally {
      setSubmitting(false);
    }
  };

  const withdraw = async () => {
    if (!confirm("Remove this from Explore / cancel the submission?")) return;
    setSubmitting(true);
    try {
      const res = await fetch("/api/community/submit", {
        method: "DELETE",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ type: "dataset", id: dataset.id }),
      });

      if (res.ok) {
        toast.success("Withdrawn.");
        router.refresh();
      } else {
        toast.error("Couldn't withdraw.");
      }
    } finally {
      setSubmitting(false);
    }
  };

  const headerActions = (
    <>
      {submission && <SubmissionBadge submission={submission} />}
      {dataset.status !== "COMPLETE" && (
        <Chip color="warning" size="sm" variant="flat">
          {dataset.status.toLowerCase()}
        </Chip>
      )}
      <Button size="sm" variant="flat" onPress={() => setEditOpen(true)}>
        Edit details
      </Button>
      {submission ? (
        <Button
          color="danger"
          isLoading={submitting}
          size="sm"
          variant="light"
          onPress={withdraw}
        >
          Withdraw
        </Button>
      ) : (
        <Button
          color="primary"
          isDisabled={dataset.status !== "COMPLETE"}
          isLoading={submitting}
          size="sm"
          onPress={submit}
        >
          Submit to Explore
        </Button>
      )}
    </>
  );

  return (
    <div className="w-full max-w-7xl mx-auto py-4 px-4">
      <DatasetDetailPage
        backLabel="← Back to your datasets"
        dataset={item}
        headerActions={headerActions}
        onBack={() => router.push("/account?view=datasets")}
      />

      {dataset.datasetType === "single_cell" &&
        dataset.status === "COMPLETE" && (
          <div className="mt-6">
            <OverlayManager scDatasetId={dataset.id} />
          </div>
        )}

      <MetadataEditModal
        heading="Edit dataset"
        initial={{
          title: dataset.title ?? "",
          description: dataset.description,
          species: dataset.species,
          disease: dataset.disease,
          tissue: dataset.tissue,
          platform: dataset.platform,
          institute: dataset.institute,
          tags: dataset.tags ?? [],
          externalLink: dataset.externalLink,
          publicationLink: dataset.publicationLink,
          metadata: dataset.metadata ?? {},
        }}
        isOpen={editOpen}
        onClose={() => setEditOpen(false)}
        onSave={saveMetadata}
      />
    </div>
  );
}
