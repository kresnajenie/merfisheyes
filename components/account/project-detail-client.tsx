"use client";

import type { MetadataDraft } from "./metadata-edit-modal";
import type { SubmissionStatus } from "./types";

import { useEffect, useMemo, useRef, useState } from "react";
import { Card, CardBody } from "@heroui/card";
import { Chip } from "@heroui/chip";
import { Button } from "@heroui/button";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import NextLink from "next/link";
import { useRouter } from "next/navigation";
import { toast } from "react-toastify";

import { MetadataEditModal } from "./metadata-edit-modal";
import { SubmissionBadge } from "./submission-badge";

interface MemberDataset {
  id: string;
  title: string | null;
  datasetType: string | null;
  thumbnailUrl: string | null;
  status: string;
  numCells: number;
  numGenes: number;
  ingestSource: string | null;
  s3BaseUrl: string | null;
}

interface ProjectData {
  id: string;
  title: string;
  description: string | null;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  tags: string[];
  thumbnailUrl: string | null;
  externalLink: string | null;
  publicationLink: string | null;
  metadata: Record<string, string> | null;
  datasets: { sortOrder: number; dataset: MemberDataset }[];
}

interface PickerDataset {
  id: string;
  title: string | null;
  datasetType: string | null;
  thumbnailUrl: string | null;
}

function viewerUrl(d: MemberDataset): string | null {
  if (d.status !== "COMPLETE") return null;
  const base = d.datasetType === "single_molecule" ? "sm-viewer" : "viewer";

  if (d.ingestSource === "s3_registered" && d.s3BaseUrl) {
    return `/${base}/from-s3?url=${encodeURIComponent(d.s3BaseUrl)}`;
  }

  return `/${base}/${d.id}`;
}

export function ProjectDetailClient({
  initialProject,
  ownDatasets,
  submission,
}: {
  initialProject: ProjectData;
  ownDatasets: PickerDataset[];
  submission: SubmissionStatus | null;
}) {
  const router = useRouter();
  const [project, setProject] = useState(initialProject);
  const [editOpen, setEditOpen] = useState(false);
  const [pickerOpen, setPickerOpen] = useState(false);
  const [thumbOpen, setThumbOpen] = useState(false);
  const [uploadingThumb, setUploadingThumb] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Keep local state in sync when router.refresh() re-runs the server page.
  useEffect(() => setProject(initialProject), [initialProject]);

  const members = useMemo(
    () => project.datasets.map((pd) => pd.dataset),
    [project],
  );
  const memberIds = useMemo(() => new Set(members.map((m) => m.id)), [members]);
  const addable = ownDatasets.filter((d) => !memberIds.has(d.id));
  // Thumbnail-reuse options: the project's own member datasets that have one.
  const datasetsWithThumb = members.filter((d) => d.thumbnailUrl);

  // ---- Thumbnail: upload a file, reuse a dataset's, or clear ----
  const MAX_THUMB = 6 * 1024 * 1024;

  const uploadThumbnail = async (file: File) => {
    if (file.size > MAX_THUMB) {
      toast.error("Image must be 6 MB or smaller.");

      return;
    }
    setUploadingThumb(true);
    try {
      const res = await fetch(`/api/projects/${project.id}/thumbnail`, {
        method: "POST",
        headers: { "Content-Type": file.type || "image/jpeg" },
        body: file,
      });

      if (res.ok) {
        toast.success("Thumbnail updated.");
        setThumbOpen(false);
        router.refresh();
      } else {
        toast.error("Couldn't upload the image.");
      }
    } finally {
      setUploadingThumb(false);
    }
  };

  const setThumbnailUrl = async (url: string | null) => {
    const res = await fetch(`/api/projects/${project.id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ thumbnailUrl: url }),
    });

    if (res.ok) {
      setThumbOpen(false);
      router.refresh();
    } else {
      toast.error("Couldn't update the thumbnail.");
    }
  };

  const saveMetadata = async (draft: MetadataDraft): Promise<boolean> => {
    const res = await fetch(`/api/projects/${project.id}`, {
      method: "PATCH",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(draft),
    });

    if (res.ok) {
      router.refresh();

      return true;
    }

    return false;
  };

  const addDataset = async (datasetId: string) => {
    const res = await fetch(`/api/projects/${project.id}/datasets`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ datasetId }),
    });

    if (res.ok) router.refresh();
    else toast.error("Couldn't add dataset.");
  };

  const removeDataset = async (datasetId: string) => {
    const res = await fetch(
      `/api/projects/${project.id}/datasets/${datasetId}`,
      { method: "DELETE" },
    );

    if (res.ok) router.refresh();
    else toast.error("Couldn't remove dataset.");
  };

  const submit = async () => {
    const res = await fetch("/api/community/submit", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type: "project", id: project.id }),
    });

    if (res.ok) {
      const j = await res.json();

      toast.success(
        j.reviewStatus === "APPROVED"
          ? "Submission updated — live on Explore."
          : "Submitted for review.",
      );
      router.refresh();
    } else {
      const j = await res.json().catch(() => ({}));

      toast.error(j.error ?? "Couldn't submit.");
    }
  };

  const withdraw = async () => {
    if (!confirm("Remove this project from Explore?")) return;
    const res = await fetch("/api/community/submit", {
      method: "DELETE",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ type: "project", id: project.id }),
    });

    if (res.ok) {
      toast.success("Withdrawn.");
      router.refresh();
    } else {
      toast.error("Couldn't withdraw.");
    }
  };

  return (
    <div className="w-full max-w-6xl mx-auto py-8 px-4">
      <NextLink
        className="text-sm text-default-500 hover:text-default-700"
        href="/account"
      >
        ← Back to your datasets
      </NextLink>

      <div className="flex flex-col sm:flex-row sm:items-start gap-4 mt-4 mb-8">
        <div className="w-full sm:w-56 shrink-0">
          {project.thumbnailUrl ? (
            // eslint-disable-next-line @next/next/no-img-element
            <img
              alt={project.title}
              className="w-full h-40 object-cover rounded-xl"
              src={project.thumbnailUrl}
            />
          ) : (
            <div className="flex h-40 w-full items-center justify-center rounded-xl bg-default-100 text-sm text-default-400">
              No thumbnail
            </div>
          )}
          <Button
            className="mt-2 w-full"
            size="sm"
            variant="flat"
            onPress={() => setThumbOpen(true)}
          >
            {project.thumbnailUrl ? "Change thumbnail" : "Add thumbnail"}
          </Button>
        </div>
        <div className="flex-1 flex flex-col gap-3">
          <div className="flex items-start justify-between gap-3">
            <div>
              <div className="flex items-center gap-2 flex-wrap">
                <h1 className="text-2xl font-bold">{project.title}</h1>
                <Chip size="sm" variant="flat">
                  Project
                </Chip>
                <SubmissionBadge submission={submission} />
              </div>
              {project.description && (
                <p className="text-default-500 mt-2 max-w-2xl">
                  {project.description}
                </p>
              )}
            </div>
            <div className="flex gap-2 shrink-0">
              <Button
                size="sm"
                variant="flat"
                onPress={() => setEditOpen(true)}
              >
                Edit
              </Button>
              {submission ? (
                <>
                  <Button size="sm" variant="flat" onPress={submit}>
                    Update
                  </Button>
                  <Button
                    color="danger"
                    size="sm"
                    variant="light"
                    onPress={withdraw}
                  >
                    Withdraw
                  </Button>
                </>
              ) : (
                <Button
                  color="primary"
                  isDisabled={members.length === 0}
                  size="sm"
                  onPress={submit}
                >
                  Submit to Explore
                </Button>
              )}
            </div>
          </div>

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
            {project.platform && (
              <Chip size="sm" variant="flat">
                {project.platform}
              </Chip>
            )}
            {project.tags.map((t) => (
              <Chip key={t} className="text-[10px]" size="sm" variant="dot">
                {t}
              </Chip>
            ))}
          </div>
        </div>
      </div>

      <div className="flex items-center justify-between mb-4">
        <h2 className="text-lg font-semibold">Datasets ({members.length})</h2>
        <Button
          color="primary"
          isDisabled={addable.length === 0}
          size="sm"
          variant="flat"
          onPress={() => setPickerOpen(true)}
        >
          Add datasets
        </Button>
      </div>

      {members.length === 0 ? (
        <p className="text-default-500 py-10 text-center">
          No datasets in this project yet. Click &quot;Add datasets&quot; to add
          some.
        </p>
      ) : (
        <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
          {members.map((d) => {
            const url = viewerUrl(d);

            return (
              <Card
                key={d.id}
                isBlurred
                className="border-none bg-background/60 dark:bg-default-100/50"
                shadow="sm"
              >
                <CardBody className="p-0 overflow-hidden">
                  {d.thumbnailUrl && (
                    <div className="h-28 w-full bg-default-200/30">
                      {/* eslint-disable-next-line @next/next/no-img-element */}
                      <img
                        alt={d.title ?? "dataset"}
                        className="w-full h-full object-cover"
                        src={d.thumbnailUrl}
                      />
                    </div>
                  )}
                  <div className="p-3 flex flex-col gap-2">
                    <div className="flex items-center gap-1">
                      <Chip
                        color={
                          d.datasetType === "single_molecule"
                            ? "secondary"
                            : "primary"
                        }
                        size="sm"
                        variant="solid"
                      >
                        {d.datasetType === "single_molecule" ? "SM" : "SC"}
                      </Chip>
                    </div>
                    <p className="text-sm font-medium line-clamp-2">
                      {d.title || "Untitled dataset"}
                    </p>
                    <div className="flex gap-2">
                      {url && (
                        <Button
                          as={NextLink}
                          color="primary"
                          href={url}
                          size="sm"
                          variant="flat"
                        >
                          Open
                        </Button>
                      )}
                      <Button
                        color="danger"
                        size="sm"
                        variant="light"
                        onPress={() => removeDataset(d.id)}
                      >
                        Remove
                      </Button>
                    </div>
                  </div>
                </CardBody>
              </Card>
            );
          })}
        </div>
      )}

      <MetadataEditModal
        showThumbnail
        heading="Edit project"
        initial={{ ...project, metadata: project.metadata ?? {} }}
        isOpen={editOpen}
        onClose={() => setEditOpen(false)}
        onSave={saveMetadata}
      />

      <Modal
        isOpen={pickerOpen}
        scrollBehavior="inside"
        size="md"
        onClose={() => setPickerOpen(false)}
      >
        <ModalContent>
          <ModalHeader>Add datasets to project</ModalHeader>
          <ModalBody className="flex flex-col gap-1">
            {addable.length === 0 ? (
              <p className="text-sm text-default-500 py-4">
                All your datasets are already in this project.
              </p>
            ) : (
              addable.map((d) => (
                <div
                  key={d.id}
                  className="flex items-center justify-between rounded-lg px-2 py-1.5 hover:bg-default-100"
                >
                  <div className="flex items-center gap-2 min-w-0">
                    <Chip
                      color={
                        d.datasetType === "single_molecule"
                          ? "secondary"
                          : "primary"
                      }
                      size="sm"
                      variant="flat"
                    >
                      {d.datasetType === "single_molecule" ? "SM" : "SC"}
                    </Chip>
                    <span className="text-sm truncate">
                      {d.title || "Untitled dataset"}
                    </span>
                  </div>
                  <Button
                    color="primary"
                    size="sm"
                    variant="flat"
                    onPress={() => addDataset(d.id)}
                  >
                    Add
                  </Button>
                </div>
              ))
            )}
          </ModalBody>
          <ModalFooter>
            <Button variant="light" onPress={() => setPickerOpen(false)}>
              Done
            </Button>
          </ModalFooter>
        </ModalContent>
      </Modal>

      <Modal
        isOpen={thumbOpen}
        scrollBehavior="inside"
        size="lg"
        onClose={() => setThumbOpen(false)}
      >
        <ModalContent>
          <ModalHeader>Project thumbnail</ModalHeader>
          <ModalBody className="flex flex-col gap-5">
            {/* Upload a file */}
            <div>
              <p className="mb-1 text-sm font-medium">Upload an image</p>
              <input
                ref={fileInputRef}
                accept="image/*"
                className="hidden"
                type="file"
                onChange={(e) => {
                  const f = e.target.files?.[0];

                  if (f) uploadThumbnail(f);
                  e.target.value = "";
                }}
              />
              <Button
                isLoading={uploadingThumb}
                size="sm"
                variant="flat"
                onPress={() => fileInputRef.current?.click()}
              >
                Choose image…
              </Button>
              <p className="mt-1 text-xs text-default-400">
                JPG or PNG, up to 6 MB.
              </p>
            </div>

            {/* Reuse one of your datasets' thumbnails */}
            <div className="border-t border-default-200 pt-4">
              <p className="mb-2 text-sm font-medium">
                Use a thumbnail from a dataset in this project
              </p>
              {datasetsWithThumb.length === 0 ? (
                <p className="text-xs text-default-400">
                  None of the datasets in this project have a thumbnail yet.
                </p>
              ) : (
                <div className="grid grid-cols-3 gap-2">
                  {datasetsWithThumb.map((d) => (
                    <button
                      key={d.id}
                      className="overflow-hidden rounded-lg border border-default-200 text-left transition-colors hover:border-primary"
                      type="button"
                      onClick={() => setThumbnailUrl(d.thumbnailUrl!)}
                    >
                      {/* eslint-disable-next-line @next/next/no-img-element */}
                      <img
                        alt={d.title ?? "dataset thumbnail"}
                        className="h-20 w-full object-cover"
                        src={d.thumbnailUrl!}
                      />
                      <span className="block truncate px-1 py-0.5 text-[10px] text-default-500">
                        {d.title || "Untitled"}
                      </span>
                    </button>
                  ))}
                </div>
              )}
            </div>
          </ModalBody>
          <ModalFooter>
            {project.thumbnailUrl && (
              <Button
                color="danger"
                size="sm"
                variant="light"
                onPress={() => setThumbnailUrl(null)}
              >
                Remove
              </Button>
            )}
            <Button variant="light" onPress={() => setThumbOpen(false)}>
              Done
            </Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </div>
  );
}
