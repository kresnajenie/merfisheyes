"use client";

import type { MetadataDraft } from "./metadata-edit-modal";
import type { SubmissionStatus } from "./types";
import type {
  CatalogDatasetItem,
  CatalogDatasetEntry,
} from "@/components/explore/types";

import { useEffect, useMemo, useRef, useState } from "react";
import { Chip } from "@heroui/chip";
import { Button } from "@heroui/button";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { useRouter } from "next/navigation";
import { toast } from "react-toastify";

import { MetadataEditModal } from "./metadata-edit-modal";
import { SubmissionBadge } from "./submission-badge";

import { DatasetDetailPage } from "@/components/explore/dataset-detail-page";

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
  genes?: string[];
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

/** Map the project + its members onto the shared Explore detail shape. */
function toItem(project: ProjectData, members: MemberDataset[]): CatalogDatasetItem {
  return {
    id: project.id,
    title: project.title,
    description: project.description,
    species: project.species,
    disease: project.disease,
    institute: project.institute,
    tissue: project.tissue,
    platform: project.platform,
    tags: project.tags ?? [],
    // Union of the member datasets' indexed genes — same list the Explore
    // detail page shows once the project is published.
    genes: [...new Set(members.flatMap((d) => d.genes ?? []))].sort(),
    thumbnailUrl: project.thumbnailUrl,
    bilCode: null,
    metadata: (project.metadata ?? {}) as Record<string, unknown>,
    externalLink: project.externalLink,
    publicationLink: project.publicationLink,
    entries: members.map((d, i) => ({
      id: d.id,
      label: d.title ?? "Untitled dataset",
      // Legacy rows store raw formats ("h5ad", "merscope") here; the app-wide
      // convention is anything not single_molecule renders as single cell.
      datasetType:
        d.datasetType === "single_molecule" ? "single_molecule" : "single_cell",
      // s3-registered datasets open via from-s3; everything else by id. A
      // non-COMPLETE member gets no link (entry renders inert).
      s3BaseUrl:
        d.status === "COMPLETE" &&
        d.ingestSource === "s3_registered" &&
        d.s3BaseUrl
          ? d.s3BaseUrl
          : null,
      datasetId: d.status === "COMPLETE" ? d.id : null,
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
  // Members are kept as their own state so add/remove can flip optimistically.
  const [members, setMembers] = useState<MemberDataset[]>(() =>
    initialProject.datasets.map((pd) => pd.dataset),
  );
  const [editOpen, setEditOpen] = useState(false);
  const [pickerOpen, setPickerOpen] = useState(false);
  const [thumbOpen, setThumbOpen] = useState(false);
  const [uploadingThumb, setUploadingThumb] = useState(false);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Keep local state in sync when router.refresh() re-runs the server page.
  useEffect(() => {
    setProject(initialProject);
    setMembers(initialProject.datasets.map((pd) => pd.dataset));
  }, [initialProject]);

  const memberIds = useMemo(() => new Set(members.map((m) => m.id)), [members]);
  const addable = ownDatasets.filter((d) => !memberIds.has(d.id));
  // Thumbnail-reuse options: the project's own member datasets that have one.
  const datasetsWithThumb = members.filter((d) => d.thumbnailUrl);

  const item = useMemo(() => toItem(project, members), [project, members]);

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

  // ---- Membership (optimistic: flip local state, settle in background) ----
  const addDataset = async (d: PickerDataset) => {
    setMembers((ms) =>
      ms.some((m) => m.id === d.id)
        ? ms
        : [
            ...ms,
            {
              ...d,
              // Picker only offers COMPLETE datasets; the refresh after the
              // API settles fills in the real ingest fields.
              status: "COMPLETE",
              numCells: 0,
              numGenes: 0,
              ingestSource: null,
              s3BaseUrl: null,
            },
          ],
    );

    const res = await fetch(`/api/projects/${project.id}/datasets`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ datasetId: d.id }),
    }).catch(() => null);

    if (res?.ok) {
      router.refresh();
    } else {
      setMembers((ms) => ms.filter((m) => m.id !== d.id));
      toast.error("Couldn't add dataset.");
    }
  };

  const removeDataset = async (datasetId: string) => {
    const removed = members.find((m) => m.id === datasetId);

    setMembers((ms) => ms.filter((m) => m.id !== datasetId));

    const res = await fetch(
      `/api/projects/${project.id}/datasets/${datasetId}`,
      { method: "DELETE" },
    ).catch(() => null);

    if (res?.ok) {
      router.refresh();
    } else {
      if (removed) setMembers((ms) => [...ms, removed]);
      toast.error("Couldn't remove dataset.");
    }
  };

  // ---- Submit / withdraw ----
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

  const headerActions = (
    <>
      <SubmissionBadge submission={submission} />
      <Button size="sm" variant="flat" onPress={() => setEditOpen(true)}>
        Edit
      </Button>
      <Button size="sm" variant="flat" onPress={() => setThumbOpen(true)}>
        {project.thumbnailUrl ? "Change thumbnail" : "Add thumbnail"}
      </Button>
      <Button
        color="primary"
        isDisabled={addable.length === 0}
        size="sm"
        variant="flat"
        onPress={() => setPickerOpen(true)}
      >
        Add datasets
      </Button>
      {submission ? (
        <>
          <Button size="sm" variant="flat" onPress={submit}>
            Update
          </Button>
          <Button color="danger" size="sm" variant="light" onPress={withdraw}>
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
    </>
  );

  const entryActions = (entry: CatalogDatasetEntry) => (
    <Button
      isIconOnly
      aria-label="Remove from project"
      className="min-w-7 w-7 h-7 bg-default-800/70 text-white backdrop-blur"
      size="sm"
      title="Remove from project"
      variant="flat"
      onPress={() => removeDataset(entry.id)}
    >
      ✕
    </Button>
  );

  return (
    <div className="w-full max-w-7xl mx-auto py-4 px-4">
      {members.length === 0 && (
        <p className="mb-2 rounded-xl bg-default-100 px-4 py-3 text-sm text-default-500">
          No datasets in this project yet. Click &quot;Add datasets&quot; to add
          some.
        </p>
      )}

      <DatasetDetailPage
        backLabel="← Back to your projects"
        dataset={item}
        entryActions={entryActions}
        headerActions={headerActions}
        onBack={() => router.push("/account?view=projects")}
      />

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
        size="lg"
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
                    onPress={() => addDataset(d)}
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
