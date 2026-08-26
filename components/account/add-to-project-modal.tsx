"use client";

import type { DatasetRow, ProjectRow } from "./types";

import { useEffect, useMemo, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";
import { Input } from "@heroui/input";
import { Chip } from "@heroui/chip";
import { Pagination } from "@heroui/pagination";

const PAGE_SIZE = 12;

interface AddToProjectModalProps {
  isOpen: boolean;
  onClose: () => void;
  dataset: DatasetRow | null;
  projects: ProjectRow[];
  /** Toggle membership; returns true on success. */
  onToggle: (
    projectId: string,
    datasetId: string,
    member: boolean,
  ) => Promise<boolean>;
  /** Create a new project (title only) and add the dataset to it. */
  onCreateAndAdd: (title: string, datasetId: string) => Promise<boolean>;
}

/** One selectable project card: thumbnail, title, count, membership state. */
function ProjectPickCard({
  project,
  selected,
  onToggle,
}: {
  project: ProjectRow;
  selected: boolean;
  onToggle: () => void;
}) {
  const thumb =
    project.thumbnailUrl ??
    project.datasets.find((d) => d.thumbnailUrl)?.thumbnailUrl ??
    null;

  return (
    <button
      className={`group relative overflow-hidden rounded-xl border text-left transition-all ${
        selected
          ? "border-primary ring-2 ring-primary/40"
          : "border-default-200 hover:border-default-400"
      }`}
      type="button"
      onClick={onToggle}
    >
      {/* Thumbnail / placeholder */}
      <div className="relative h-24 w-full bg-default-200/30">
        <div className="h-full w-full bg-gradient-to-br from-default-100 to-default-200/40" />
        {thumb && (
          <img
            alt={project.title}
            className="absolute inset-0 h-full w-full object-cover"
            src={thumb}
            onError={(e) => {
              e.currentTarget.style.display = "none";
            }}
          />
        )}
        {/* Membership check */}
        <div
          className={`absolute top-2 right-2 flex h-6 w-6 items-center justify-center rounded-full transition-colors ${
            selected
              ? "bg-primary text-white"
              : "bg-default-800/50 text-white/60 group-hover:text-white"
          }`}
        >
          <svg
            className="h-3.5 w-3.5"
            fill="none"
            stroke="currentColor"
            strokeWidth={3}
            viewBox="0 0 24 24"
          >
            <path d="M5 13l4 4L19 7" strokeLinecap="round" />
          </svg>
        </div>
      </div>
      <div className="flex flex-col gap-1 p-2.5">
        <span className="truncate text-sm font-medium" title={project.title}>
          {project.title}
        </span>
        <Chip className="text-[10px]" size="sm" variant="flat">
          {project.datasetCount} dataset{project.datasetCount === 1 ? "" : "s"}
        </Chip>
      </div>
    </button>
  );
}

export function AddToProjectModal({
  isOpen,
  onClose,
  dataset,
  projects,
  onToggle,
  onCreateAndAdd,
}: AddToProjectModalProps) {
  const [newTitle, setNewTitle] = useState("");
  const [creating, setCreating] = useState(false);
  const [query, setQuery] = useState("");
  const [page, setPage] = useState(1);

  // Fresh search each time the modal opens for a dataset.
  useEffect(() => {
    if (isOpen) {
      setQuery("");
      setPage(1);
    }
  }, [isOpen]);

  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase();

    if (!q) return projects;

    return projects.filter((p) => p.title.toLowerCase().includes(q));
  }, [projects, query]);

  const pages = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE));
  const clampedPage = Math.min(page, pages);
  const visible = filtered.slice(
    (clampedPage - 1) * PAGE_SIZE,
    clampedPage * PAGE_SIZE,
  );

  if (!dataset) return null;

  const memberOf = (p: ProjectRow) =>
    p.datasets.some((d) => d.id === dataset.id);

  // Membership updates optimistically in the parent, so the card flips
  // instantly — no need to disable/spin while the request settles.
  const toggle = (p: ProjectRow) => {
    void onToggle(p.id, dataset.id, !memberOf(p));
  };

  const create = async () => {
    const title = newTitle.trim();

    if (!title) return;
    setCreating(true);
    const ok = await onCreateAndAdd(title, dataset.id);

    setCreating(false);
    if (ok) setNewTitle("");
  };

  return (
    <Modal isOpen={isOpen} scrollBehavior="inside" size="5xl" onClose={onClose}>
      <ModalContent className="max-h-[85vh]">
        <ModalHeader className="flex flex-col gap-1">
          <span>Add to project</span>
          <span className="text-xs font-normal text-default-400">
            {dataset.title || "Untitled dataset"}
          </span>
        </ModalHeader>
        <ModalBody className="flex flex-col gap-4">
          {/* Create a new project — stays on top */}
          <div className="flex items-end gap-2">
            <Input
              label="New project"
              placeholder="Project title"
              size="sm"
              value={newTitle}
              onKeyDown={(e) => {
                if (e.key === "Enter") create();
              }}
              onValueChange={setNewTitle}
            />
            <Button
              color="primary"
              isDisabled={!newTitle.trim()}
              isLoading={creating}
              size="sm"
              onPress={create}
            >
              Create
            </Button>
          </div>

          {projects.length === 0 ? (
            <p className="text-sm text-default-500">
              No projects yet — create one above.
            </p>
          ) : (
            <>
              <Input
                classNames={{ inputWrapper: "bg-default-100" }}
                placeholder={`Search ${projects.length} projects…`}
                value={query}
                onValueChange={(v) => {
                  setQuery(v);
                  setPage(1);
                }}
              />

              {visible.length === 0 ? (
                <p className="py-8 text-center text-sm text-default-500">
                  No projects match “{query}”.
                </p>
              ) : (
                <div className="grid grid-cols-2 gap-3 sm:grid-cols-3 md:grid-cols-4">
                  {visible.map((p) => (
                    <ProjectPickCard
                      key={p.id}
                      project={p}
                      selected={memberOf(p)}
                      onToggle={() => toggle(p)}
                    />
                  ))}
                </div>
              )}

              {pages > 1 && (
                <div className="flex justify-center">
                  <Pagination
                    showControls
                    page={clampedPage}
                    size="sm"
                    total={pages}
                    onChange={setPage}
                  />
                </div>
              )}
            </>
          )}
        </ModalBody>
        <ModalFooter>
          <Button variant="light" onPress={onClose}>
            Done
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
