"use client";

import type { DatasetRow, ProjectRow } from "./types";

import { useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Checkbox } from "@heroui/checkbox";
import { Button } from "@heroui/button";
import { Input } from "@heroui/input";

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

  if (!dataset) return null;

  const memberOf = (p: ProjectRow) =>
    p.datasets.some((d) => d.id === dataset.id);

  // Membership updates optimistically in the parent, so the checkbox flips
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
    <Modal isOpen={isOpen} scrollBehavior="inside" size="md" onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          <span>Add to project</span>
          <span className="text-xs font-normal text-default-400">
            {dataset.title || "Untitled dataset"}
          </span>
        </ModalHeader>
        <ModalBody className="flex flex-col gap-3">
          {projects.length === 0 ? (
            <p className="text-sm text-default-500">
              No projects yet — create one below.
            </p>
          ) : (
            <div className="flex flex-col gap-2">
              {projects.map((p) => (
                <div
                  key={p.id}
                  className="flex items-center justify-between rounded-lg px-2 py-1.5 hover:bg-default-100"
                >
                  <Checkbox
                    isSelected={memberOf(p)}
                    onValueChange={() => toggle(p)}
                  >
                    <span className="text-sm">{p.title}</span>
                  </Checkbox>
                  <span className="text-xs text-default-400">
                    {p.datasetCount} dataset{p.datasetCount === 1 ? "" : "s"}
                  </span>
                </div>
              ))}
            </div>
          )}

          <div className="flex items-end gap-2 border-t border-default-200 pt-3">
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
