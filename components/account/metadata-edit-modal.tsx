"use client";

import type { EditableMetadata } from "./types";

import { useEffect, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Input, Textarea } from "@heroui/input";
import { Button } from "@heroui/button";

export interface MetadataDraft extends EditableMetadata {
  title: string;
  /** Only editable for projects (datasets set their thumbnail in the viewer). */
  thumbnailUrl: string | null;
}

interface MetadataEditModalProps {
  isOpen: boolean;
  onClose: () => void;
  heading: string;
  /** Initial values to seed the form. */
  initial: Partial<MetadataDraft>;
  /** When true, show the Thumbnail URL field (projects). */
  showThumbnail?: boolean;
  /**
   * Persist the draft. Return true on success (modal closes), false to keep it
   * open. The parent builds the request from the returned payload.
   */
  onSave: (draft: MetadataDraft) => Promise<boolean>;
}

function toForm(initial: Partial<MetadataDraft>) {
  return {
    title: initial.title ?? "",
    thumbnailUrl: initial.thumbnailUrl ?? "",
    description: initial.description ?? "",
    species: initial.species ?? "",
    disease: initial.disease ?? "",
    tissue: initial.tissue ?? "",
    platform: initial.platform ?? "",
    institute: initial.institute ?? "",
    tags: (initial.tags ?? []).join(", "),
    externalLink: initial.externalLink ?? "",
    publicationLink: initial.publicationLink ?? "",
  };
}

export function MetadataEditModal({
  isOpen,
  onClose,
  heading,
  initial,
  showThumbnail = false,
  onSave,
}: MetadataEditModalProps) {
  const [form, setForm] = useState(() => toForm(initial));
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Reseed whenever the modal reopens for a (possibly different) target.
  useEffect(() => {
    if (isOpen) {
      setForm(toForm(initial));
      setError(null);
    }
  }, [isOpen]);

  const set = (key: keyof ReturnType<typeof toForm>) => (val: string) =>
    setForm((prev) => ({ ...prev, [key]: val }));

  const handleSave = async () => {
    if (!form.title.trim()) {
      setError("Title is required");

      return;
    }
    setSaving(true);
    setError(null);

    const draft: MetadataDraft = {
      title: form.title.trim(),
      thumbnailUrl: form.thumbnailUrl.trim() || null,
      description: form.description.trim() || null,
      species: form.species.trim() || null,
      disease: form.disease.trim() || null,
      tissue: form.tissue.trim() || null,
      platform: form.platform.trim() || null,
      institute: form.institute.trim() || null,
      tags: form.tags
        .split(",")
        .map((t) => t.trim())
        .filter(Boolean),
      externalLink: form.externalLink.trim() || null,
      publicationLink: form.publicationLink.trim() || null,
    };

    const ok = await onSave(draft);

    setSaving(false);
    if (ok) onClose();
    else setError("Couldn't save. Please try again.");
  };

  return (
    <Modal isOpen={isOpen} scrollBehavior="inside" size="2xl" onClose={onClose}>
      <ModalContent>
        <ModalHeader>{heading}</ModalHeader>
        <ModalBody className="flex flex-col gap-4">
          {error && (
            <div className="text-danger text-sm bg-danger-50 dark:bg-danger-100/10 px-4 py-2 rounded-lg">
              {error}
            </div>
          )}

          <Input
            isRequired
            label="Title"
            value={form.title}
            onValueChange={set("title")}
          />

          <Textarea
            label="Description"
            minRows={3}
            value={form.description}
            onValueChange={set("description")}
          />

          {showThumbnail && (
            <Input
              description="Direct URL to a thumbnail image"
              label="Thumbnail URL"
              value={form.thumbnailUrl}
              onValueChange={set("thumbnailUrl")}
            />
          )}

          <div className="grid grid-cols-2 gap-4">
            <Input
              label="Species"
              value={form.species}
              onValueChange={set("species")}
            />
            <Input
              label="Disease"
              value={form.disease}
              onValueChange={set("disease")}
            />
            <Input
              label="Tissue"
              value={form.tissue}
              onValueChange={set("tissue")}
            />
            <Input
              label="Platform"
              value={form.platform}
              onValueChange={set("platform")}
            />
            <Input
              label="Institute"
              value={form.institute}
              onValueChange={set("institute")}
            />
            <Input
              description="Comma-separated"
              label="Tags"
              value={form.tags}
              onValueChange={set("tags")}
            />
          </div>

          <Input
            description="Original source page (publication, data repository)"
            label="Source Link"
            value={form.externalLink}
            onValueChange={set("externalLink")}
          />
          <Input
            description="Link to the publication (DOI, PubMed, journal page)"
            label="Publication Link"
            value={form.publicationLink}
            onValueChange={set("publicationLink")}
          />
        </ModalBody>
        <ModalFooter>
          <Button variant="light" onPress={onClose}>
            Cancel
          </Button>
          <Button color="primary" isLoading={saving} onPress={handleSave}>
            Save
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
