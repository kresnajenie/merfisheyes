"use client";

import { useEffect, useState } from "react";
import {
  Modal,
  ModalContent,
  ModalHeader,
  ModalBody,
  ModalFooter,
} from "@heroui/modal";
import { Button } from "@heroui/button";
import { Checkbox } from "@heroui/checkbox";
import { Input } from "@heroui/input";
import { toast } from "react-toastify";

const CONFIRM_WORD = "delete";

interface DeleteDatasetModalProps {
  /** The dataset to delete; null keeps the modal closed. */
  dataset: {
    id: string;
    title: string | null;
    datasetType: string | null;
    onExplore: boolean;
  } | null;
  onClose: () => void;
  onDeleted: (deletedIds: string[]) => void;
}

/**
 * Confirms and performs the permanent deletion of an owned dataset. A cell
 * dataset with a molecule overlay offers to delete that dataset along with it.
 */
export function DeleteDatasetModal({
  dataset,
  onClose,
  onDeleted,
}: DeleteDatasetModalProps) {
  const [typed, setTyped] = useState("");
  const [overlayId, setOverlayId] = useState<string | null>(null);
  const [withOverlay, setWithOverlay] = useState(false);
  const [deleting, setDeleting] = useState(false);

  useEffect(() => {
    setTyped("");
    setOverlayId(null);
    setWithOverlay(false);
    if (!dataset || dataset.datasetType === "single_molecule") return;

    let cancelled = false;

    fetch(`/api/datasets/${dataset.id}/overlay`)
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => !cancelled && setOverlayId(j?.smDatasetId ?? null))
      .catch(() => {});

    return () => {
      cancelled = true;
    };
  }, [dataset]);

  const confirm = async () => {
    if (!dataset) return;
    setDeleting(true);
    try {
      const res = await fetch(
        `/api/ingest/${dataset.id}${withOverlay ? "?withOverlay=1" : ""}`,
        { method: "DELETE" },
      );
      const body = await res.json().catch(() => ({}));

      if (res.ok) {
        toast.success(
          body.deleted?.length > 1 ? "Datasets deleted." : "Dataset deleted.",
        );
        onDeleted(body.deleted ?? [dataset.id]);
        onClose();
      } else {
        toast.error(body.message ?? "Couldn't delete the dataset.");
      }
    } finally {
      setDeleting(false);
    }
  };

  return (
    <Modal isOpen={dataset !== null} onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          <span>Delete dataset</span>
          <span className="text-xs font-normal text-default-400">
            {dataset?.title || "Untitled dataset"}
          </span>
        </ModalHeader>
        {dataset?.onExplore ? (
          <ModalBody>
            <p className="text-sm">
              This dataset is on Explore or awaiting review. Withdraw it from
              Explore first, then delete it.
            </p>
          </ModalBody>
        ) : (
          <ModalBody className="flex flex-col gap-4">
            <p className="text-sm">
              This permanently deletes the dataset and its files. Links to it
              will stop working. This can&apos;t be undone.
            </p>
            {overlayId && (
              <Checkbox
                isSelected={withOverlay}
                size="sm"
                onValueChange={setWithOverlay}
              >
                Also delete its linked single molecule dataset
              </Checkbox>
            )}
            <Input
              label={`Type "${CONFIRM_WORD}" to confirm`}
              labelPlacement="outside"
              placeholder={CONFIRM_WORD}
              size="sm"
              value={typed}
              onValueChange={setTyped}
            />
          </ModalBody>
        )}
        <ModalFooter>
          <Button variant="light" onPress={onClose}>
            Cancel
          </Button>
          {!dataset?.onExplore && (
            <Button
              color="danger"
              isDisabled={typed.trim().toLowerCase() !== CONFIRM_WORD}
              isLoading={deleting}
              onPress={confirm}
            >
              Delete
            </Button>
          )}
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
