"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useState } from "react";
import {
  Button,
  Input,
  Modal,
  ModalBody,
  ModalContent,
  ModalFooter,
  ModalHeader,
  Progress,
} from "@heroui/react";
import NextLink from "next/link";
import { useSession } from "next-auth/react";

import { AuthModal } from "@/components/auth-modal";
import { saveLocalBundle, type SavedBundle } from "@/lib/local/save-bundle";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";

/**
 * "Upload & Save" for a folder of already-processed parts opened as a local
 * preview. Everything in it is uploaded as it is — there is nothing to
 * process — and the saved cells keep the molecules over them and the image
 * under them.
 */
export function BundleUploadModal({
  isOpen,
  onClose,
  cells,
}: {
  isOpen: boolean;
  onClose: () => void;
  /** The opened cells, when the preview has any. */
  cells: StandardizedDataset | null;
}) {
  const { data: session } = useSession();
  const bundle = useLocalBundleStore((s) => s.bundle);
  const [name, setName] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  const [progress, setProgress] = useState(0);
  const [message, setMessage] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [saved, setSaved] = useState<SavedBundle | null>(null);
  const [authOpen, setAuthOpen] = useState(false);

  if (!bundle) return null;

  const title = name ?? bundle.name.replace(/\.zarr$/, "");
  const parts = [
    bundle.cells &&
      (bundle.cells.kind === "spatialdata"
        ? "Cells and their image (SpatialData store)"
        : "Cells"),
    bundle.molecules !== undefined && "Molecules",
    bundle.image !== undefined && "Image",
  ].filter((p): p is string => !!p);

  const save = async () => {
    if (!session?.user) {
      setAuthOpen(true);

      return;
    }
    setBusy(true);
    setError(null);
    try {
      setSaved(
        await saveLocalBundle({
          bundle,
          name: title.trim() || "Untitled dataset",
          cells,
          onProgress: (fraction, text) => {
            setProgress(fraction * 100);
            setMessage(text);
          },
        }),
      );
    } catch (e) {
      setError(e instanceof Error ? e.message : "Upload failed");
    } finally {
      setBusy(false);
    }
  };

  return (
    <>
      <Modal
        isDismissable={!busy}
        isOpen={isOpen}
        size="md"
        onClose={busy ? undefined : onClose}
      >
        <ModalContent>
          <ModalHeader className="flex flex-col gap-1">
            <span>Upload &amp; Save</span>
            <span className="text-xs font-normal text-default-400">
              Uploaded as it is — nothing needs processing
            </span>
          </ModalHeader>
          <ModalBody className="gap-4">
            {saved ? (
              <div className="flex flex-col gap-3" data-testid="bundle-saved">
                <p className="text-sm text-success">
                  Saved.{" "}
                  {parts.length > 1 &&
                    "The parts are grouped in one project and open overlaid."}
                </p>
                <Button
                  as={NextLink}
                  color="primary"
                  href={saved.viewerUrl}
                  onPress={onClose}
                >
                  Open the saved dataset
                </Button>
              </div>
            ) : (
              <>
                <Input
                  isDisabled={busy}
                  label={parts.length > 1 ? "Project name" : "Dataset name"}
                  value={title}
                  onValueChange={setName}
                />
                <div className="rounded-xl border border-default-200 bg-default-100/50 p-3 text-sm">
                  <p className="text-default-500">This will upload:</p>
                  <ul className="mt-1 list-disc pl-5">
                    {parts.map((p) => (
                      <li key={p}>{p}</li>
                    ))}
                  </ul>
                </div>
                {busy && (
                  <div className="flex flex-col gap-1">
                    <Progress
                      aria-label="Upload progress"
                      color="primary"
                      value={progress}
                    />
                    <p className="text-xs text-default-500">{message}</p>
                  </div>
                )}
                {error && <p className="text-sm text-danger">{error}</p>}
              </>
            )}
          </ModalBody>
          {!saved && (
            <ModalFooter>
              <Button isDisabled={busy} variant="light" onPress={onClose}>
                Cancel
              </Button>
              <Button
                color="primary"
                isDisabled={!title.trim()}
                isLoading={busy}
                onPress={save}
              >
                {session?.user ? "Upload" : "Sign in to upload"}
              </Button>
            </ModalFooter>
          )}
        </ModalContent>
      </Modal>
      <AuthModal
        isOpen={authOpen}
        onAuthenticated={() => setAuthOpen(false)}
        onClose={() => setAuthOpen(false)}
      />
    </>
  );
}
