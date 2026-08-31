"use client";

import { Button } from "@heroui/button";
import { useSession } from "next-auth/react";
import { useState } from "react";
import { createPortal } from "react-dom";

import { OverlayManager } from "@/components/overlay-manager";

import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";
import { glassButton } from "@/components/primitives";

/**
 * Floating "Overlay" control in the single-cell viewer. Shown only to the
 * dataset's owner; opens the shared OverlayManager so they can attach, change,
 * or remove a single-molecule overlay while viewing.
 */
export function ViewerOverlayControl({ scDatasetId }: { scDatasetId: string }) {
  const { data: session } = useSession();
  const ownerId = useViewerRegistrationStore((s) => s.ownerId);
  const registered = useViewerRegistrationStore((s) => s.registered);
  const [open, setOpen] = useState(false);

  const isOwner =
    registered && !!session?.user?.id && session.user.id === ownerId;

  if (!isOwner) return null;

  return (
    <>
      <div className="absolute bottom-4 right-4 z-[var(--z-rail)]" data-ui-overlay>
        <Button
          className={glassButton()}
          size="sm"
          variant="light"
          onPress={() => setOpen(true)}
        >
          Overlay
        </Button>
      </div>

      {open &&
        createPortal(
          <div
            className="fixed inset-0 z-[var(--z-modal-top)] flex items-center justify-center bg-black/70 p-6 backdrop-blur-md"
            role="presentation"
            onClick={(e) => {
              if (e.target === e.currentTarget) setOpen(false);
            }}
          >
            <div
              className="glass-panel w-full max-w-md rounded-2xl p-6"
              role="dialog"
            >
              <div className="flex items-center justify-between">
                <h2 className="text-lg font-semibold">Manage overlay</h2>
                <button
                  className="text-default-400 hover:text-foreground"
                  type="button"
                  onClick={() => setOpen(false)}
                >
                  ✕
                </button>
              </div>
              <div className="mt-4">
                <OverlayManager scDatasetId={scDatasetId} />
              </div>
            </div>
          </div>,
          document.body,
        )}
    </>
  );
}
