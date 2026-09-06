"use client";

import type { ViewerConfig } from "@/lib/utils/viewer-config";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { Tooltip } from "@heroui/tooltip";
import { useSession } from "next-auth/react";
import { useState } from "react";
import { toast } from "react-toastify";

import { glassButton, glassPanel } from "@/components/primitives";
import { useLabelledMoleculeVisualizationStore } from "@/lib/stores/labelledMoleculeVisualizationStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useViewerRegistrationStore } from "@/lib/stores/viewerRegistrationStore";

// Same rail geometry as visualization-controls.tsx / single-molecule-controls.tsx.
const buttonBaseClass = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

/**
 * Top-right cluster: camera, hide-UI, share — matching the single-cell and
 * single-molecule viewers.
 *
 * The camera flyout is trimmed relative to `CameraPanel`: this scene has no
 * rotation or flip support yet, so only the controls that actually do
 * something are offered.
 */
export default function LabelledMoleculeTopControls() {
  const {
    resetView,
    camera,
    selectedScale,
    setSelectedScale,
    unselectedScale,
    setUnselectedScale,
    showMeshes,
    setShowMeshes,
    meshMode,
    setMeshMode,
    meshOpacity,
    setMeshOpacity,
  } = useLabelledMoleculeVisualizationStore();
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);
  const [isCameraOpen, setIsCameraOpen] = useState(false);
  const [saving, setSaving] = useState(false);

  const { data: session } = useSession();
  const { dbId, ownerId, adminOwned, viewerConfig } =
    useViewerRegistrationStore();
  const userId = (session?.user as { id?: string } | undefined)?.id;
  const isAdmin = !!(session?.user as { isAdmin?: boolean } | undefined)
    ?.isAdmin;
  // Same gate the single-cell camera panel uses.
  const canSave =
    !!dbId && ((!!ownerId && ownerId === userId) || (adminOwned && isAdmin));

  const saveDefaults = async () => {
    if (!dbId) return;
    setSaving(true);
    try {
      const merged: ViewerConfig = {
        ...(viewerConfig ?? { version: 1 }),
        version: 1,
        ...(camera ? { camera } : {}),
      };
      const res = await fetch(`/api/ingest/mine/${dbId}`, {
        method: "PATCH",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ viewerConfig: merged }),
      });

      if (!res.ok) {
        toast.error("Couldn't save defaults.");

        return;
      }
      useViewerRegistrationStore.getState().set({ viewerConfig: merged });
      toast.success("Saved as this dataset's defaults.");
    } catch {
      toast.error("Couldn't save defaults.");
    } finally {
      setSaving(false);
    }
  };

  const handleShare = async () => {
    try {
      await navigator.clipboard.writeText(window.location.href);
      toast.success("Link copied");
    } catch {
      toast.error("Couldn't copy link");
    }
  };

  return (
    <div
      data-ui-overlay
      className="absolute top-4 right-4 z-[var(--z-rail)] flex flex-row gap-2"
    >
      <Tooltip content="Camera controls" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${isCameraOpen ? "" : glassButton()}`}
          color={isCameraOpen ? "primary" : "default"}
          variant={isCameraOpen ? "shadow" : "light"}
          onPress={() => setIsCameraOpen((v) => !v)}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M6.827 6.175A2.31 2.31 0 015.186 7.23c-.38.054-.757.112-1.134.175C2.999 7.58 2.25 8.507 2.25 9.574V18a2.25 2.25 0 002.25 2.25h15A2.25 2.25 0 0021.75 18V9.574c0-1.067-.75-1.994-1.802-2.169a47.865 47.865 0 00-1.134-.175 2.31 2.31 0 01-1.64-1.055l-.822-1.316a2.192 2.192 0 00-1.736-1.039 48.774 48.774 0 00-5.232 0 2.192 2.192 0 00-1.736 1.039l-.821 1.316z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
            <path
              d="M16.5 12.75a4.5 4.5 0 11-9 0 4.5 4.5 0 019 0z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      <Tooltip content="Hide UI for screenshot (H)" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${glassButton()}`}
          color="default"
          variant="light"
          onPress={() => setHideUi(true)}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M3.98 8.223A10.477 10.477 0 001.934 12C3.226 16.338 7.244 19.5 12 19.5c.993 0 1.953-.138 2.863-.395M6.228 6.228A10.45 10.45 0 0112 4.5c4.756 0 8.773 3.162 10.065 7.498a10.523 10.523 0 01-4.293 5.774M6.228 6.228L3 3m3.228 3.228l3.65 3.65m7.894 7.894L21 21m-3.228-3.228l-3.65-3.65m0 0a3 3 0 10-4.243-4.243m4.242 4.242L9.88 9.88"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      <Tooltip content="Copy link to this view" placement="bottom">
        <Button
          className={`${buttonBaseClass} ${glassButton()}`}
          color="default"
          variant="light"
          onPress={handleShare}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={1.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M7.217 10.907a2.25 2.25 0 100 2.186m0-2.186c.18.324.283.696.283 1.093s-.103.77-.283 1.093m0-2.186l9.566-5.314m-9.566 7.5l9.566 5.314m0 0a2.25 2.25 0 103.935 2.186 2.25 2.25 0 00-3.935-2.186zm0-12.814a2.25 2.25 0 103.933-2.185 2.25 2.25 0 00-3.933 2.185z"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </Button>
      </Tooltip>

      {isCameraOpen && (
        <div className={`absolute top-16 right-0 w-[240px] ${glassPanel()}`}>
          <div className="p-4 space-y-4">
            <Slider
              label="Selected size"
              maxValue={4}
              minValue={0.1}
              size="sm"
              step={0.05}
              value={selectedScale}
              onChange={(v) => setSelectedScale(Number(v))}
            />

            <Slider
              label="Unselected size"
              maxValue={4}
              minValue={0}
              size="sm"
              step={0.05}
              value={unselectedScale}
              onChange={(v) => setUnselectedScale(Number(v))}
            />

            <Button
              className="w-full"
              size="sm"
              variant="ghost"
              onPress={resetView}
            >
              Reset view
            </Button>

            {/* Cell segmentation surfaces. They follow the cell menu: every
                cell when nothing is selected, otherwise just the selection. */}
            <div className="space-y-2 border-t border-default-200/40 pt-3">
              <div className="flex items-center justify-between">
                <span className="text-xs text-default-500">Cell meshes</span>
                <Button
                  color={showMeshes ? "primary" : "default"}
                  size="sm"
                  variant={showMeshes ? "flat" : "light"}
                  onPress={() => setShowMeshes(!showMeshes)}
                >
                  {showMeshes ? "On" : "Off"}
                </Button>
              </div>

              {showMeshes && (
                <>
                  <div className="flex gap-1">
                    {(
                      [
                        ["wireframe", "Wireframe"],
                        ["translucent", "Translucent"],
                      ] as const
                    ).map(([mode, label]) => (
                      <Button
                        key={mode}
                        className="flex-1"
                        color={meshMode === mode ? "primary" : "default"}
                        size="sm"
                        variant={meshMode === mode ? "flat" : "light"}
                        onPress={() => setMeshMode(mode)}
                      >
                        {label}
                      </Button>
                    ))}
                  </div>

                  {meshMode === "translucent" && (
                    <Slider
                      label="Mesh opacity"
                      maxValue={1}
                      minValue={0.02}
                      size="sm"
                      step={0.02}
                      value={meshOpacity}
                      onChange={(v) => setMeshOpacity(Number(v))}
                    />
                  )}
                </>
              )}
            </div>

            {canSave && (
              <Button
                className="w-full"
                color="primary"
                isDisabled={saving}
                size="sm"
                variant="flat"
                onPress={saveDefaults}
              >
                {saving ? "Saving…" : "Save current view as default"}
              </Button>
            )}

            <p className="text-[11px] text-default-500">
              Rotate and flip aren&apos;t available for this dataset type yet.
            </p>
          </div>
        </div>
      )}
    </div>
  );
}
