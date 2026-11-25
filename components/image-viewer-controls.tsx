"use client";

import { Button } from "@heroui/button";
import { Tooltip } from "@heroui/tooltip";
import { useRef, useState } from "react";

import { glassButton } from "@/components/primitives";
import { CameraPanel } from "@/components/camera-panel";
import { ViewerTopBar } from "@/components/viewer-top-bar";
import {
  usePanelId,
  usePanelVisualizationStore,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

const BTN = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

/**
 * The buttons an image shown on its own has, where the cell and molecule
 * viewers have theirs: Split screen on the left, and the same bar at the top
 * right (components/viewer-top-bar.tsx). Its channels are set from the image
 * layer's own panel. The camera panel has what applies to a flat image:
 * reset, rotation and flip (ImageScene reads them from the panel's
 * visualization store).
 *
 * `share`: false where the page's link would not open it elsewhere (an image
 * read from disk).
 */
export function ImageViewerControls({ share = true }: { share?: boolean }) {
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const panelId = usePanelId();
  const topRightRef = useRef<HTMLDivElement>(null);
  const [isCameraOpen, setIsCameraOpen] = useState(false);
  const {
    sceneRotation,
    setSceneRotation,
    flipX,
    setFlipX,
    flipY,
    setFlipY,
    resetCamera,
  } = usePanelVisualizationStore();

  return (
    <>
      {!isSplitMode && !panelId && (
        // Below the Image button, where the other viewers' rail starts
        <div
          data-ui-overlay
          className="absolute top-28 left-4 z-[var(--z-rail)] flex flex-col gap-2"
        >
          <Tooltip content="Split screen" placement="right">
            <Button
              aria-label="Split screen"
              className={`${BTN} ${glassButton()}`}
              color="default"
              variant="light"
              onPress={enableSplit}
            >
              <svg
                className="w-5 h-5"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M9 4.5v15m6-15v15M4.5 19.5h15a1.5 1.5 0 001.5-1.5V6a1.5 1.5 0 00-1.5-1.5h-15A1.5 1.5 0 003 6v12a1.5 1.5 0 001.5 1.5z"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
            </Button>
          </Tooltip>
        </div>
      )}

      <ViewerTopBar
        barRef={topRightRef}
        camera={{
          open: isCameraOpen,
          onToggle: () => setIsCameraOpen((open) => !open),
        }}
        shareUnavailable={
          share ? undefined : "An image opened from disk has no link to share"
        }
      >
        {isCameraOpen && (
          <CameraPanel
            controlsRef={topRightRef}
            flipX={flipX}
            flipY={flipY}
            placement="top-right"
            sceneRotation={sceneRotation}
            setFlipX={setFlipX}
            setFlipY={setFlipY}
            setSceneRotation={setSceneRotation}
            viewerKind="image"
            onClose={() => setIsCameraOpen(false)}
            onResetCamera={resetCamera}
          />
        )}
      </ViewerTopBar>
    </>
  );
}
