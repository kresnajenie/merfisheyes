"use client";

import { Button } from "@heroui/button";
import { Tooltip } from "@heroui/tooltip";
import { useRef } from "react";

import { glassButton } from "@/components/primitives";
import { SharePanel, useShareLink } from "@/components/share-panel";
import { usePanelId } from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

const BTN = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

/**
 * The buttons an image shown on its own has, where the cell and molecule
 * viewers have theirs: Split screen on the left, Hide UI and Share at the top
 * right. Its channels are set from the image layer's own panel. No camera
 * panel: the image is flat and is neither rotated nor flipped.
 *
 * `share`: false where the page's link would not open it elsewhere (an image
 * read from disk).
 */
export function ImageViewerControls({ share = true }: { share?: boolean }) {
  const { isSplitMode, enableSplit } = useSplitScreenStore();
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);
  const panelId = usePanelId();
  const topRightRef = useRef<HTMLDivElement>(null);
  const link = useShareLink();

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

      <div
        ref={topRightRef}
        data-ui-overlay
        className="absolute top-4 right-4 z-[var(--z-rail)] flex flex-row gap-2"
      >
        {/* Hide UI — strips overlays for screenshotting (H key) */}
        <Tooltip content="Hide UI for screenshot (H)" placement="bottom">
          <Button
            aria-label="Hide UI"
            className={`${BTN} ${glassButton()}`}
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

        {share && (
          <Tooltip
            content="Share this view (link + QR code)"
            placement="bottom"
          >
            <Button
              aria-label="Share"
              className={`${BTN} ${glassButton()}`}
              color="default"
              variant="light"
              onPress={link.toggle}
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
        )}

        {link.url && (
          <SharePanel
            controlsRef={topRightRef}
            url={link.url}
            onClose={link.close}
          />
        )}
      </div>
    </>
  );
}
