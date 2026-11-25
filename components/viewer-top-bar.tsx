"use client";

import type { ReactNode, RefObject } from "react";

import { Button } from "@heroui/button";
import { Tooltip } from "@heroui/tooltip";

import { DatasetInfoButton } from "@/components/dataset-info-panel";
import { glassButton } from "@/components/primitives";
import { SharePanel, useShareLink } from "@/components/share-panel";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

const BUTTON = "w-14 h-14 min-w-0 rounded-full font-medium text-xs";

/** One round button of the bar; `active` while the flyout it opens is up. */
export function TopBarButton({
  label,
  tooltip,
  active = false,
  onPress,
  children,
}: {
  label: string;
  tooltip: string;
  active?: boolean;
  onPress: () => void;
  children: ReactNode;
}) {
  return (
    <Tooltip content={tooltip} placement="bottom">
      <Button
        aria-label={label}
        className={`${BUTTON} ${active ? "" : glassButton()}`}
        color={active ? "primary" : "default"}
        variant={active ? "shadow" : "light"}
        onPress={onPress}
      >
        {children}
      </Button>
    </Tooltip>
  );
}

const Icon = ({ children }: { children: ReactNode }) => (
  <svg
    className="w-5 h-5"
    fill="none"
    stroke="currentColor"
    strokeWidth={1.5}
    viewBox="0 0 24 24"
  >
    {children}
  </svg>
);

/**
 * The buttons at the top right of every viewer (cells, molecules, an image on
 * its own), always these and in this order: Camera, Hide UI, About this
 * dataset, Share. A viewer adds what only it has through `before` and
 * `afterHide`, and puts the flyouts its buttons open in `children`, which are
 * anchored to the bar.
 *
 * `shareUnavailable`: why this view has no link (a file read from disk). The
 * button stays where it is, disabled, and says so.
 */
export function ViewerTopBar({
  barRef,
  camera,
  before,
  afterHide,
  shareUnavailable,
  children,
}: {
  barRef: RefObject<HTMLDivElement>;
  camera: { open: boolean; onToggle: () => void };
  before?: ReactNode;
  afterHide?: ReactNode;
  shareUnavailable?: string;
  children?: ReactNode;
}) {
  const setHideUi = useSplitScreenStore((s) => s.setHideUi);
  const share = useShareLink();

  return (
    <div
      ref={barRef}
      data-ui-overlay
      className="absolute top-4 right-4 z-[var(--z-rail)] flex flex-row gap-2"
      data-testid="viewer-top-bar"
    >
      {before}

      <TopBarButton
        active={camera.open}
        label="Camera controls"
        tooltip="Camera controls"
        onPress={camera.onToggle}
      >
        <Icon>
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
        </Icon>
      </TopBarButton>

      {/* Strips the overlays for a screenshot (H key) */}
      <TopBarButton
        label="Hide UI"
        tooltip="Hide UI for screenshot (H)"
        onPress={() => setHideUi(true)}
      >
        <Icon>
          <path
            d="M3.98 8.223A10.477 10.477 0 001.934 12C3.226 16.338 7.244 19.5 12 19.5c.993 0 1.953-.138 2.863-.395M6.228 6.228A10.45 10.45 0 0112 4.5c4.756 0 8.773 3.162 10.065 7.498a10.523 10.523 0 01-4.293 5.774M6.228 6.228L3 3m3.228 3.228l3.65 3.65m7.894 7.894L21 21m-3.228-3.228l-3.65-3.65m0 0a3 3 0 10-4.243-4.243m4.242 4.242L9.88 9.88"
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </Icon>
      </TopBarButton>

      {afterHide}

      {/* What the store says about itself; absent when it says nothing */}
      <DatasetInfoButton />

      {/* The link to this view, with a QR code */}
      <Tooltip
        content={shareUnavailable ?? "Share this view (link + QR code)"}
        placement="bottom"
      >
        {/* The span keeps the tooltip on a disabled button */}
        <span>
          <Button
            aria-label="Share"
            className={`${BUTTON} ${glassButton()}`}
            color="default"
            isDisabled={!!shareUnavailable}
            variant="light"
            onPress={share.toggle}
          >
            <Icon>
              <path
                d="M7.217 10.907a2.25 2.25 0 100 2.186m0-2.186c.18.324.283.696.283 1.093s-.103.77-.283 1.093m0-2.186l9.566-5.314m-9.566 7.5l9.566 5.314m0 0a2.25 2.25 0 103.935 2.186 2.25 2.25 0 00-3.935-2.186zm0-12.814a2.25 2.25 0 103.933-2.185 2.25 2.25 0 00-3.933 2.185z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </Icon>
          </Button>
        </span>
      </Tooltip>

      {children}

      {share.url && (
        <SharePanel
          controlsRef={barRef}
          url={share.url}
          onClose={share.close}
        />
      )}
    </div>
  );
}
