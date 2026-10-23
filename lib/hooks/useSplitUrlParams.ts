"use client";

import type { PanelType } from "@/lib/stores/splitScreenStore";

import { useEffect } from "react";
import { useSearchParams } from "next/navigation";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

/**
 * The split screen's share of the viewer URL (`split`, `splitS3Url`,
 * `splitType`, `sync`): read once on mount into the split store, and written
 * back whenever the split state changes. Every full-page viewer calls this,
 * so a split opened on one surface restores the same way on any other.
 */
export function useSplitUrlParams(): void {
  const searchParams = useSearchParams();
  const {
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
    enableSplit,
    setRightPanel,
    setRightPanelS3,
    setSyncEnabled,
    setSyncFromUrl,
  } = useSplitScreenStore();

  // Read split params from URL on mount
  useEffect(() => {
    const splitId = searchParams.get("split");
    const splitS3Url = searchParams.get("splitS3Url");
    const splitType = searchParams.get("splitType") as PanelType | null;

    if (splitS3Url && splitType) {
      enableSplit();
      setRightPanelS3(decodeURIComponent(splitS3Url), splitType);
    } else if (splitId && splitType) {
      enableSplit();
      setRightPanel(splitId, splitType);
    }
    if (searchParams.get("sync") === "1") {
      setSyncEnabled(true);
      setSyncFromUrl(true);
    }
  }, []);

  // Write split params to URL when split state changes. Written at once
  // with history.replaceState, as lib/utils/url-state-writer.ts writes the
  // panels' view state: both read the address as it is when they write, so
  // neither undoes the other. (router.replace applied later, and a panel's
  // image settings written in between were lost — or the split was.)
  useEffect(() => {
    const url = new URL(window.location.href);
    const newParams = url.searchParams;

    if (isSplitMode && rightPanelType) {
      if (rightPanelS3Url) {
        newParams.set("splitS3Url", encodeURIComponent(rightPanelS3Url));
        newParams.delete("split");
      } else if (rightPanelDatasetId) {
        newParams.set("split", rightPanelDatasetId);
        newParams.delete("splitS3Url");
      }
      newParams.set("splitType", rightPanelType);
      if (syncEnabled) {
        newParams.set("sync", "1");
      } else {
        newParams.delete("sync");
      }
    } else if (!isSplitMode) {
      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      newParams.delete("sync");
    } else {
      return;
    }
    if (url.href !== window.location.href) {
      window.history.replaceState(window.history.state, "", url.href);
    }
  }, [
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
  ]);
}
