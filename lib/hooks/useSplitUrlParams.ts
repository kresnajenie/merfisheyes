"use client";

import type { PanelType } from "@/lib/stores/splitScreenStore";

import { useEffect } from "react";
import { useRouter, useSearchParams } from "next/navigation";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";

/**
 * The split screen's share of the viewer URL (`split`, `splitS3Url`,
 * `splitType`, `sync`): read once on mount into the split store, and written
 * back whenever the split state changes. Every full-page viewer calls this,
 * so a split opened on one surface restores the same way on any other.
 */
export function useSplitUrlParams(): void {
  const router = useRouter();
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

  // Write split params to URL when split state changes
  useEffect(() => {
    // Use window.location.search as base to avoid stale Next.js searchParams
    const newParams = new URLSearchParams(window.location.search);

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
      router.replace(`?${newParams.toString()}`, { scroll: false });
    } else if (!isSplitMode) {
      newParams.delete("split");
      newParams.delete("splitS3Url");
      newParams.delete("splitType");
      newParams.delete("sync");
      const paramStr = newParams.toString();

      router.replace(paramStr ? `?${paramStr}` : window.location.pathname, {
        scroll: false,
      });
    }
  }, [
    isSplitMode,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    syncEnabled,
  ]);
}
