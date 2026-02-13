"use client";

import { useMemo } from "react";

import { PanelContext } from "@/lib/contexts/PanelContext";
import { createVisualizationStoreInstance } from "@/lib/stores/createVisualizationStore";
import { createDatasetStoreInstance } from "@/lib/stores/createDatasetStore";
import { createSingleMoleculeStoreInstance } from "@/lib/stores/createSingleMoleculeStore";
import { createSingleMoleculeVisualizationStoreInstance } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { ResizableDivider } from "./resizable-divider";
import { SplitPanelPicker } from "./split-panel-picker";
import { SplitPanelViewer } from "./split-panel-viewer";

interface SplitScreenContainerProps {
  children: React.ReactNode;
}

export function SplitScreenContainer({ children }: SplitScreenContainerProps) {
  const { isSplitMode, dividerPosition, rightPanelDatasetId, rightPanelS3Url, rightPanelType, closeSplit } =
    useSplitScreenStore();

  // Create stable store instances for the right panel
  const rightPanelStores = useMemo(
    () => ({
      panelId: "right",
      visualizationStore: createVisualizationStoreInstance(),
      datasetStore: createDatasetStoreInstance(),
      singleMoleculeStore: createSingleMoleculeStoreInstance(),
      singleMoleculeVisualizationStore:
        createSingleMoleculeVisualizationStoreInstance(),
    }),
    [],
  );

  if (!isSplitMode) {
    // Single panel mode — render children directly in a fixed container
    return (
      <div className="fixed inset-0 overflow-hidden">
        {children}
      </div>
    );
  }

  // Split mode — two panels with divider
  return (
    <div className="fixed inset-0 flex overflow-hidden">
      {/* Left panel — uses global stores (no PanelProvider) */}
      <div
        className="relative overflow-hidden h-full"
        style={{ width: `${dividerPosition}%` }}
      >
        {children}
      </div>

      <ResizableDivider />

      {/* Right panel — uses scoped stores via PanelProvider */}
      <PanelContext.Provider value={rightPanelStores}>
        <div
          className="relative overflow-hidden h-full"
          style={{ width: `${100 - dividerPosition}%` }}
        >
          {/* Close button */}
          <button
            className="absolute top-4 right-4 z-[60] p-2 rounded-full bg-black/50 backdrop-blur-sm border border-white/20 hover:bg-white/10 transition-colors"
            onClick={closeSplit}
          >
            <svg
              className="w-4 h-4 text-white"
              fill="none"
              stroke="currentColor"
              strokeWidth={2}
              viewBox="0 0 24 24"
            >
              <path
                d="M6 18L18 6M6 6l12 12"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </button>

          {(rightPanelDatasetId || rightPanelS3Url) && rightPanelType ? (
            <SplitPanelViewer
              datasetId={rightPanelDatasetId}
              s3Url={rightPanelS3Url}
              type={rightPanelType}
            />
          ) : (
            <SplitPanelPicker />
          )}
        </div>
      </PanelContext.Provider>
    </div>
  );
}
