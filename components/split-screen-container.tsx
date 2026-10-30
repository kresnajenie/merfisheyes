"use client";

import { useMemo } from "react";

import { ResizableDivider } from "./resizable-divider";
import { SplitPanelPicker } from "./split-panel-picker";
import { SplitPanelViewer } from "./split-panel-viewer";

import { PanelContext } from "@/lib/contexts/PanelContext";
import { createVisualizationStoreInstance } from "@/lib/stores/createVisualizationStore";
import { createDatasetStoreInstance } from "@/lib/stores/createDatasetStore";
import { createSingleMoleculeStoreInstance } from "@/lib/stores/createSingleMoleculeStore";
import { createSingleMoleculeVisualizationStoreInstance } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import { createLabelledMoleculeVisualizationStoreInstance } from "@/lib/stores/createLabelledMoleculeVisualizationStore";
import { createLayersStoreInstance } from "@/lib/stores/createLayersStore";
import { canSync, useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { useSyncCapabilities } from "@/lib/hooks/useSyncCapabilities";
import { useSyncVisualization } from "@/lib/hooks/useSyncVisualization";
import { useSyncSingleMoleculeVisualization } from "@/lib/hooks/useSyncSingleMoleculeVisualization";
import { useSyncLabelledMoleculeVisualization } from "@/lib/hooks/useSyncLabelledMoleculeVisualization";
import { useSyncImageChannels } from "@/lib/hooks/useSyncImageChannels";
import { useSplitUrlParams } from "@/lib/hooks/useSplitUrlParams";

/** A viewer page that can split: the split's state is kept in the URL. */
export function SplitScreenPage({ children }: { children: React.ReactNode }) {
  useSplitUrlParams();

  return <SplitScreenContainer>{children}</SplitScreenContainer>;
}

interface SplitScreenContainerProps {
  children: React.ReactNode;
}

export function SplitScreenContainer({ children }: SplitScreenContainerProps) {
  const {
    isSplitMode,
    dividerPosition,
    rightPanelDatasetId,
    rightPanelS3Url,
    rightPanelType,
    closeSplit,
  } = useSplitScreenStore();

  // Create stable store instances for the right panel
  const rightPanelStores = useMemo(
    () => ({
      panelId: "right",
      visualizationStore: createVisualizationStoreInstance(),
      datasetStore: createDatasetStoreInstance(),
      singleMoleculeStore: createSingleMoleculeStoreInstance(),
      singleMoleculeVisualizationStore:
        createSingleMoleculeVisualizationStoreInstance(),
      labelledMoleculeVisualizationStore:
        createLabelledMoleculeVisualizationStoreInstance(),
      layersStore: createLayersStoreInstance(),
    }),
    [],
  );

  // What the two panels can sync, from what each has open; each hook syncs
  // its own kind and the button shows when any can, so the two agree
  const syncable = useSyncCapabilities(rightPanelStores);

  // Bidirectional sync between left and right visualization stores
  useSyncVisualization(
    rightPanelStores.visualizationStore,
    rightPanelStores.datasetStore,
    syncable.cells,
  );
  useSyncSingleMoleculeVisualization(
    rightPanelStores.singleMoleculeVisualizationStore,
    rightPanelStores.singleMoleculeStore,
    syncable.molecules,
  );
  useSyncLabelledMoleculeVisualization(
    rightPanelStores.labelledMoleculeVisualizationStore,
    syncable.labelled,
  );
  useSyncImageChannels(rightPanelStores.layersStore, syncable.image);

  // Always render the same DOM tree so {children} (ThreeScene) never unmounts.
  // When not in split mode, the left panel takes 100% width and the right panel
  // is hidden with display:none (preserving the left panel's WebGL context).
  return (
    <div className="fixed inset-0 flex overflow-hidden">
      {/* Left panel — uses global stores (no PanelProvider) */}
      <div
        className="relative overflow-hidden h-full"
        style={{ width: isSplitMode ? `${dividerPosition}%` : "100%" }}
      >
        {children}
      </div>

      {isSplitMode && (
        <>
          <ResizableDivider canSync={canSync(syncable)} />

          {/* Right panel — uses scoped stores via PanelProvider */}
          <PanelContext.Provider value={rightPanelStores}>
            <div
              className="relative overflow-hidden h-full"
              style={{ width: `${100 - dividerPosition}%` }}
            >
              {/* Close button */}
              <button
                aria-label="Close split view"
                className="absolute top-4 right-4 z-[60] p-2 rounded-full bg-black/50 backdrop-blur-sm border border-white/20 hover:bg-white/10 transition-colors"
                data-ui-overlay
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
        </>
      )}
    </div>
  );
}
