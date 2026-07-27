"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useState } from "react";
import { usePathname } from "next/navigation";

import { Navbar } from "@/components/navbar";
import { ViewerNavbar } from "@/components/viewer-navbar";
import { UploadSettingsModal } from "@/components/upload-settings-modal";
import { SingleMoleculeUploadModal } from "@/components/single-molecule-upload-modal";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";

export function NavbarWrapper() {
  const [isUploadModalOpen, setIsUploadModalOpen] = useState(false);
  const pathname = usePathname();

  // Check if we're in single molecule viewer based on path
  const isSingleMolecule = pathname?.startsWith("/sm-viewer");

  // Both viewer families get the compact top-left brand pill instead of the
  // centered nav bar.
  const isViewer =
    pathname?.startsWith("/viewer") || pathname?.startsWith("/sm-viewer");

  // Get dataset from appropriate store
  const cellDatasetId = useDatasetStore((state) => state.currentDatasetId);
  const getCellDataset = useDatasetStore((state) => state.getCurrentDataset);

  const smDatasetId = useSingleMoleculeStore((state) => state.currentDatasetId);
  const getSmDataset = useSingleMoleculeStore(
    (state) => state.getCurrentDataset,
  );

  const currentDatasetId = isSingleMolecule ? smDatasetId : cellDatasetId;
  const hasDataset = currentDatasetId !== null;

  // Helper to ensure we only pass StandardizedDataset to UploadSettingsModal
  const getCellDatasetAsStandardized = (): StandardizedDataset | null => {
    const dataset = getCellDataset();

    // Type guard: check if it's a StandardizedDataset
    if (dataset && "spatial" in dataset && "embeddings" in dataset) {
      return dataset as StandardizedDataset;
    }

    return null;
  };

  const isPlayback = useVisualizationStore((s) => s.celltypePlayback);

  if (isPlayback) return null;

  return (
    <>
      <div data-ui-overlay>
        {isViewer ? (
          <ViewerNavbar onUploadClick={() => setIsUploadModalOpen(true)} />
        ) : (
          <Navbar onUploadClick={() => setIsUploadModalOpen(true)} />
        )}
      </div>
      {isSingleMolecule ? (
        <SingleMoleculeUploadModal
          dataset={getSmDataset()}
          isOpen={isUploadModalOpen}
          onClose={() => setIsUploadModalOpen(false)}
        />
      ) : (
        <UploadSettingsModal
          dataset={getCellDatasetAsStandardized()}
          isOpen={isUploadModalOpen}
          onClose={() => setIsUploadModalOpen(false)}
        />
      )}
    </>
  );
}
