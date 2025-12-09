"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useState } from "react";
import { usePathname } from "next/navigation";

import { Navbar } from "@/components/navbar";
import { UploadSettingsModal } from "@/components/upload-settings-modal";
import { SingleMoleculeUploadModal } from "@/components/single-molecule-upload-modal";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";

export function NavbarWrapper() {
  const [isUploadModalOpen, setIsUploadModalOpen] = useState(false);
  const pathname = usePathname();

  // Check if we're in single molecule viewer based on path
  const isSingleMolecule = pathname?.startsWith("/sm-viewer");

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

  return (
    <>
      <Navbar
        onUploadClick={() => setIsUploadModalOpen(true)}
      />
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
