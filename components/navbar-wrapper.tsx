"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";

import { useState } from "react";
import { usePathname } from "next/navigation";

import { Navbar } from "@/components/navbar";
import { ViewerNavbar } from "@/components/viewer-navbar";
import { UploadSettingsModal } from "@/components/upload-settings-modal";
import { SingleMoleculeUploadModal } from "@/components/single-molecule-upload-modal";
import { BundleUploadModal } from "@/components/bundle-upload-modal";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useVisualizationStore } from "@/lib/stores/visualizationStore";

export function NavbarWrapper() {
  const [isUploadModalOpen, setIsUploadModalOpen] = useState(false);
  const pathname = usePathname();

  // Check if we're in single molecule viewer based on path
  const isSingleMolecule = pathname?.startsWith("/sm-viewer");

  // Every viewer family gets the compact top-left brand pill instead of the
  // centered nav bar.
  const isViewer =
    pathname?.startsWith("/viewer") ||
    pathname?.startsWith("/sm-viewer") ||
    pathname?.startsWith("/lm-viewer") ||
    pathname?.startsWith("/spatialdata-viewer") ||
    pathname?.startsWith("/image-viewer");

  // Get dataset from appropriate store
  const cellDatasetId = useDatasetStore((state) => state.currentDatasetId);
  const getCellDataset = useDatasetStore((state) => state.getCurrentDataset);

  const smDatasetId = useSingleMoleculeStore((state) => state.currentDatasetId);
  const getSmDataset = useSingleMoleculeStore(
    (state) => state.getCurrentDataset,
  );

  const currentDatasetId = isSingleMolecule ? smDatasetId : cellDatasetId;

  // A folder of already-processed parts opened as a local preview is saved
  // as it is, all parts together — not through the per-dataset dialogs, which
  // process what they are given.
  const bundleAnchorId = useLocalBundleStore((s) =>
    s.bundle ? s.anchorId : undefined,
  );
  const isBundle =
    pathname === "/image-viewer/from-local"
      ? bundleAnchorId !== undefined
      : !!pathname?.endsWith("/from-local") &&
        bundleAnchorId != null &&
        bundleAnchorId === currentDatasetId;
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
      {isBundle ? (
        <BundleUploadModal
          cells={isSingleMolecule ? null : getCellDatasetAsStandardized()}
          isOpen={isUploadModalOpen}
          onClose={() => setIsUploadModalOpen(false)}
        />
      ) : isSingleMolecule ? (
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
