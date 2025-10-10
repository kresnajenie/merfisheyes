"use client";

import { Navbar } from "@/components/navbar";
import { UploadSettingsModal } from "@/components/upload-settings-modal";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useState } from "react";

export function NavbarWrapper() {
  const [isUploadModalOpen, setIsUploadModalOpen] = useState(false);
  const currentDatasetId = useDatasetStore((state) => state.currentDatasetId);
  const getCurrentDataset = useDatasetStore((state) => state.getCurrentDataset);

  const dataset = getCurrentDataset();
  const hasDataset = currentDatasetId !== null && dataset !== null;

  return (
    <>
      <Navbar
        onUploadClick={hasDataset ? () => setIsUploadModalOpen(true) : undefined}
      />
      <UploadSettingsModal
        isOpen={isUploadModalOpen}
        onClose={() => setIsUploadModalOpen(false)}
        dataset={dataset}
      />
    </>
  );
}
