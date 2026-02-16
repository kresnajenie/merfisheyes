"use client";

import { Button } from "@heroui/react";

import { usePanelSingleMoleculeVisualizationStore } from "@/lib/hooks/usePanelStores";

export function ViewModeToggle() {
  const { viewMode, setViewMode } = usePanelSingleMoleculeVisualizationStore();

  return (
    <div className="absolute top-20 right-4 z-50">
      <Button
        className="font-semibold"
        color={viewMode === "2D" ? "primary" : "default"}
        size="sm"
        variant={viewMode === "2D" ? "solid" : "bordered"}
        onPress={() => setViewMode(viewMode === "2D" ? "3D" : "2D")}
      >
        {viewMode === "2D" ? "2D View" : "3D View"}
        <span className="ml-2 text-xs opacity-70">
          (Switch to {viewMode === "2D" ? "3D" : "2D"})
        </span>
      </Button>
    </div>
  );
}
