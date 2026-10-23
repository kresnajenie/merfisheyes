import { create } from "zustand";

/**
 * What the right panel shows, as its full-page viewer would: cells
 * (`/viewer`), molecules (`/sm-viewer`), labelled molecules (`/lm-viewer`), a
 * SpatialData store (`/spatialdata-viewer`) or an image (`/image-viewer`).
 */
export type PanelType = "cell" | "sm" | "lm" | "spatialdata" | "image";

/** The panel for a dataset of this `datasetType` (lib/datasets/kind.ts). */
export function panelTypeForDataset(
  datasetType: string | null | undefined,
): PanelType {
  if (datasetType === "single_molecule") return "sm";
  if (datasetType === "spatialdata") return "spatialdata";
  if (datasetType === "image") return "image";

  return "cell";
}

/** The panel a viewer link opens, by its first path segment; null if none. */
export function panelTypeForViewer(segment: string): PanelType | null {
  const panels: Record<string, PanelType> = {
    viewer: "cell",
    "sm-viewer": "sm",
    "lm-viewer": "lm",
    "spatialdata-viewer": "spatialdata",
    "image-viewer": "image",
  };

  return panels[segment] ?? null;
}

interface SplitScreenState {
  isSplitMode: boolean;
  rightPanelDatasetId: string | null;
  rightPanelS3Url: string | null; // For from-s3 style loading
  rightPanelType: PanelType | null;
  dividerPosition: number; // percentage, default 50
  syncEnabled: boolean;
  syncFromUrl: boolean; // true when sync was restored from URL (needs settling)

  enableSplit: () => void;
  closeSplit: () => void;
  setRightPanel: (id: string | null, type: PanelType | null) => void;
  setRightPanelS3: (url: string, type: PanelType) => void;
  setDividerPosition: (pos: number) => void;
  toggleSync: () => void;
  setSyncEnabled: (enabled: boolean) => void;
  setSyncFromUrl: (fromUrl: boolean) => void;

  // Global "hide UI" mode — strips all overlay UI for clean screenshots.
  hideUi: boolean;
  setHideUi: (hidden: boolean) => void;
  toggleHideUi: () => void;
}

export const useSplitScreenStore = create<SplitScreenState>((set) => ({
  isSplitMode: false,
  rightPanelDatasetId: null,
  rightPanelS3Url: null,
  rightPanelType: null,
  dividerPosition: 50,
  syncEnabled: false,
  syncFromUrl: false,

  enableSplit: () => set({ isSplitMode: true }),

  closeSplit: () =>
    set({
      isSplitMode: false,
      rightPanelDatasetId: null,
      rightPanelS3Url: null,
      rightPanelType: null,
      dividerPosition: 50,
      syncEnabled: false,
      syncFromUrl: false,
    }),

  setRightPanel: (id, type) =>
    set({
      rightPanelDatasetId: id,
      rightPanelS3Url: null,
      rightPanelType: type,
    }),

  setRightPanelS3: (url, type) =>
    set({
      rightPanelDatasetId: null,
      rightPanelS3Url: url,
      rightPanelType: type,
    }),

  setDividerPosition: (pos) =>
    set({ dividerPosition: Math.max(25, Math.min(75, pos)) }),

  toggleSync: () => set((state) => ({ syncEnabled: !state.syncEnabled })),
  setSyncEnabled: (enabled) => set({ syncEnabled: enabled }),
  setSyncFromUrl: (fromUrl) => set({ syncFromUrl: fromUrl }),

  hideUi: false,
  setHideUi: (hidden) => set({ hideUi: hidden }),
  toggleHideUi: () => set((state) => ({ hideUi: !state.hideUi })),
}));
