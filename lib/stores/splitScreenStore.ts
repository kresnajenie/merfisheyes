import { create } from "zustand";

export type PanelType = "cell" | "sm";

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
}));
