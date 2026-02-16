import { create } from "zustand";

export type PanelType = "cell" | "sm";

interface SplitScreenState {
  isSplitMode: boolean;
  rightPanelDatasetId: string | null;
  rightPanelS3Url: string | null; // For from-s3 style loading
  rightPanelType: PanelType | null;
  dividerPosition: number; // percentage, default 50

  enableSplit: () => void;
  closeSplit: () => void;
  setRightPanel: (id: string | null, type: PanelType | null) => void;
  setRightPanelS3: (url: string, type: PanelType) => void;
  setDividerPosition: (pos: number) => void;
}

export const useSplitScreenStore = create<SplitScreenState>((set) => ({
  isSplitMode: false,
  rightPanelDatasetId: null,
  rightPanelS3Url: null,
  rightPanelType: null,
  dividerPosition: 50,

  enableSplit: () => set({ isSplitMode: true }),

  closeSplit: () =>
    set({
      isSplitMode: false,
      rightPanelDatasetId: null,
      rightPanelS3Url: null,
      rightPanelType: null,
      dividerPosition: 50,
    }),

  setRightPanel: (id, type) =>
    set({ rightPanelDatasetId: id, rightPanelS3Url: null, rightPanelType: type }),

  setRightPanelS3: (url, type) =>
    set({ rightPanelDatasetId: null, rightPanelS3Url: url, rightPanelType: type }),

  setDividerPosition: (pos) =>
    set({ dividerPosition: Math.max(25, Math.min(75, pos)) }),
}));
