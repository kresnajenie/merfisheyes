import type { TourConfig } from "@/lib/tour/types";

import { create } from "zustand";

export type TourPhase = "idle" | "fly" | "wait" | "orbit" | "zoom-out";

interface TourState {
  config: TourConfig | null;
  configName: string | null;
  isPlaying: boolean;
  currentStopIndex: number;
  currentPhase: TourPhase;

  setConfig: (config: TourConfig | null, configName?: string | null) => void;
  start: () => void;
  stop: () => void;
  setProgress: (stopIndex: number, phase: TourPhase) => void;
  finish: () => void;
}

export const useTourStore = create<TourState>((set, get) => ({
  config: null,
  configName: null,
  isPlaying: false,
  currentStopIndex: 0,
  currentPhase: "idle",

  setConfig: (config, configName = null) =>
    set({
      config,
      configName,
      isPlaying: false,
      currentStopIndex: 0,
      currentPhase: "idle",
    }),

  start: () => {
    const { config, isPlaying } = get();

    if (!config || isPlaying) return;
    set({ isPlaying: true, currentStopIndex: 0, currentPhase: "fly" });
  },

  stop: () => set({ isPlaying: false, currentPhase: "idle" }),

  setProgress: (stopIndex, phase) =>
    set({ currentStopIndex: stopIndex, currentPhase: phase }),

  finish: () => set({ isPlaying: false, currentPhase: "idle" }),
}));
