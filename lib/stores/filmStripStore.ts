import { create } from "zustand";

/**
 * Whether a film strip of sections takes the bottom edge of the page, so the
 * footer and the scene's bottom chrome sit above it.
 */
export const useFilmStripStore = create<{
  shown: boolean;
  setShown: (shown: boolean) => void;
}>((set) => ({
  shown: false,
  setShown: (shown) => set({ shown }),
}));
