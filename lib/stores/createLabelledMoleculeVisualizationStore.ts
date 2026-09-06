import { createStore } from "zustand";

import {
  findLowestAvailableSlot,
  getColorForSlot,
} from "@/lib/utils/gene-color-palette";

export type LmViewMode = "2D" | "3D";

/** How molecules that fail the filter are drawn. */
export type UnselectedMode = "grey" | "hidden";

/** The three menus, each backed by an obs column of the same name. */
export type LmMenu = "gene" | "domain" | "cell";

/** Gene lists are long enough to need ordering; the other two are short. */
export type LmSortBy = "alpha" | "count";

export const LM_MENUS: LmMenu[] = ["gene", "domain", "cell"];

export interface LabelledMoleculeVisualizationState {
  /** Which menu drives point colour. Every menu filters regardless. */
  colorBy: LmMenu;
  /** Which menu's panel is open, or null when closed. */
  openMenu: LmMenu | null;
  /** Ordering of the gene list. Counts first, as in the molecule viewer. */
  geneSortBy: LmSortBy;
  geneSortDir: "asc" | "desc";

  /**
   * Checked values per menu. Empty means "no constraint from this menu" — a
   * molecule has exactly one value per column, so checking several can only
   * mean OR; the menus then AND together.
   */
  selections: Record<LmMenu, Set<string>>;

  /** Per-value colour overrides, keyed by menu then value. */
  colorOverrides: Record<LmMenu, Record<string, string>>;
  /** Per-value size multipliers, keyed by menu then value. */
  sizeOverrides: Record<LmMenu, Record<string, number>>;
  /**
   * Palette slot per selected gene. Genes are coloured on selection rather
   * than from a fixed 213-entry palette, so the few that are checked stay
   * visually distinct.
   */
  geneColorSlots: Map<string, number>;

  unselectedMode: UnselectedMode;
  unselectedAlpha: number;
  globalScale: number;
  globalAlpha: number;
  viewMode: LmViewMode;
  searchTerm: Record<LmMenu, string>;
  /**
   * Incremented to ask the scene to re-frame the cloud. A nonce rather than a
   * boolean so repeated resets each fire, and so the scene needs no imperative
   * handle back into this store.
   */
  resetViewNonce: number;

  setColorBy: (menu: LmMenu) => void;
  setOpenMenu: (menu: LmMenu | null) => void;
  setGeneSort: (by: LmSortBy, dir: "asc" | "desc") => void;
  toggleValue: (menu: LmMenu, value: string) => void;
  setSelection: (menu: LmMenu, values: Set<string>) => void;
  clearMenu: (menu: LmMenu) => void;
  clearAll: () => void;
  setColorOverride: (menu: LmMenu, value: string, color: string) => void;
  setSizeOverride: (menu: LmMenu, value: string, size: number) => void;
  setUnselectedMode: (mode: UnselectedMode) => void;
  setUnselectedAlpha: (alpha: number) => void;
  setGlobalScale: (scale: number) => void;
  setGlobalAlpha: (alpha: number) => void;
  setViewMode: (mode: LmViewMode) => void;
  setSearchTerm: (menu: LmMenu, term: string) => void;
  resetView: () => void;
  /** Bulk apply, for restoring a shared URL in one commit. */
  applyUrlState: (patch: Partial<LabelledMoleculeVisualizationState>) => void;
  reset: () => void;
}

const emptySelections = (): Record<LmMenu, Set<string>> => ({
  gene: new Set(),
  domain: new Set(),
  cell: new Set(),
});

const emptyRecords = <T>(): Record<LmMenu, Record<string, T>> => ({
  gene: {},
  domain: {},
  cell: {},
});

const initialState = () => ({
  colorBy: "cell" as LmMenu,
  openMenu: null as LmMenu | null,
  geneSortBy: "count" as LmSortBy,
  geneSortDir: "desc" as "asc" | "desc",
  selections: emptySelections(),
  colorOverrides: emptyRecords<string>(),
  sizeOverrides: emptyRecords<number>(),
  geneColorSlots: new Map<string, number>(),
  unselectedMode: "grey" as UnselectedMode,
  unselectedAlpha: 0.2,
  globalScale: 1.0,
  globalAlpha: 1.0,
  viewMode: "3D" as LmViewMode,
  searchTerm: { gene: "", domain: "", cell: "" } as Record<LmMenu, string>,
  resetViewNonce: 0,
});

export function createLabelledMoleculeVisualizationStoreInstance() {
  return createStore<LabelledMoleculeVisualizationState>((set) => ({
    ...initialState(),

    setColorBy: (menu) => set({ colorBy: menu }),
    setOpenMenu: (menu) => set({ openMenu: menu }),

    setGeneSort: (by, dir) => set({ geneSortBy: by, geneSortDir: dir }),

    toggleValue: (menu, value) =>
      set((s) => {
        const next = new Set(s.selections[menu]);

        next.has(value) ? next.delete(value) : next.add(value);

        const patch: Partial<LabelledMoleculeVisualizationState> = {
          selections: { ...s.selections, [menu]: next },
        };

        // Genes take their colour from a slot assigned at selection time.
        if (menu === "gene") {
          const slots = new Map(s.geneColorSlots);

          if (next.has(value)) {
            if (!slots.has(value)) {
              slots.set(
                value,
                findLowestAvailableSlot(new Set(slots.values())),
              );
            }
          } else {
            slots.delete(value);
          }
          patch.geneColorSlots = slots;
        }

        return patch;
      }),

    setSelection: (menu, values) =>
      set((s) => {
        const patch: Partial<LabelledMoleculeVisualizationState> = {
          selections: { ...s.selections, [menu]: new Set(values) },
        };

        if (menu === "gene") {
          const slots = new Map<string, number>();

          for (const v of values) {
            slots.set(v, findLowestAvailableSlot(new Set(slots.values())));
          }
          patch.geneColorSlots = slots;
        }

        return patch;
      }),

    clearMenu: (menu) =>
      set((s) => ({
        selections: { ...s.selections, [menu]: new Set() },
        ...(menu === "gene" ? { geneColorSlots: new Map() } : {}),
      })),

    clearAll: () =>
      set({ selections: emptySelections(), geneColorSlots: new Map() }),

    setColorOverride: (menu, value, color) =>
      set((s) => ({
        colorOverrides: {
          ...s.colorOverrides,
          [menu]: { ...s.colorOverrides[menu], [value]: color },
        },
      })),

    setSizeOverride: (menu, value, size) =>
      set((s) => ({
        sizeOverrides: {
          ...s.sizeOverrides,
          [menu]: { ...s.sizeOverrides[menu], [value]: Math.max(0, size) },
        },
      })),

    setUnselectedMode: (mode) => set({ unselectedMode: mode }),
    setUnselectedAlpha: (alpha) =>
      set({ unselectedAlpha: Math.min(1, Math.max(0, alpha)) }),
    setGlobalScale: (scale) => set({ globalScale: Math.max(0, scale) }),
    setGlobalAlpha: (alpha) =>
      set({ globalAlpha: Math.min(1, Math.max(0, alpha)) }),
    setViewMode: (mode) => set({ viewMode: mode }),
    setSearchTerm: (menu, term) =>
      set((s) => ({ searchTerm: { ...s.searchTerm, [menu]: term } })),

    resetView: () => set((s) => ({ resetViewNonce: s.resetViewNonce + 1 })),

    applyUrlState: (patch) => set(patch),

    reset: () => set(initialState()),
  }));
}

/**
 * Colour for one value of a menu: an explicit override wins, then the
 * selection-time slot for genes, then the dataset's own palette.
 */
export function resolveValueColor(
  menu: LmMenu,
  value: string,
  overrides: Record<string, string>,
  geneColorSlots: Map<string, number>,
  datasetPalette: Record<string, string> | null,
): string {
  const override = overrides[value];

  if (override) return override;

  if (menu === "gene") {
    const slot = geneColorSlots.get(value);

    if (slot !== undefined) return getColorForSlot(slot);
  }

  return datasetPalette?.[value] ?? "#808080";
}
