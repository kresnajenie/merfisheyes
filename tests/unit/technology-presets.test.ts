import { describe, expect, it } from "vitest";

import {
  TECHNOLOGY_PRESETS,
  applyTechnologyPreset,
  sizeLocked,
  technologyOf,
  technologyPresetOf,
} from "@/lib/config/technology-presets";
import { COLORMAPS } from "@/lib/utils/colormaps";
import {
  decodeCellVizState,
  encodeCellVizState,
} from "@/lib/utils/url-visualization-state";
import { createVisualizationStoreInstance } from "@/lib/stores/createVisualizationStore";

const bins = { cellFootprint: { shape: "square" as const, size: 8 } };
const spots = { cellFootprint: { shape: "circle" as const, size: 55 } };


describe("technology presets", () => {
  it("names square bins Visium HD, round spots Visium, and nothing else", () => {
    expect(technologyOf(bins)).toBe("Visium HD");
    expect(technologyOf(spots)).toBe("Visium");
    expect(technologyOf({ cellFootprint: null })).toBeNull();
    expect(technologyPresetOf({ cellFootprint: null })).toEqual({});
  });

  it("locks spots and bins to their real size, by their preset", () => {
    expect(sizeLocked(spots)).toBe(true);
    expect(sizeLocked(bins)).toBe(true);
    // Dots have no real size to keep
    expect(sizeLocked({ cellFootprint: null })).toBe(false);
    expect(sizeLocked(null)).toBe(false);
  });

  it("only names colormaps the viewer has, and opacities in 0-1", () => {
    for (const preset of Object.values(TECHNOLOGY_PRESETS)) {
      if (preset.colormap) expect(Object.keys(COLORMAPS)).toContain(preset.colormap);
      if (preset.cellOpacity !== undefined) {
        expect(preset.cellOpacity).toBeGreaterThan(0);
        expect(preset.cellOpacity).toBeLessThanOrEqual(1);
      }
    }
  });

  it("opens Visium HD with its preset", () => {
    const store = createVisualizationStoreInstance();

    applyTechnologyPreset(technologyPresetOf(bins), store.getState());
    const s = store.getState();

    expect(s.colormap).toBe("viridis");
    expect(s.hideZeroExpression).toBe(true);
    expect(s.alphaScale).toBe(0.8);
  });

  it("leaves the preset's own values out of the link, and keeps changes", () => {
    const store = createVisualizationStoreInstance();
    const preset = technologyPresetOf(bins);

    applyTechnologyPreset(preset, store.getState());
    expect(encodeCellVizState(store.getState(), preset)).toBeNull();

    // Back to the built-in defaults: still a change from the preset
    store.getState().setColormap("bwr");
    store.getState().setHideZeroExpression(false);
    store.getState().setAlphaScale(1);
    const decoded = decodeCellVizState(
      encodeCellVizState(store.getState(), preset)!,
    );

    expect(decoded).toMatchObject({ cm: "bwr", hz: false, op: 1 });
  });

  it("without a preset, writes only what differs from the defaults", () => {
    const store = createVisualizationStoreInstance();

    expect(encodeCellVizState(store.getState())).toBeNull();
    store.getState().setAlphaScale(0.4);
    expect(
      decodeCellVizState(encodeCellVizState(store.getState())!),
    ).toEqual({ op: 0.4 });
  });
});
