import { describe, expect, it } from "vitest";

import {
  decodeSMVizState,
  encodeSMVizState,
  toBase64Url,
} from "@/lib/utils/url-visualization-state";

const gene = (g: string) => ({
  gene: g,
  color: "#00d8ff",
  localScale: 1,
  showAssigned: true,
  showUnassigned: true,
  unassignedColor: "#00d8ff",
  unassignedLocalScale: 1,
  colorSynced: true,
});

const base = {
  selectedGenes: new Map([["Gad1", gene("Gad1")]]),
  geneDataCache: new Map([["Gad1", gene("Gad1")]]),
  globalScale: 1,
  viewMode: "2D" as const,
  showAssigned: true,
  showUnassigned: true,
};

describe("molecule viewer link state", () => {
  it("keeps molecule opacity in the link, and leaves it out at 1", () => {
    const at = decodeSMVizState(
      encodeSMVizState({ ...base, moleculeOpacity: 0.35 })!,
    );

    expect(at?.op).toBe(0.35);
    expect(
      decodeSMVizState(encodeSMVizState({ ...base, moleculeOpacity: 1 })!)?.op,
    ).toBeUndefined();
  });

  it("reads a link written by an agent (merfisheyes-mcp meyes_open_view)", () => {
    const state = { genes: [["Slc17a7", "#00d8ff", 1, true]], op: 0.35 };

    expect(decodeSMVizState(toBase64Url(state))).toEqual(state);
  });

  it("refuses an opacity that is not a number", () => {
    expect(decodeSMVizState(toBase64Url({ op: "faint" }))).toBeNull();
  });
});
