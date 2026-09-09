import type { StandardizedDataset } from "../StandardizedDataset";

/**
 * Which labelled-molecule dataset each panel currently has open.
 *
 * The cell and single-molecule viewers keep their datasets in zustand stores,
 * so their sync hooks can reach across panels through those. The LM viewer
 * holds its dataset in component state — it loads exactly one, and a store
 * bought nothing — so selection sync needs somewhere to look up the *other*
 * panel's dataset to intersect against its vocabulary.
 *
 * Deliberately not a store: nothing renders off it. It is read inside a
 * subscription callback at the moment a selection changes, so a plain map is
 * enough and avoids a re-render on every dataset swap.
 */
const datasets = new Map<string, StandardizedDataset | null>();

/** Panel id for the left/only viewer, which has no PanelContext. */
export const LEFT_PANEL = "left";

export function setLmDataset(
  panelId: string | null,
  dataset: StandardizedDataset | null,
) {
  datasets.set(panelId ?? LEFT_PANEL, dataset);
}

export function getLmDataset(
  panelId: string | null,
): StandardizedDataset | null {
  return datasets.get(panelId ?? LEFT_PANEL) ?? null;
}

export function clearLmDataset(panelId: string | null) {
  datasets.delete(panelId ?? LEFT_PANEL);
}
