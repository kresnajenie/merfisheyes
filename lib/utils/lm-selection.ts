import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type {
  LabelledMoleculeVisualizationState,
  LmMenu,
} from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { LM_MENUS } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

/**
 * Values of `menu` the dataset actually has.
 *
 * Null means "cannot tell yet" — a column that has not finished lazy-loading
 * reports no values, and treating that as an empty vocabulary would drop every
 * selection rather than keep the ones the dataset does have.
 */
function vocabulary(
  dataset: StandardizedDataset | null,
  menu: LmMenu,
): Set<string> | null {
  const col = dataset?.clusters?.find((c) => c.column === menu);

  if (!col?.uniqueValues?.length) return null;

  return new Set(col.uniqueValues);
}

/**
 * Carry a selection onto another embryo, keeping only values it has.
 *
 * Used when switching datasets in place and when syncing a split panel. The
 * intersection is what lets one rule cover the whole series: measured pairwise
 * across all 45 spiralia datasets, gene labels overlap 96% within a stage and
 * 94% across, while cell overlaps 80% within a stage and 0% across, and domain
 * 67% and 0%. No cell or domain name is common to all 45 — a 2-cell embryo has
 * AB/CD where a 24-cell has 1a1 — so genes carry everywhere and the other two
 * carry within a developmental stage and simply fall away across one.
 *
 * Returns null when nothing would change, so a caller can skip the write.
 */
export function intersectSelectionsWithDataset(
  state: Pick<
    LabelledMoleculeVisualizationState,
    "selections" | "hiddenValues" | "geneColorSlots"
  >,
  dataset: StandardizedDataset | null,
): Partial<LabelledMoleculeVisualizationState> | null {
  const selections = {} as Record<LmMenu, Set<string>>;
  const hiddenValues = {} as Record<LmMenu, Set<string>>;
  let changed = false;

  for (const menu of LM_MENUS) {
    const vocab = vocabulary(dataset, menu);
    const current = state.selections[menu];
    const hidden = state.hiddenValues[menu];

    if (!vocab) {
      selections[menu] = current;
      hiddenValues[menu] = hidden;
      continue;
    }

    const keptSel = new Set([...current].filter((v) => vocab.has(v)));
    const keptHid = new Set([...hidden].filter((v) => vocab.has(v)));

    if (keptSel.size !== current.size || keptHid.size !== hidden.size) {
      changed = true;
    }
    selections[menu] = keptSel;
    hiddenValues[menu] = keptHid;
  }

  if (!changed) return null;

  // Keep each surviving gene's colour slot. Rebuilding them would recolour
  // genes that carried over perfectly well, which is the opposite of what
  // carrying a selection across embryos is for.
  const geneColorSlots = new Map(
    [...state.geneColorSlots].filter(([gene]) => selections.gene.has(gene)),
  );

  return { selections, hiddenValues, geneColorSlots };
}
