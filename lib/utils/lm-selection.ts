import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type {
  LabelledMoleculeVisualizationState,
  LmMenu,
} from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { LM_MENUS } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

/**
 * Everything that has to travel together for a selection to look the same on
 * another embryo. Colours are part of it: carrying which genes are selected
 * without carrying what colour each one is gives two panels the same genes in
 * different colours, which is worse than not carrying them at all.
 */
export type LmCarryFields = Pick<
  LabelledMoleculeVisualizationState,
  "selections" | "hiddenValues" | "geneColorSlots" | "colorOverrides"
>;

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
 */
export function carrySelections(
  source: LmCarryFields,
  dataset: StandardizedDataset | null,
): LmCarryFields {
  const selections = {} as LmCarryFields["selections"];
  const hiddenValues = {} as LmCarryFields["hiddenValues"];

  for (const menu of LM_MENUS) {
    const vocab = vocabulary(dataset, menu);

    if (!vocab) {
      selections[menu] = source.selections[menu];
      hiddenValues[menu] = source.hiddenValues[menu];
      continue;
    }

    selections[menu] = new Set(
      [...source.selections[menu]].filter((v) => vocab.has(v)),
    );
    hiddenValues[menu] = new Set(
      [...source.hiddenValues[menu]].filter((v) => vocab.has(v)),
    );
  }

  return {
    selections,
    hiddenValues,
    // Surviving genes keep their slot. Rebuilding the slots would recolour
    // genes that carried over perfectly well — and across a split it would
    // show the same gene in two colours.
    geneColorSlots: new Map(
      [...source.geneColorSlots].filter(([gene]) => selections.gene.has(gene)),
    ),
    // Copied whole, not pruned: an override for a value this embryo lacks is
    // inert, and keeping it means the colour returns if you go back.
    colorOverrides: source.colorOverrides,
  };
}

/** Whether two carried states would draw identically. */
export function sameCarry(a: LmCarryFields, b: LmCarryFields): boolean {
  const setsEqual = (x: Set<string>, y: Set<string>) =>
    x.size === y.size && [...x].every((v) => y.has(v));

  for (const menu of LM_MENUS) {
    if (!setsEqual(a.selections[menu], b.selections[menu])) return false;
    if (!setsEqual(a.hiddenValues[menu], b.hiddenValues[menu])) return false;
  }

  if (a.geneColorSlots.size !== b.geneColorSlots.size) return false;
  for (const [gene, slot] of a.geneColorSlots) {
    if (b.geneColorSlots.get(gene) !== slot) return false;
  }

  return LM_MENUS.every((menu) => {
    const x = a.colorOverrides[menu];
    const y = b.colorOverrides[menu];
    const kx = Object.keys(x);

    return (
      kx.length === Object.keys(y).length && kx.every((k) => x[k] === y[k])
    );
  });
}
