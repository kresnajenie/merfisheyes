/**
 * Picks up to 3 genes to auto-select when a single molecule dataset loads.
 *
 * Priority order:
 * 1. Priority genes (Aqp4, Gad2, Cldn5) — case-insensitive match
 * 2. Remaining slots filled with genes that don't start with a digit
 * 3. If all genes start with a digit, fall back to the first 3
 */

const PRIORITY_GENES = ["Igfbpl1", "Drd1", "Th"];
const AUTO_SELECT_COUNT = 3;

export function pickDefaultGenes(uniqueGenes: string[]): string[] {
  if (uniqueGenes.length === 0) return [];

  const selected: string[] = [];

  // 1. Add priority genes that exist (case-insensitive)
  for (const priority of PRIORITY_GENES) {
    if (selected.length >= AUTO_SELECT_COUNT) break;

    const match = uniqueGenes.find(
      (g) => g.toLowerCase() === priority.toLowerCase(),
    );

    if (match && !selected.includes(match)) {
      selected.push(match);
    }
  }

  // 2. Fill remaining with non-numeric-start genes
  if (selected.length < AUTO_SELECT_COUNT) {
    for (const gene of uniqueGenes) {
      if (selected.length >= AUTO_SELECT_COUNT) break;
      if (selected.includes(gene)) continue;
      if (!/^\d/.test(gene)) {
        selected.push(gene);
      }
    }
  }

  // 3. If still not enough (all start with digits), take first available
  if (selected.length < AUTO_SELECT_COUNT) {
    for (const gene of uniqueGenes) {
      if (selected.length >= AUTO_SELECT_COUNT) break;
      if (!selected.includes(gene)) {
        selected.push(gene);
      }
    }
  }

  return selected;
}
