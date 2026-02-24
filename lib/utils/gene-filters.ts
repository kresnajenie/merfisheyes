/**
 * Shared utility for filtering control/unassigned genes
 * Used by both single cell and single molecule data processing
 */

/** Suffix appended to gene names for unassigned molecules */
export const UNASSIGNED_SUFFIX = " (unassigned)";

/**
 * Determines if a cell_id value represents an unassigned molecule.
 * Returns true for null, undefined, empty string, "UNASSIGNED" (case-insensitive),
 * or negative numbers.
 */
export function isUnassignedCell(cellId: unknown): boolean {
  if (cellId === null || cellId === undefined) return true;

  if (typeof cellId === "string") {
    const trimmed = cellId.trim();

    if (trimmed === "") return true;
    if (trimmed.toLowerCase() === "unassigned") return true;

    // Check if string represents a negative number
    const num = Number(trimmed);

    if (!isNaN(num) && num < 0) return true;

    return false;
  }

  if (typeof cellId === "number") {
    return cellId < 0;
  }

  return false;
}

/**
 * Determines if a gene should be filtered out based on its name or feature type
 * Filters control probes, negative controls, unassigned genes, etc.
 */
export function shouldFilterGene(gene: string, featureType?: string): boolean {
  // Strip unassigned suffix before checking patterns
  const baseName = gene.endsWith(UNASSIGNED_SUFFIX)
    ? gene.slice(0, -UNASSIGNED_SUFFIX.length)
    : gene;

  const g = baseName.toLowerCase();
  const t = (featureType || "").toLowerCase();
  const patterns = [
    /negative[\s_-]*control/,
    /neg[\s_-]*ctrl/,
    /unassigned/,
    /deprecated/,
    /codeword/,
    /blank/,
    /negcontrol/,
  ];

  return patterns.some((rx) => rx.test(g) || rx.test(t));
}
