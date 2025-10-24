/**
 * Shared utility for filtering control/unassigned genes
 * Used by both single cell and single molecule data processing
 */

/**
 * Determines if a gene should be filtered out based on its name or feature type
 * Filters control probes, negative controls, unassigned genes, etc.
 */
export function shouldFilterGene(gene: string, featureType?: string): boolean {
  const g = gene.toLowerCase();
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
