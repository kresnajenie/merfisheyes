/**
 * Human-readable labels for the free-form `metadata` bag carried by catalog
 * rows, projects and datasets. Shared by the detail page and the card's info
 * modal so the two can't drift; unknown keys fall through to the raw key.
 */
export const METADATA_LABELS: Record<string, string> = {
  authors: "Authors",
  investigator: "Investigator",
  institution: "Institution",
  coInvestigators: "Co-Investigators",
  funding: "Funding",
  publicationYear: "Year",
  license: "License",
  age: "Age",
  sex: "Sex",
  genotype: "Genotype",
  technique: "Technique",
  citation: "Citation",
};

/** Label + rendered value for each non-empty entry of a metadata bag. */
export function metadataEntries(
  metadata: Record<string, unknown> | null | undefined,
): { key: string; label: string; value: string }[] {
  return Object.entries(metadata ?? {})
    .filter(([, value]) => Boolean(value))
    .map(([key, value]) => ({
      key,
      label: METADATA_LABELS[key] || key,
      value: Array.isArray(value) ? value.join(", ") : String(value),
    }));
}
