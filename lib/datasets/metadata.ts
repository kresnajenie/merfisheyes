/**
 * Catalog-style display metadata shared by owned Datasets, Projects, and the
 * CatalogDataset records they're submitted into. One source of truth for the
 * editable field list + a sanitizer, so the account editor, the project editor,
 * and the community-submit flow all agree on the shape.
 */

export interface DatasetMetadata {
  description: string | null;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  tags: string[];
  externalLink: string | null;
  publicationLink: string | null;
}

const STRING_FIELDS: {
  key: keyof Omit<DatasetMetadata, "tags">;
  max: number;
}[] = [
  { key: "description", max: 5000 },
  { key: "species", max: 200 },
  { key: "disease", max: 200 },
  { key: "tissue", max: 200 },
  { key: "platform", max: 200 },
  { key: "institute", max: 500 },
  { key: "externalLink", max: 1000 },
  { key: "publicationLink", max: 1000 },
];

/**
 * Pull the metadata fields out of an arbitrary request body, trimming strings
 * (empty → null) and normalizing tags. Only keys present in `body` are
 * returned, so this can drive a partial (PATCH) update. Unknown keys ignored.
 */
export function sanitizeMetadataPatch(
  body: Record<string, unknown>,
): Partial<DatasetMetadata> {
  const out: Partial<DatasetMetadata> = {};

  for (const { key, max } of STRING_FIELDS) {
    if (key in body) {
      const raw = body[key];

      if (raw === null) {
        out[key] = null;
      } else if (typeof raw === "string") {
        const trimmed = raw.trim().slice(0, max);

        out[key] = trimmed || null;
      }
    }
  }

  if ("tags" in body) {
    const raw = body.tags;

    if (Array.isArray(raw)) {
      out.tags = raw
        .filter((t): t is string => typeof t === "string")
        .map((t) => t.trim())
        .filter(Boolean)
        .slice(0, 30);
    } else if (typeof raw === "string") {
      // Accept a comma-separated string too (matches the admin form).
      out.tags = raw
        .split(",")
        .map((t) => t.trim())
        .filter(Boolean)
        .slice(0, 30);
    }
  }

  return out;
}

/**
 * Flexible display-metadata bag: the dynamic fields the Explore detail page
 * renders (mirrors CatalogDataset.metadata). These are stored in the `metadata`
 * Json column on Project / Dataset / CatalogDataset, keyed exactly like this so
 * DatasetDetailPage's METADATA_LABELS pick up the human labels. Single source of
 * truth for the editor form + the sanitizer.
 */
export const METADATA_BAG_FIELDS: {
  key: string;
  label: string;
  max: number;
  /** Render as a multi-line textarea (e.g. citation). */
  multiline?: boolean;
}[] = [
  { key: "authors", label: "Authors", max: 2000 },
  { key: "investigator", label: "Investigator (PI)", max: 300 },
  { key: "coInvestigators", label: "Co-investigators", max: 1000 },
  { key: "institution", label: "Institution", max: 500 },
  { key: "age", label: "Age", max: 200 },
  { key: "sex", label: "Sex", max: 100 },
  { key: "genotype", label: "Genotype", max: 300 },
  { key: "technique", label: "Technique", max: 300 },
  { key: "funding", label: "Funding", max: 1000 },
  { key: "publicationYear", label: "Publication year", max: 20 },
  { key: "license", label: "License", max: 200 },
  { key: "citation", label: "Citation", max: 3000, multiline: true },
];

const BAG_KEYS = new Set(METADATA_BAG_FIELDS.map((f) => f.key));
const BAG_MAX = new Map(METADATA_BAG_FIELDS.map((f) => [f.key, f.max]));

/**
 * Clean a metadata-bag object: keep only known keys, trim, drop empties, cap
 * length. Returns a plain object suitable for the `metadata` Json column.
 */
export function sanitizeMetadataBag(raw: unknown): Record<string, string> {
  const out: Record<string, string> = {};

  if (!raw || typeof raw !== "object") return out;

  for (const [key, value] of Object.entries(raw as Record<string, unknown>)) {
    if (!BAG_KEYS.has(key) || typeof value !== "string") continue;
    const trimmed = value.trim().slice(0, BAG_MAX.get(key) ?? 500);

    if (trimmed) out[key] = trimmed;
  }

  return out;
}
