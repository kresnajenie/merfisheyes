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
