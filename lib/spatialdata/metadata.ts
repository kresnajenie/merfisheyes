/**
 * A store's dataset metadata, read from its root `zarr.json` attributes
 * (docs/specs/metadata.md §4; schema docs/specs/schemas/metadata.schema.json).
 * Up to four blocks:
 *
 *   - `metadata_general`: SpatialData-DB's, labels only (Title Case keys)
 *   - `metadata_<technology>`: SpatialData-DB's per-technology block
 *   - `meyes_metadata`: ours, the same facts with ontology ids, under the
 *     query index's column names
 *   - `meyes_provenance`: ours, who wrote it and what validated it
 *
 * A SpatialData-DB store carries only the first two (and their
 * `sddb_provenance`). Where both state a fact, ours wins (it has the id);
 * theirs fills what ours lacks. The read is the root attributes the store
 * was opened with: nothing is fetched here.
 */
import { CONVENTIONS } from "./conventions";

export type MetadataValue = string | number | boolean;

/** One fact about the dataset, as the info panel shows it. */
export interface MetadataFact {
  /** The index column it is (`organism`, `n_obs`, ...). */
  key: string;
  label: string;
  value: MetadataValue;
  /** Its ontology id (`UBERON:0000955`), when the store gives one. */
  id: string | null;
  /** Where the value or id can be looked up. */
  href: string | null;
  /** Which block said it: ours, or SpatialData-DB's labels. */
  from: "meyes" | "general";
  /** Where the producer says the value came from (`bil-api`, `curator`). */
  source: string | null;
}

export interface MetadataGroup {
  title: string;
  facts: MetadataFact[];
}

/** One section of a store holding several (`meyes_sections`). */
export interface MetadataSection {
  id: string;
  sampleId: string | null;
  donorId: string | null;
  timepoint: string | null;
  sex: string | null;
  nObs: number | null;
  nTranscripts: number | null;
  matrixHash: string | null;
}

export type ValidationStatus = "pass" | "fail" | "not_validated";

export interface MetadataValidation {
  status: ValidationStatus;
  at: string | null;
  /** Ontology → the version checked against. */
  sources: Record<string, string>;
  errors: string[];
}

export interface MetadataProvenance {
  /** Who wrote the metadata: our pipeline, or SpatialData-DB's converter. */
  writer: "merfisheyes" | "spatialdata-db";
  canonicalKey: string | null;
  /** Our pipeline's git sha, or SpatialData-DB's container version. */
  version: string | null;
  input: string | null;
  /** ISO 8601. */
  writtenAt: string | null;
}

/** A block of the root attributes, every key as written. */
export interface MetadataBlock {
  name: string;
  entries: [string, unknown][];
}

export interface StoreMetadata {
  title: string | null;
  description: string | null;
  groups: MetadataGroup[];
  sections: MetadataSection[];
  /** `mh1:<64 hex>`: equal hashes are the same counts (metadata.md §7.7). */
  matrixHash: string | null;
  validation: MetadataValidation | null;
  provenance: MetadataProvenance | null;
  /** Every block, all its keys: nothing the panel does not know is dropped. */
  blocks: MetadataBlock[];
}

/** Our block and its provenance: a rule of our format, not a producer's name. */
const OURS = "meyes_metadata";
const OUR_PROVENANCE = "meyes_provenance";
/** `metadata_<technology>` (the schema's patternProperties). */
const TECHNOLOGY_BLOCK = /^metadata_[a-z0-9_]+$/;

type Rec = Record<string, unknown>;

const isRecord = (v: unknown): v is Rec =>
  typeof v === "object" && v !== null && !Array.isArray(v);

const str = (v: unknown): string | null =>
  typeof v === "string" && v.trim() !== "" ? v : null;

const num = (v: unknown): number | null =>
  typeof v === "number" && Number.isFinite(v) ? v : null;

const present = (v: unknown): v is MetadataValue =>
  (typeof v === "string" && v.trim() !== "") ||
  (typeof v === "number" && Number.isFinite(v)) ||
  typeof v === "boolean";

/** An ontology term's page: OLS knows every source we pin. */
export function ontologyHref(id: string): string | null {
  return /^[A-Za-z][A-Za-z0-9]*:[A-Za-z0-9_]+$/.test(id)
    ? `https://www.ebi.ac.uk/ols4/search?q=${encodeURIComponent(id)}`
    : null;
}

const httpHref = (v: MetadataValue) =>
  typeof v === "string" && /^https?:\/\//.test(v) ? v : null;

const doiHref = (v: MetadataValue) =>
  typeof v === "string"
    ? `https://doi.org/${v.replace(/^(https?:\/\/(dx\.)?doi\.org\/|doi:)/i, "")}`
    : null;

interface FactSpec {
  key: string;
  label: string;
  /** Our column holding its ontology id. */
  id?: string;
  href?: (v: MetadataValue) => string | null;
}

/** The facts shown up front, grouped, in order; the rest is under all fields. */
const GROUPS: { title: string; facts: FactSpec[] }[] = [
  {
    title: "Sample",
    facts: [
      { key: "organism", label: "Organism", id: "organism_id" },
      { key: "tissue", label: "Tissue", id: "tissue_id" },
      { key: "disease", label: "Disease", id: "disease_id" },
      {
        key: "development_stage",
        label: "Development stage",
        id: "meyes_development_stage_id",
      },
      { key: "meyes_sex", label: "Sex", id: "meyes_sex_id" },
      { key: "timepoint", label: "Timepoint" },
      { key: "genotype", label: "Genotype" },
      { key: "condition", label: "Condition" },
      { key: "treatment", label: "Treatment" },
      { key: "sample_id", label: "Sample ID" },
      { key: "donor_id", label: "Donor" },
      { key: "preservation_method", label: "Preservation" },
    ],
  },
  {
    title: "Technology",
    facts: [
      { key: "assay", label: "Assay" },
      { key: "technology", label: "Technology" },
      { key: "vendor", label: "Vendor" },
      { key: "product", label: "Product" },
      { key: "modality", label: "Modality" },
      { key: "chemistry_version", label: "Chemistry" },
      { key: "segmentation_method", label: "Segmentation" },
    ],
  },
  {
    title: "Data",
    facts: [
      { key: "n_obs", label: "Cells" },
      { key: "n_features", label: "Features" },
      { key: "n_transcripts", label: "Transcripts" },
      { key: "meyes_n_sections", label: "Sections" },
      { key: "meyes_dims", label: "Dimensions" },
      { key: "default_table", label: "Default table" },
    ],
  },
  {
    title: "Source",
    facts: [
      { key: "license", label: "License" },
      { key: "doi", label: "DOI", href: doiHref },
      { key: "publication_date", label: "Published" },
      { key: "dataset_url", label: "Source page", href: httpHref },
      { key: "canonical_key", label: "Key" },
      { key: "study_id", label: "Study" },
      { key: "collection_name", label: "Collection" },
    ],
  },
];

/** A licence id or label → its label and page, from the pinned vocabulary. */
function license(value: MetadataValue): {
  label: MetadataValue;
  href: string | null;
  id: string | null;
} {
  const licenses = CONVENTIONS.metadata.licenses;
  const byId = typeof value === "string" ? licenses[value] : undefined;

  if (byId) return { label: byId.label, href: byId.url, id: value as string };
  const byLabel = Object.entries(licenses).find(([, l]) => l.label === value);

  return byLabel
    ? { label: value, href: byLabel[1].url, id: byLabel[0] }
    : { label: value, href: null, id: null };
}

function readFact(
  spec: FactSpec,
  ours: Rec,
  general: Rec,
  sources: Rec,
): MetadataFact | null {
  const generalKey = CONVENTIONS.metadata.sddb.generalKeys[spec.key];
  const mine = ours[spec.key];
  const theirs = generalKey ? general[generalKey] : undefined;
  const from = present(mine) ? "meyes" : present(theirs) ? "general" : null;

  if (!from) return null;
  let value = (from === "meyes" ? mine : theirs) as MetadataValue;
  let id = from === "meyes" && spec.id ? str(ours[spec.id]) : null;
  let href = id ? ontologyHref(id) : (spec.href?.(value) ?? null);

  if (spec.key === "license") {
    const l = license(value);

    value = l.label;
    id = l.id;
    href = l.href;
  }
  if (spec.key === "timepoint" && from === "meyes") {
    const kind = str(ours.timepoint_kind);

    if (kind) value = `${value} (${kind})`;
  }

  return {
    key: spec.key,
    label: spec.label,
    value,
    id,
    href,
    from,
    source: str(sources[spec.key]),
  };
}

function readSections(v: unknown): MetadataSection[] {
  if (!Array.isArray(v)) return [];

  return v.filter(isRecord).flatMap((s) => {
    const id = str(s.id) ?? (typeof s.id === "number" ? String(s.id) : null);

    return id === null
      ? []
      : [
          {
            id,
            sampleId: str(s.sample_id),
            donorId: str(s.donor_id),
            timepoint: str(s.timepoint),
            sex: str(s.meyes_sex),
            nObs: num(s.n_obs),
            nTranscripts: num(s.n_transcripts),
            matrixHash: str(s.meyes_matrix_hash),
          },
        ];
  });
}

const stringRecord = (v: unknown): Record<string, string> =>
  isRecord(v)
    ? Object.fromEntries(
        Object.entries(v).filter(
          (e): e is [string, string] => typeof e[1] === "string",
        ),
      )
    : {};

function readValidation(provenance: Rec): MetadataValidation | null {
  const v = provenance.validation;

  if (!isRecord(v)) return null;
  const status = v.status;

  return {
    status: status === "pass" || status === "fail" ? status : "not_validated",
    at: str(v.at),
    sources: stringRecord(v.sources),
    errors: Array.isArray(v.errors)
      ? v.errors.filter((e): e is string => typeof e === "string")
      : [],
  };
}

function readProvenance(
  ours: Rec | null,
  sddb: Rec | null,
): MetadataProvenance | null {
  if (ours) {
    return {
      writer: "merfisheyes",
      canonicalKey: str(ours.canonical_key),
      version: str(ours.pipeline),
      input: str(ours.input),
      writtenAt: str(ours.written_at),
    };
  }
  if (sddb) {
    return {
      writer: "spatialdata-db",
      canonicalKey: str(sddb.canonical_key),
      version: str(sddb.container_version),
      input: null,
      writtenAt: str(sddb.converted_at),
    };
  }

  return null;
}

/**
 * The store's metadata from its root attributes, or null when it carries
 * none (a store we did not write and SpatialData-DB did not convert).
 * Never throws: any block may be missing, any key absent or null.
 */
export function readStoreMetadata(rootAttrs: unknown): StoreMetadata | null {
  if (!isRecord(rootAttrs)) return null;
  const generalName = CONVENTIONS.tables.curatedDefaultTable.attribute;
  const sddbName = CONVENTIONS.metadata.sddb.provenanceAttribute;
  const block = (name: string) =>
    isRecord(rootAttrs[name]) ? (rootAttrs[name] as Rec) : null;
  const general = block(generalName);
  const ours = block(OURS);
  const ourProvenance = block(OUR_PROVENANCE);
  const sddbProvenance = block(sddbName);
  const technology = Object.keys(rootAttrs)
    .filter((k) => k !== generalName && TECHNOLOGY_BLOCK.test(k) && block(k))
    .sort();

  if (
    !general &&
    !ours &&
    !ourProvenance &&
    !sddbProvenance &&
    technology.length === 0
  ) {
    return null;
  }
  const sources = isRecord(ourProvenance?.sources) ? ourProvenance.sources : {};
  const groups = GROUPS.map(({ title, facts }) => ({
    title,
    facts: facts.flatMap(
      (f) => readFact(f, ours ?? {}, general ?? {}, sources) ?? [],
    ),
  })).filter((g) => g.facts.length > 0);
  const blocks: MetadataBlock[] = [
    [OURS, ours],
    [generalName, general],
    ...technology.map((k) => [k, block(k)] as const),
    [OUR_PROVENANCE, ourProvenance],
    [sddbName, sddbProvenance],
  ].flatMap(([name, value]) =>
    value ? [{ name: name as string, entries: Object.entries(value) }] : [],
  );

  return {
    title: str(ours?.title),
    description: str(ours?.meyes_description),
    groups,
    sections: readSections(ours?.meyes_sections),
    matrixHash: str(ours?.meyes_matrix_hash),
    validation: ourProvenance ? readValidation(ourProvenance) : null,
    provenance: readProvenance(ourProvenance, sddbProvenance),
    blocks,
  };
}
