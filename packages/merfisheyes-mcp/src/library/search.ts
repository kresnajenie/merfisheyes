import type { GeneStat, LibraryIndex, Row } from "./data.ts";

/**
 * Search over the published index (MCP.md pillar B item 3): structured
 * filters over every live dataset, answered as a count, a breakdown by the
 * facets that matter, notes on what was resolved or left out, and a page of
 * rows. Pure: an index and filters in, a result out.
 *
 * Rules the answers keep:
 * - A tissue or disease includes its descendants ("brain" finds the
 *   hippocampus) through `ontology_rollup.json` and the grouping columns.
 * - A filter on a field few rows record (age, sex, genotype, condition,
 *   treatment) leaves out the rows that do not record it, and says how many.
 * - Expression and measured genes come from the gene index, which covers our
 *   datasets only for now; the note says so.
 * - Re-releases of the same data (`same_section_relation = rerelease`)
 *   collapse to the newest release, the others listed with it.
 */

export interface GeneFilter {
  gene: string;
  /** Fraction of cells detecting it, 0–1; 0 means "measured". */
  min_fraction: number;
}

export interface Filters {
  organism?: string[];
  tissue?: string[];
  disease?: string[];
  technology?: string[];
  source?: string[];
  kind?: string[];
  metadata_from?: string[];
  modality?: string;
  sex?: string;
  age_days_min?: number;
  age_days_max?: number;
  genotype?: string;
  condition?: string;
  treatment?: string;
  text?: string;
  min_cells?: number;
  genes?: GeneFilter[];
  collapse_rereleases?: boolean;
}

export type Sort = "newest" | "cells" | "title";

export interface Hit {
  id: string;
  title: string | null;
  organism: string | null;
  tissue: string | null;
  tissue_organ: string | null;
  disease: string | null;
  technology: string | null;
  cells: number | null;
  obs_unit: string | null;
  genes_measured: number | null;
  timepoint: string | null;
  age_days: number | null;
  sex: string | null;
  genotype: string | null;
  condition: string | null;
  treatment: string | null;
  source: string | null;
  metadata_from: string | null;
  format: string | null;
  kind: string | null;
  study_id: string | null;
  doi: string | null;
  bil_code: string | null;
  page_url: string | null;
  viewer_url: string | null;
  also_published_as?: { id: string; release: string | null }[];
  gene_stats?: {
    gene: string;
    fraction_detected: number | null;
    total_counts: number | null;
  }[];
}

export interface Breakdown {
  [facet: string]: { value: string; count: number }[];
}

export interface SearchResult {
  release: string;
  total: number;
  notes: string[];
  breakdown: Breakdown;
  rows: Row[];
  geneStats: Map<string, GeneStat[]>;
  collapsed: Map<string, Row[]>;
}

/** Thrown for a filter value the index has no term for; the message names
 * the closest ones. */
export class UnknownValue extends Error {}

const lc = (v: unknown) =>
  v === null || v === undefined ? "" : String(v).toLowerCase().trim();
const str = (v: unknown) => (v === null || v === undefined ? null : String(v));
const num = (v: unknown) => (typeof v === "number" ? v : null);

const ORGANISM_ALIASES: Record<string, string> = {
  "mus musculus": "mouse",
  "homo sapiens": "human",
  "macaca mulatta": "macaque",
  "rhesus macaque": "macaque",
  "macaca fascicularis": "crab-eating macaque",
  marmoset: "white-tufted-ear marmoset",
  "callithrix jacchus": "white-tufted-ear marmoset",
  "common marmoset": "white-tufted-ear marmoset",
  "rattus norvegicus": "rat",
  "danio rerio": "zebrafish",
  mice: "mouse",
  humans: "human",
  macaques: "macaque",
  rats: "rat",
  "rhesus monkey": "macaque",
  monkey: "macaque",
};

const DISEASE_ALIASES: Record<string, string> = {
  healthy: "normal",
  control: "normal",
  "no disease": "normal",
};

/** Edit distance, for "did you mean" on a misspelt value. */
function distance(a: string, b: string): number {
  let prev = Array.from({ length: b.length + 1 }, (_, j) => j);

  for (let i = 1; i <= a.length; i++) {
    const cur = [i];

    for (let j = 1; j <= b.length; j++)
      cur[j] = Math.min(
        prev[j] + 1,
        cur[j - 1] + 1,
        prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1),
      );
    prev = cur;
  }

  return prev[b.length];
}

/** The known values nearest a query: containing it (most common first),
 * else within a few edits. */
function closest(names: Map<string, number>, q: string, n = 8): string[] {
  const sub = [...names.entries()].filter(
    ([k]) => k.includes(q) || (q.length > 3 && q.includes(k)),
  );

  if (sub.length)
    return sub
      .sort((a, b) => b[1] - a[1])
      .slice(0, n)
      .map(([k]) => k);
  const max = Math.max(2, Math.floor(q.length / 3));
  const near = (k: string) =>
    Math.min(distance(q, k), ...k.split(/[\s,/-]+/).map((w) => distance(q, w)));

  return [...names.keys()]
    .map((k) => [k, near(k)] as const)
    .filter(([, d]) => d <= max)
    .sort((a, b) => a[1] - b[1])
    .slice(0, n)
    .map(([k]) => k);
}

function countBy(rows: Row[], get: (r: Row) => unknown[]): Map<string, number> {
  const m = new Map<string, number>();

  for (const r of rows)
    for (const v of get(r)) {
      const k = lc(v);

      if (k) m.set(k, (m.get(k) ?? 0) + 1);
    }

  return m;
}

/**
 * A tissue or disease, by name or ontology id, with everything below it.
 * Returns a row test and what it was resolved to.
 */
function termMatcher(
  index: LibraryIndex,
  kind: "tissue" | "disease",
  query: string,
): { test: (r: Row) => boolean; said: string } {
  let q = lc(query);

  if (kind === "disease") q = DISEASE_ALIASES[q] ?? q;
  const prefix = kind === "tissue" ? "UBERON" : "MONDO";
  const groupCols =
    kind === "tissue"
      ? (["tissue_organ", "tissue_system"] as const)
      : (["disease_class", "disease_group"] as const);
  const ids = new Set<string>();
  const names = new Set<string>();

  if (/^[a-z]+:\d+$/.test(q)) ids.add(q.toUpperCase());
  else {
    for (const [id, t] of Object.entries(index.rollup))
      if (id.startsWith(prefix) && lc(t.name) === q) ids.add(id);
    for (const r of index.rows) {
      if (lc(r[kind]) === q) {
        names.add(q);
        if (r[`${kind}_id`]) ids.add(String(r[`${kind}_id`]));
      }
      for (const g of groupCols) if (lc(r[g]) === q) names.add(q);
    }
  }
  if (!ids.size && !names.size) {
    const known = countBy(index.rows, (r) => [
      r[kind],
      ...groupCols.map((g) => r[g]),
    ]);

    for (const [id, t] of Object.entries(index.rollup))
      if (id.startsWith(prefix) && !known.has(lc(t.name)))
        known.set(lc(t.name), 0);
    const near = closest(known, q);

    throw new UnknownValue(
      `No ${kind} named "${query}" in the index.${near.length ? ` Closest: ${near.join(", ")}.` : ""} Call meyes_index_values with column="${kind}" to list them.`,
    );
  }
  const idCol = `${kind}_id` as const;
  const test = (r: Row) => {
    if (names.has(lc(r[kind]))) return true;
    for (const g of groupCols) if (names.has(lc(r[g]))) return true;
    const id = str(r[idCol]);

    if (!id) return false;
    if (ids.has(id)) return true;
    const up = index.rollup[id]?.ancestors_present;

    return !!up && up.some((a) => ids.has(a));
  };
  const idList = [...ids].slice(0, 4).join(", ");

  return {
    test,
    said: `${kind} "${query}" → ${idList || [...names].join(", ")}${ids.size > 4 ? ", …" : ""}, with everything recorded below it`,
  };
}

function technologyMatcher(index: LibraryIndex, query: string) {
  const q = lc(query);
  const exact = index.technologies.filter(
    (t) =>
      lc(t.id) === q ||
      lc(t.label) === q ||
      (t.legacy_assay ?? []).some((a) => lc(a) === q),
  );
  const pick = exact.length
    ? exact
    : index.technologies.filter(
        (t) => lc(t.id).includes(q) || lc(t.label).includes(q),
      );
  const ids = new Set(pick.map((t) => t.id));

  // Rows can carry ids the vocabulary does not list yet.
  for (const r of index.rows)
    if (lc(r.technology) === q) ids.add(String(r.technology));
  if (!ids.size) {
    const known = countBy(index.rows, (r) => [r.technology]);

    throw new UnknownValue(
      `No technology "${query}" in the index. Known: ${[...known.keys()].join(", ")}.`,
    );
  }

  return {
    test: (r: Row) => ids.has(String(r.technology)),
    said: `technology "${query}" → ${[...ids].join(", ")}`,
  };
}

function organismMatcher(index: LibraryIndex, query: string) {
  let q = lc(query);

  q = ORGANISM_ALIASES[q] ?? q;
  const ok = index.rows.some(
    (r) => lc(r.organism) === q || lc(r.organism_id) === q,
  );

  if (!ok) {
    const known = countBy(index.rows, (r) => [r.organism]);
    const near = closest(known, q);

    throw new UnknownValue(
      `No organism "${query}" in the index.${near.length ? ` Closest: ${near.join(", ")}.` : ""} Known: ${[
        ...known.entries(),
      ]
        .sort((a, b) => b[1] - a[1])
        .slice(0, 20)
        .map(([k]) => k)
        .join(", ")}.`,
    );
  }

  return (r: Row) => lc(r.organism) === q || lc(r.organism_id) === q;
}

/** Version-aware comparison of release labels ("Space Ranger 4.0.1" >
 * "Space Ranger 3.0.0"). */
function compareRelease(a: Row, b: Row): number {
  const ra = lc(a.same_section_release);
  const rb = lc(b.same_section_release);

  if (ra !== rb) return ra.localeCompare(rb, undefined, { numeric: true });

  return lc(a.published_at).localeCompare(lc(b.published_at));
}

/** Contains, case-insensitive; every word of the query must appear. */
function contains(hay: string, query: string): boolean {
  return lc(query)
    .split(/\s+/)
    .filter(Boolean)
    .every((w) => hay.includes(w));
}

const TEXT_FIELDS: (keyof Row)[] = [
  "title",
  "meyes_description",
  "sample_title",
  "sample_id",
  "study_id",
  "collection_name",
  "doi",
  "tissue",
  "disease",
  "genotype",
  "condition",
  "treatment",
  "meyes_bil_code",
  "meyes_investigator",
  "meyes_institution",
  "canonical_key",
];

export const BREAKDOWN_FACETS = [
  "organism",
  "tissue_organ",
  "technology",
  "disease",
  "meyes_source",
  "meyes_metadata_from",
] as const;

function breakdown(rows: Row[], top = 8): Breakdown {
  const out: Breakdown = {};

  for (const f of BREAKDOWN_FACETS) {
    const m = new Map<string, number>();

    for (const r of rows) {
      const k = str(r[f]) ?? "(not recorded)";

      m.set(k, (m.get(k) ?? 0) + 1);
    }
    const sorted = [...m.entries()].sort((a, b) => b[1] - a[1]);
    const head = sorted
      .slice(0, top)
      .map(([value, count]) => ({ value, count }));
    const rest = sorted.slice(top).reduce((s, [, c]) => s + c, 0);

    if (rest)
      head.push({ value: `(${sorted.length - top} others)`, count: rest });
    out[f] = head;
  }

  return out;
}

/** The filters applied, in order, with notes on what each did. */
export async function search(
  index: LibraryIndex,
  f: Filters,
  sort: Sort = "newest",
): Promise<SearchResult> {
  const notes: string[] = [];
  let rows = index.rows;
  const keep = (test: (r: Row) => boolean) => {
    rows = rows.filter(test);
  };
  /** A filter on a field few rows record: say how many had none. */
  const sparse = (
    label: string,
    has: (r: Row) => boolean,
    test: (r: Row) => boolean,
  ) => {
    const missing = rows.filter((r) => !has(r)).length;

    rows = rows.filter((r) => has(r) && test(r));
    if (missing)
      notes.push(
        `${missing.toLocaleString("en-US")} datasets matching the other filters record no ${label} and were left out.`,
      );
  };
  const anyOf = (
    vals: string[] | undefined,
    mk: (v: string) => (r: Row) => boolean,
  ) => {
    if (!vals?.length) return;
    const tests = vals.map(mk);

    keep((r) => tests.some((t) => t(r)));
  };

  anyOf(f.organism, (v) => organismMatcher(index, v));
  anyOf(f.tissue, (v) => {
    const m = termMatcher(index, "tissue", v);

    notes.push(m.said + ".");

    return m.test;
  });
  anyOf(f.disease, (v) => {
    const m = termMatcher(index, "disease", v);

    notes.push(m.said + ".");

    return m.test;
  });
  anyOf(f.technology, (v) => {
    const m = technologyMatcher(index, v);

    notes.push(m.said + ".");

    return m.test;
  });
  anyOf(f.source, (v) => (r) => lc(r.meyes_source) === lc(v));
  anyOf(f.kind, (v) => (r) => lc(r.meyes_kind) === lc(v));
  anyOf(f.metadata_from, (v) => (r) => lc(r.meyes_metadata_from) === lc(v));
  if (f.modality) {
    const m = lc(f.modality);

    keep((r) => lc(r.modality).split("|").includes(m));
  }
  if (f.min_cells !== undefined) {
    const min = f.min_cells;

    keep((r) => (num(r.n_obs) ?? 0) >= min);
  }
  if (f.text) {
    const t = f.text;

    keep((r) => contains(TEXT_FIELDS.map((k) => lc(r[k])).join(" "), t));
  }
  if (f.age_days_min !== undefined || f.age_days_max !== undefined) {
    const lo = f.age_days_min ?? -Infinity;
    const hi = f.age_days_max ?? Infinity;

    sparse(
      "age in days",
      (r) => lc(r.timepoint_kind) === "age" && num(r.timepoint_days) !== null,
      (r) => {
        const d = num(r.timepoint_days)!;

        return d >= lo && d <= hi;
      },
    );
  }
  if (f.sex) {
    const s = lc(f.sex);

    sparse(
      "sex",
      (r) => !!r.meyes_sex,
      (r) => lc(r.meyes_sex) === s,
    );
  }
  for (const k of ["genotype", "condition", "treatment"] as const) {
    const v = f[k];

    if (v)
      sparse(
        k,
        (r) => !!r[k],
        (r) => contains(lc(r[k]), v),
      );
  }

  const geneStats = new Map<string, GeneStat[]>();

  if (f.genes?.length) {
    // One range read per gene, all at once.
    const [gi, ...found] = await Promise.all([
      index.geneIndex(),
      ...f.genes.map((g) => index.measured(g.gene)),
    ]);

    for (const [i, g] of f.genes.entries()) {
      const m = found[i];
      const stats = new Map(
        (gi.bySymbol.get(lc(g.gene)) ?? []).map((s) => [s.meyes_id, s]),
      );

      if (!m.entities.length)
        notes.push(
          `No dataset in the index measures "${g.gene}" (symbols and Ensembl ids are matched whole, case-insensitive).`,
        );
      rows = rows.filter((r) => m.rows.has(r.row));
      if (g.min_fraction > 0) {
        const unknown = rows.filter(
          (r) => stats.get(String(r.meyes_id))?.fraction_obs_detected == null,
        ).length;

        rows = rows.filter(
          (r) =>
            (stats.get(String(r.meyes_id))?.fraction_obs_detected ?? -1) >=
            g.min_fraction,
        );
        notes.push(
          `"${g.gene}" expressed: detected in ≥ ${(g.min_fraction * 100).toFixed(1)} % of cells (fraction_obs_detected).` +
            (unknown
              ? ` ${unknown.toLocaleString("en-US")} datasets that measure it have no per-dataset counts in the index yet (SpatialData-DB's, and ours with no store) and were left out.`
              : ""),
        );
      } else
        notes.push(
          `"${g.gene}" measured: in the dataset's features (${m.entities.length} gene entr${m.entities.length === 1 ? "y" : "ies"} across organisms and annotations).`,
        );
      for (const r of rows) {
        const s = stats.get(String(r.meyes_id));

        if (!s) continue;
        const list = geneStats.get(String(r.meyes_id)) ?? [];

        list.push({ ...s, symbol: g.gene });
        geneStats.set(String(r.meyes_id), list);
      }
    }
  }

  const collapsed = new Map<string, Row[]>();

  if (f.collapse_rereleases !== false) {
    const groups = new Map<string, Row[]>();

    for (const r of rows)
      if (r.same_section_relation === "rerelease" && r.same_section_group) {
        const g = String(r.same_section_group);

        groups.set(g, [...(groups.get(g) ?? []), r]);
      }
    const drop = new Set<Row>();

    for (const members of groups.values()) {
      if (members.length < 2) continue;
      members.sort(compareRelease);
      const newest = members[members.length - 1];
      const older = members.slice(0, -1);

      older.forEach((m) => drop.add(m));
      collapsed.set(String(newest.meyes_id), older);
    }
    if (drop.size) {
      rows = rows.filter((r) => !drop.has(r));
      notes.push(
        `${drop.size} re-release${drop.size === 1 ? "" : "s"} of the same data folded into the newest (listed under "also published as"; collapse_rereleases=false shows them all).`,
      );
    }
  }

  const sorted = [...rows];

  if (sort === "cells")
    sorted.sort((a, b) => (num(b.n_obs) ?? -1) - (num(a.n_obs) ?? -1));
  else if (sort === "title")
    sorted.sort((a, b) => lc(a.title).localeCompare(lc(b.title)));
  else
    sorted.sort((a, b) =>
      lc(b.publication_date ?? b.published_at).localeCompare(
        lc(a.publication_date ?? a.published_at),
      ),
    );

  return {
    release: index.release,
    total: sorted.length,
    notes,
    breakdown: breakdown(sorted),
    rows: sorted,
    geneStats,
    collapsed,
  };
}

/** A row as a tool returns it. */
export function hit(
  r: Row,
  res?: Pick<SearchResult, "geneStats" | "collapsed">,
): Hit {
  const id = String(r.meyes_id);
  const h: Hit = {
    id,
    title: str(r.title),
    organism: str(r.organism),
    tissue: str(r.tissue),
    tissue_organ: str(r.tissue_organ),
    disease: str(r.disease),
    technology: str(r.technology),
    cells: num(r.n_obs),
    obs_unit: str(r.obs_unit),
    genes_measured: num(r.n_features),
    timepoint: str(r.timepoint),
    age_days: lc(r.timepoint_kind) === "age" ? num(r.timepoint_days) : null,
    sex: str(r.meyes_sex),
    genotype: str(r.genotype),
    condition: str(r.condition),
    treatment: str(r.treatment),
    source: str(r.meyes_source),
    metadata_from: str(r.meyes_metadata_from),
    format: str(r.meyes_format),
    kind: str(r.meyes_kind),
    study_id: str(r.study_id),
    doi: str(r.doi),
    bil_code: str(r.meyes_bil_code),
    page_url: str(r.meyes_page_url),
    viewer_url: str(r.meyes_viewer_url),
  };
  const older = res?.collapsed.get(id);

  if (older?.length)
    h.also_published_as = older.map((o) => ({
      id: String(o.meyes_id),
      release: str(o.same_section_release),
    }));
  const g = res?.geneStats.get(id);

  if (g?.length)
    h.gene_stats = g.map((s) => ({
      gene: s.symbol,
      fraction_detected: s.fraction_obs_detected,
      total_counts: s.total_counts,
    }));

  return h;
}

/** Columns `meyes_index_values` lists. */
export const VALUE_COLUMNS = [
  "organism",
  "tissue",
  "tissue_organ",
  "tissue_system",
  "disease",
  "disease_class",
  "disease_group",
  "technology",
  "vendor",
  "modality",
  "obs_unit",
  "meyes_source",
  "meyes_kind",
  "meyes_metadata_from",
  "meyes_format",
  "meyes_sex",
  "development_stage",
  "timepoint",
  "genotype",
  "condition",
  "treatment",
  "license",
  "study_id",
  "meyes_investigator",
  "meyes_institution",
] as const;

/** Distinct values of a column among `rows`, most common first. */
export function values(
  rows: Row[],
  column: (typeof VALUE_COLUMNS)[number],
  containing?: string,
): { value: string; count: number }[] {
  const m = new Map<string, number>();
  const q = containing ? lc(containing) : "";

  for (const r of rows) {
    const v = str(r[column]);

    if (v === null || (q && !lc(v).includes(q))) continue;
    m.set(v, (m.get(v) ?? 0) + 1);
  }

  return [...m.entries()]
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
    .map(([value, count]) => ({ value, count }));
}
