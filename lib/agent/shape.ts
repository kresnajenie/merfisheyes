import type {
  CatalogDataset,
  CatalogDatasetEntry,
  Dataset,
} from "@prisma/client";

import { datasetKind, viewerHref } from "@/lib/datasets/kind";
import { catalogPageHref, entryHref } from "@/lib/explore/links";
import { issuer } from "@/lib/oauth/config";
import { readSidecars, summarise } from "@/lib/agent/sidecars";

/**
 * What the agent API returns for a catalog row, a dataset, a project — one
 * shape per thing, absolute URLs, no gene arrays (those are paged on their
 * own route).
 */

const abs = (path: string | null) => (path ? `${issuer()}${path}` : null);

type EntryLike = Pick<
  CatalogDatasetEntry,
  "label" | "datasetType" | "s3BaseUrl" | "datasetId"
>;

export function entryOut(e: EntryLike) {
  return {
    label: e.label,
    kind: datasetKind(e.datasetType),
    datasetType: e.datasetType,
    datasetId: e.datasetId,
    url: e.s3BaseUrl,
    viewerUrl: abs(entryHref(e)),
  };
}

type CardLike = Pick<
  CatalogDataset,
  | "id"
  | "title"
  | "species"
  | "tissue"
  | "platform"
  | "bilCode"
  | "numCells"
  | "numGenes"
  | "source"
  | "license"
> & {
  entries: EntryLike[];
  matchedGenes?: string[];
  description?: string | null;
};

export function cardOut(c: CardLike) {
  return {
    id: c.id,
    title: c.title,
    species: c.species,
    tissue: c.tissue,
    platform: c.platform,
    bilCode: c.bilCode,
    cells: c.numCells,
    genes: c.numGenes,
    source: c.source,
    license: c.license,
    ...(c.matchedGenes ? { matchedGenes: c.matchedGenes } : {}),
    pageUrl: abs(catalogPageHref(c)),
    entries: c.entries.map(entryOut),
  };
}

/** A store's sidecar files, once per folder, for the entries that are stores. */
async function sidecarsFor(urls: (string | null)[]) {
  const unique = [...new Set(urls.filter((u): u is string => !!u))].slice(0, 5);
  const read = await Promise.all(unique.map((u) => readSidecars(u)));

  return read
    .filter(
      (s): s is NonNullable<typeof s> =>
        !!s && !!(s.catalog || s.stats || s.qc),
    )
    .map((s) => ({
      folder: s.base,
      summary: summarise(s),
      files: {
        catalog: s.catalog,
        stats: s.stats ? { ...s.stats, packages: undefined } : null,
        qc: s.qc,
      },
    }));
}

export async function catalogDetailOut(
  row: CatalogDataset & { entries: CatalogDatasetEntry[] },
) {
  return {
    type: "catalog" as const,
    ...cardOut(row),
    description: row.description,
    disease: row.disease,
    institute: row.institute,
    tags: row.tags,
    geneCount: row.genes.length,
    metadata: row.metadata,
    externalLink: row.externalLink,
    publicationLink: row.publicationLink,
    licenseNoncommercial: row.licenseNoncommercial,
    sidecars: await sidecarsFor(row.entries.map((e) => e.s3BaseUrl)),
  };
}

type DatasetLike = Pick<
  Dataset,
  | "id"
  | "title"
  | "datasetType"
  | "status"
  | "numCells"
  | "numGenes"
  | "species"
  | "tissue"
  | "platform"
  | "s3BaseUrl"
  | "createdAt"
  | "completedAt"
>;

export function datasetCardOut(d: DatasetLike) {
  return {
    id: d.id,
    title: d.title,
    kind: datasetKind(d.datasetType),
    datasetType: d.datasetType,
    status: d.status,
    cells: d.numCells,
    genes: d.numGenes,
    species: d.species,
    tissue: d.tissue,
    platform: d.platform,
    url: d.s3BaseUrl,
    viewerUrl: d.status === "COMPLETE" ? abs(viewerHref(d)) : null,
    createdAt: d.createdAt,
    completedAt: d.completedAt,
  };
}

export async function datasetDetailOut(d: Dataset, owned: boolean) {
  return {
    type: "dataset" as const,
    ...datasetCardOut(d),
    owned,
    description: d.description,
    disease: d.disease,
    institute: d.institute,
    tags: d.tags,
    geneCount: d.genes.length,
    metadata: d.metadata,
    externalLink: d.externalLink,
    publicationLink: d.publicationLink,
    sidecars: await sidecarsFor([d.s3BaseUrl]),
  };
}

/** Genes matching a substring, paged. */
export function genesPage(
  genes: string[],
  q: string,
  p: { limit: number; offset: number },
) {
  const needle = q.trim().toLowerCase();
  const hits = needle
    ? genes.filter((g) => g.toLowerCase().includes(needle))
    : genes;

  return { hits: hits.slice(p.offset, p.offset + p.limit), total: hits.length };
}
