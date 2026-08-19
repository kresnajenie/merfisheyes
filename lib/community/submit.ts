/**
 * Community submission helper. Turns an owned Dataset (or a Project of owned
 * datasets) into a single community CatalogDataset that goes through admin
 * review. Re-submitting updates the same row (tracked via sourceDatasetId /
 * sourceProjectId) rather than creating duplicates.
 */
import { prisma } from "@/lib/prisma";

// The Dataset fields this helper reads. Kept narrow so callers can pass a
// `select`-ed subset.
export interface SubmitDataset {
  id: string;
  title: string | null;
  datasetType: string | null;
  ingestSource: string | null;
  s3BaseUrl: string | null;
  thumbnailUrl: string | null;
  description: string | null;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  tags: string[];
  externalLink: string | null;
  publicationLink: string | null;
  metadata: unknown;
  numCells: number | null;
  numGenes: number | null;
  /** Gene symbols (filled at processing time) — carried into the catalog row
   * so the submission is gene-searchable on Explore. */
  genes?: string[];
}

export interface SubmitProject {
  id: string;
  title: string;
  thumbnailUrl: string | null;
  description: string | null;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  tags: string[];
  externalLink: string | null;
  publicationLink: string | null;
  metadata: unknown;
  // Ordered member datasets (by ProjectDataset.sortOrder).
  datasets: SubmitDataset[];
}

interface EntryInput {
  label: string;
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
  sortOrder: number;
}

/**
 * Build a CatalogDatasetEntry input from a Dataset, mirroring how the viewer
 * resolves a dataset's URL: S3-registered rows load from their base URL, all
 * others load by id.
 */
export function entryForDataset(d: SubmitDataset): {
  label: string;
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
  sortOrder: number;
} {
  const isS3Registered =
    d.ingestSource === "s3_registered" && Boolean(d.s3BaseUrl);

  return {
    label: d.title || "Dataset",
    datasetType: d.datasetType || "single_cell",
    s3BaseUrl: isS3Registered ? d.s3BaseUrl!.replace(/\/+$/, "") : null,
    datasetId: isS3Registered ? null : d.id,
    sortOrder: 0,
  };
}

type Source =
  | { kind: "dataset"; dataset: SubmitDataset }
  | { kind: "project"; project: SubmitProject };

/**
 * Create or update the community CatalogDataset for an owned dataset/project.
 * Idempotent per source: re-submitting replaces the row's metadata + entries.
 */
export async function upsertCommunitySubmission({
  userId,
  source,
}: {
  userId: string;
  source: Source;
}) {
  const sourceDatasetId = source.kind === "dataset" ? source.dataset.id : null;
  const sourceProjectId = source.kind === "project" ? source.project.id : null;

  const src = source.kind === "dataset" ? source.dataset : source.project;

  // Column fields on CatalogDataset (distinct from the `metadata` Json bag).
  const fields = {
    title: src.title || "Dataset",
    description: src.description ?? null,
    species: src.species ?? null,
    disease: src.disease ?? null,
    tissue: src.tissue ?? null,
    platform: src.platform ?? null,
    institute: src.institute ?? null,
    tags: src.tags ?? [],
    externalLink: src.externalLink ?? null,
    publicationLink: src.publicationLink ?? null,
    thumbnailUrl: src.thumbnailUrl ?? null,
    // Flexible display-metadata bag (investigator, authors, citation, …),
    // rendered by the Explore detail page.
    metadata: (src.metadata ?? {}) as object,
    numCells: source.kind === "dataset" ? source.dataset.numCells : null,
    numGenes: source.kind === "dataset" ? source.dataset.numGenes : null,
    // Union of member genes for a project; makes the catalog row
    // gene-searchable on Explore.
    genes:
      source.kind === "dataset"
        ? (source.dataset.genes ?? [])
        : Array.from(
            new Set(source.project.datasets.flatMap((d) => d.genes ?? [])),
          ),
  };

  const entries: EntryInput[] =
    source.kind === "dataset"
      ? [entryForDataset(source.dataset)]
      : source.project.datasets.map((d, i) => ({
          ...entryForDataset(d),
          sortOrder: i,
        }));

  const existing = await prisma.catalogDataset.findFirst({
    where: {
      isCommunity: true,
      ...(sourceDatasetId
        ? { sourceDatasetId }
        : { sourceProjectId: sourceProjectId! }),
    },
    select: { id: true, reviewStatus: true },
  });

  const includeEntries = {
    entries: { orderBy: { sortOrder: "asc" as const } },
  };

  if (!existing) {
    return prisma.catalogDataset.create({
      data: {
        ...fields,
        isCommunity: true,
        reviewStatus: "PENDING",
        isPublished: false,
        submittedById: userId,
        createdBy: userId,
        sourceDatasetId,
        sourceProjectId,
        entries: { create: entries },
      },
      include: includeEntries,
    });
  }

  // Re-submit: replace entries and refresh metadata.
  await prisma.catalogDatasetEntry.deleteMany({
    where: { catalogId: existing.id },
  });

  // An owner editing their already-live listing keeps it live; otherwise the
  // edit re-enters the review queue.
  const keepApproved = existing.reviewStatus === "APPROVED";

  return prisma.catalogDataset.update({
    where: { id: existing.id },
    data: {
      ...fields,
      submittedById: userId,
      ...(keepApproved
        ? {}
        : {
            reviewStatus: "PENDING",
            isPublished: false,
            reviewNote: null,
          }),
      entries: { create: entries },
    },
    include: includeEntries,
  });
}
