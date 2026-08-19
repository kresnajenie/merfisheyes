export interface CatalogDatasetEntry {
  id: string;
  label: string;
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
  thumbnailUrl: string | null;
  sortOrder: number;
}

export interface CatalogDatasetItem {
  id: string;
  title: string;
  description: string | null;
  species: string | null;
  disease: string | null;
  institute: string | null;
  tissue: string | null;
  platform: string | null;
  tags: string[];
  /**
   * Full gene list. Only present on single-dataset payloads (the detail
   * page); list/card payloads omit it for size and carry `matchedGenes`
   * instead.
   */
  genes?: string[];
  /** Genes that matched the active gene search, capped, for card chips. */
  matchedGenes?: string[];
  thumbnailUrl: string | null;
  bilCode: string | null;
  metadata: Record<string, unknown> | null;
  externalLink: string | null;
  publicationLink: string | null;
  entries: CatalogDatasetEntry[];
  isPublished: boolean;
  isFeatured: boolean;
  isBil: boolean;
  isInternal: boolean;
  isCommunity: boolean;
  /** Set when a community row was submitted from a Project — such rows open
   *  the detail page even with a single entry. Optional: not every payload
   *  (e.g. account-card adapters) carries it. */
  sourceProjectId?: string | null;
  sortOrder: number;
  numCells: number | null;
  numGenes: number | null;
}

export interface ExploreFilters {
  species: string[];
  tissues: string[];
  platforms: string[];
}

export interface ExploreApiResponse {
  items: CatalogDatasetItem[];
  total: number;
  page: number;
  limit: number;
  /** Present only when the request asked for `include=meta`. */
  featured?: CatalogDatasetItem[];
  filters?: ExploreFilters;
}
