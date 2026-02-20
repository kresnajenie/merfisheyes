export interface CatalogDatasetEntry {
  id: string;
  label: string;
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
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
  thumbnailUrl: string | null;
  externalLink: string | null;
  entries: CatalogDatasetEntry[];
  isPublished: boolean;
  isFeatured: boolean;
  isBil: boolean;
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
  featured: CatalogDatasetItem[];
  bil: CatalogDatasetItem[];
  filters: ExploreFilters;
}
