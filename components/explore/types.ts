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
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
  externalLink: string | null;
  isPublished: boolean;
  isFeatured: boolean;
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
  filters: ExploreFilters;
}
