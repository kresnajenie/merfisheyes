export interface SubmissionStatus {
  catalogId: string;
  reviewStatus: "PENDING" | "APPROVED" | "REJECTED";
  isPublished: boolean;
  reviewNote: string | null;
}

/** Owner-editable catalog-style metadata shared by datasets and projects. */
export interface EditableMetadata {
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

/** Flexible display-metadata bag (investigator, authors, citation, age, …). */
export type MetadataBag = Record<string, string>;

export interface DatasetRow extends EditableMetadata {
  id: string;
  title: string | null;
  status: string;
  datasetType: string | null;
  numCells: number;
  numGenes: number;
  viewCount: number;
  ingestSource: string | null;
  adminOwned: boolean;
  thumbnailUrl: string | null;
  metadata: MetadataBag | null;
  createdAt: string;
  viewerUrl: string | null;
  projectNames: string[];
  submission: SubmissionStatus | null;
}

export interface ProjectDatasetSummary {
  id: string;
  title: string | null;
  datasetType: string | null;
  thumbnailUrl: string | null;
  status: string;
}

export interface ProjectRow extends EditableMetadata {
  id: string;
  title: string;
  thumbnailUrl: string | null;
  metadata: MetadataBag | null;
  datasetCount: number;
  datasets: ProjectDatasetSummary[];
  createdAt: string;
  updatedAt: string;
  submission: SubmissionStatus | null;
}

export const EMPTY_METADATA: EditableMetadata = {
  description: null,
  species: null,
  disease: null,
  tissue: null,
  platform: null,
  institute: null,
  tags: [],
  externalLink: null,
  publicationLink: null,
};
