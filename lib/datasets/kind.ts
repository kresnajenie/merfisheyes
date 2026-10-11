/**
 * The three kinds of dataset, and where each one opens.
 *
 * `Dataset.datasetType` is free text: older rows hold the raw format ("h5ad",
 * "merscope", ...) or nothing, and all of those are single cell.
 */
export type DatasetKind = "single_cell" | "single_molecule" | "image";

export function datasetKind(datasetType?: string | null): DatasetKind {
  return datasetType === "single_molecule" || datasetType === "image"
    ? datasetType
    : "single_cell";
}

const VIEWER: Record<DatasetKind, string> = {
  single_cell: "/viewer",
  single_molecule: "/sm-viewer",
  image: "/image-viewer",
};

/** `/viewer`, `/sm-viewer` or `/image-viewer`. */
export function viewerBase(datasetType?: string | null): string {
  return VIEWER[datasetKind(datasetType)];
}

/**
 * The page a dataset opens at: by URL when it was registered from a public
 * S3 link, by id when it lives in the app's own bucket.
 */
export function viewerHref(dataset: {
  id?: string | null;
  datasetType?: string | null;
  s3BaseUrl?: string | null;
}): string | null {
  const base = viewerBase(dataset.datasetType);

  if (dataset.s3BaseUrl) {
    const url = dataset.s3BaseUrl.replace(/\/+$/, "");

    return `${base}/from-s3?url=${encodeURIComponent(url)}`;
  }

  return dataset.id ? `${base}/${dataset.id}` : null;
}

export const KIND_LABEL: Record<DatasetKind, string> = {
  single_cell: "Single cell",
  single_molecule: "Single molecule",
  image: "Image",
};

/** Two-letter badge used in dense lists. */
export const KIND_BADGE: Record<DatasetKind, string> = {
  single_cell: "SC",
  single_molecule: "SM",
  image: "IMG",
};
