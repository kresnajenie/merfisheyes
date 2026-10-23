/**
 * The three kinds of dataset, and where each one opens.
 *
 * `Dataset.datasetType` is free text: older rows hold the raw format ("h5ad",
 * "merscope", ...) or nothing, and all of those are single cell. So is
 * "spatialdata", which only differs in the viewer it opens in.
 */
export type DatasetKind = "single_cell" | "single_molecule" | "image";

export function datasetKind(datasetType?: string | null): DatasetKind {
  return datasetType === "single_molecule" || datasetType === "image"
    ? datasetType
    : "single_cell";
}

/**
 * A SpatialData store saved as it is: cells like any other single cell
 * dataset (it takes overlays, is badged the same), but read by its own viewer.
 */
const SPATIALDATA = "spatialdata";

/**
 * The type to put on a card entry: the kind, except that what routing needs
 * to know is kept.
 */
export function entryType(datasetType?: string | null): string {
  return datasetType === SPATIALDATA ? SPATIALDATA : datasetKind(datasetType);
}

const VIEWER: Record<DatasetKind, string> = {
  single_cell: "/viewer",
  single_molecule: "/sm-viewer",
  image: "/image-viewer",
};

/** `/viewer`, `/sm-viewer` or `/image-viewer`. */
export function viewerBase(datasetType?: string | null): string {
  return datasetType === SPATIALDATA
    ? "/spatialdata-viewer"
    : VIEWER[datasetKind(datasetType)];
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

/**
 * A dataset in the app's own bucket (`datasets/{id}/`): `ds_…`, or `sm_…`
 * for one uploaded from the single molecule card. The id's prefix is where
 * the upload began, not what the dataset is: a molecule table processed by
 * scripts/prepare.py keeps its `sm_` id and is a SpatialData store.
 */
export function isAppDatasetId(id: string | null | undefined): id is string {
  return !!id && (id.startsWith("ds_") || id.startsWith("sm_"));
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
