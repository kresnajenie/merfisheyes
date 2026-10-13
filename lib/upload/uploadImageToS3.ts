import {
  uploadZarrFolderToS3,
  type ZarrFolderUploadOptions,
  type ZarrUploadFile,
} from "./uploadZarrFolderToS3";

/**
 * Upload an OME-Zarr image folder to S3 as an image dataset, as it is: an
 * image needs no processing, the viewer reads the zarr directly. See
 * uploadZarrFolderToS3 for how the many small files go up.
 */

/** Path inside the image folder, e.g. "zarr.json" or "0/c/0/3/7". */
export type ImageUploadFile = ZarrUploadFile;

export type ImageUploadOptions = Omit<
  ZarrFolderUploadOptions,
  "platform" | "lastKey" | "numCells" | "numGenes"
>;

/**
 * True when a root `zarr.json` describes an OME-Zarr multiscale image the
 * viewer can draw (Zarr v3, axes c/y/x) rather than some other zarr store.
 */
export function isOmeZarrImage(rootJson: unknown): boolean {
  const multiscale = (
    rootJson as {
      attributes?: { ome?: { multiscales?: { axes?: { name?: string }[] }[] } };
    } | null
  )?.attributes?.ome?.multiscales?.[0];

  return multiscale?.axes?.map((a) => a.name).join("") === "cyx";
}

export function uploadImageToS3(
  opts: ImageUploadOptions,
): Promise<{ datasetId: string }> {
  return uploadZarrFolderToS3({
    ...opts,
    platform: "image",
    lastKey: "zarr.json",
  });
}
