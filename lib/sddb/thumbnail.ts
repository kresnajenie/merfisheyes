/**
 * A SpatialData-DB card's thumbnail. SpatialData-DB serves one image per
 * dataset at a URL made from the dataset's id alone: the first 16 characters
 * of the catalog `uid` (the last four are the version, and the image stays at
 * the same URL across versions). Not every dataset has one; a missing image
 * is a 404, which the cards hide.
 */
export function sddbThumbnailUrl(uid: string): string {
  return `https://sdata-db.com/t/${uid.slice(0, 16)}.webp`;
}

/** The image a card shows: its own thumbnail, else SpatialData-DB's. */
export function cardThumbnailUrl(dataset: {
  thumbnailUrl: string | null;
  source?: string;
  sddbUid?: string | null;
}): string | null {
  if (dataset.thumbnailUrl) return dataset.thumbnailUrl;

  return dataset.source === "sddb" && dataset.sddbUid
    ? sddbThumbnailUrl(dataset.sddbUid)
    : null;
}
