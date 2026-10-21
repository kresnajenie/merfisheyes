import { viewerHref } from "@/lib/datasets/kind";

/**
 * Where the Explore surfaces send a visitor. One place for both links so the
 * grid card, the featured carousel, the detail page and the info modal can't
 * disagree.
 */

/**
 * Every catalog row has a standalone page, and that is where a card goes:
 * `/explore/bil/[code]` for a BIL row, `/explore/[id]` for everything else.
 * The page lists the row's entries (cells, molecules, images) and opens the
 * viewer from there.
 */
export function catalogPageHref(dataset: {
  id: string;
  bilCode?: string | null;
}): string {
  return dataset.bilCode
    ? `/explore/bil/${dataset.bilCode}`
    : `/explore/${dataset.id}`;
}

/**
 * Where an entry opens: the viewer for its kind (`/viewer`, `/sm-viewer`,
 * `/image-viewer`, `/spatialdata-viewer`), by public URL when it has one,
 * else by dataset id. Null when the entry has neither.
 */
export function entryHref(entry: {
  datasetType: string;
  s3BaseUrl: string | null;
  datasetId: string | null;
}): string | null {
  return viewerHref({
    id: entry.datasetId,
    datasetType: entry.datasetType,
    s3BaseUrl: entry.s3BaseUrl,
  });
}
