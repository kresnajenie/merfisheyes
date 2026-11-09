/**
 * What a pasted link holds as an image, for an overlay: one rule for the
 * server (which writes a saved dataset's overlay) and the browser (which
 * shows what was found). Plain `fetch`, so it runs in both — the same shape
 * as `moleculeSource.ts` next to it.
 */
import { imageMultiscale } from "@/lib/spatialdata/ome";
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";

export { dataUrlOf } from "./moleculeSource";

/** An image a link offers, by the URL `mapping.json`'s `image` would name. */
export interface ImageTarget {
  /** The image group itself — what the viewer opens (openOmeZarrImage). */
  url: string;
  /**
   * Its element name, when it came from inside a store (`images/<name>`).
   * Absent for an image that is the link itself.
   */
  name?: string;
}

/**
 * The image elements of a store, from its consolidated metadata: a node at
 * `images/<name>` carrying multiscale metadata, in either layout
 * (lib/spatialdata/ome.ts).
 */
function storeImageNames(
  root: Record<string, any> | null | undefined,
): string[] {
  const nodes: Record<string, any> =
    root?.consolidated_metadata?.metadata ?? {};

  return Object.entries(nodes)
    .filter(
      ([key, node]) =>
        /^images\/[^/]+$/.test(key) && !!imageMultiscale(node?.attributes),
    )
    .map(([key]) => key.slice("images/".length))
    .sort();
}

/**
 * Every image the link offers, in the order the viewer would list them:
 *   - an OME-Zarr image → itself, one entry
 *   - a SpatialData store, or a folder holding one → its `images/` elements,
 *     whatever else the store holds. A store with cells and molecules too
 *     lends only its image, the way a pasted molecule link lends only its
 *     points.
 * Empty when there is no image, or the URL can't be read (private, blocked
 * by CORS, a typo).
 *
 * An element's URL is the group inside the store, so the viewer opens it as
 * a standalone image and places it by its own transform (the one into
 * `global`). A store whose cells sit in some other coordinate system would
 * need the cells to place it, which an overlay by link has no way to know;
 * such an image lands where its own transform says.
 */
export async function imageUrlTargets(url: string): Promise<ImageTarget[]> {
  const base = url.trim().replace(/\/+$/, "");
  const target = await resolveZarrUrl(base);

  if (!target) return [];
  if (target.kind === "image") return [{ url: target.url }];

  return storeImageNames(target.root).map((name) => ({
    url: `${target.url}/images/${name}`,
    name,
  }));
}
