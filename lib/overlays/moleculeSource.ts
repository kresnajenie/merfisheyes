/**
 * What a pasted link holds as molecules, for an overlay: one rule for the
 * server (which writes a saved dataset's overlay) and the browser (which
 * checks one over a local preview). Plain `fetch`, so it runs in both.
 */
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";

/** How the viewer opens the molecules (DatasetLinkConfig.linkFormat). */
export type MoleculeLinkFormat = "store" | "molecules";

/**
 * The public data URL a pasted link names: the link itself, or, for a
 * merfisheyes viewer link opened from S3 (`…/from-s3?url=…`), the URL it
 * opened. Null for anything that is not an http(s) URL.
 */
export function dataUrlOf(input: string): string | null {
  let url: URL;

  try {
    url = new URL(input.trim());
  } catch {
    return null;
  }
  if (url.protocol !== "http:" && url.protocol !== "https:") return null;
  const inner = url.searchParams.get("url");

  if (inner && /\/from-s3\/?$/.test(url.pathname)) return dataUrlOf(inner);

  return url.href.replace(/\/+$/, "");
}

/**
 * Whether a store's root names a points element in the format the viewer
 * draws (`points:version`), from its consolidated metadata.
 */
export function rootHasPoints(
  root: Record<string, any> | null | undefined,
): boolean {
  const nodes: Record<string, any> =
    root?.consolidated_metadata?.metadata ?? {};

  return Object.entries(nodes).some(
    ([key, node]) =>
      /^points\/[^/]+$/.test(key) &&
      node?.attributes?.["points:version"] != null,
  );
}

/**
 * How the viewer would open the molecules at a public URL: a SpatialData
 * store with a points element, or a folder holding one (`"store"`; only its
 * points are drawn, whatever cells it has too), or a merfisheyes molecule
 * folder (`"molecules"`, its manifest.json.gz). Null when there are none it
 * can draw, or the URL can't be read (private, not allowed by CORS, a typo).
 */
export async function moleculeUrlFormat(
  url: string,
): Promise<MoleculeLinkFormat | null> {
  const base = url.trim().replace(/\/+$/, "");
  const target = await resolveZarrUrl(base);

  if (target) {
    return target.kind === "spatialdata" && rootHasPoints(target.root)
      ? "store"
      : null;
  }
  try {
    const res = await fetch(`${base}/manifest.json.gz`);

    return res.ok ? "molecules" : null;
  } catch {
    return null;
  }
}
