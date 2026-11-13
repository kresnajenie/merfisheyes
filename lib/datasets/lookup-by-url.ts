import { isLocalUrl } from "@/lib/utils/local-url";

/**
 * The site's record of a dataset opened by URL (GET /api/datasets/by-url):
 * its id, owner, viewer defaults — or null when it is unregistered, the
 * lookup fails, or the URL is on this computer (`isLocalUrl`), whose address
 * is never sent to the site. The one way the viewers ask.
 */
export async function lookupByUrl(
  url: string,
): Promise<Record<string, any> | null> {
  if (isLocalUrl(url)) return null;
  try {
    const r = await fetch(
      `/api/datasets/by-url?url=${encodeURIComponent(url)}`,
    );

    return r.ok ? await r.json() : null;
  } catch {
    return null;
  }
}
