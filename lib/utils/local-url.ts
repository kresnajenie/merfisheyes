/**
 * Whether a URL points at this computer (127.0.0.1, ::1, localhost): data a
 * person opens from their own machine — `merfisheyes-mcp view`, which serves
 * a folder on 127.0.0.1 for the viewer to read. Such an address is never sent
 * to the site: no registration lookup, no claim banner, no analytics.
 */
export function isLocalUrl(url: string | null | undefined): boolean {
  if (!url) return false;
  let host: string;

  try {
    host = new URL(url).hostname.toLowerCase();
  } catch {
    return false;
  }

  return (
    host === "localhost" ||
    host.endsWith(".localhost") ||
    host === "[::1]" ||
    /^127(\.\d{1,3}){3}$/.test(host)
  );
}

/** A page address without the parameters that point at this computer
 * (`?url=http://127.0.0.1:…`, `splitS3Url=…`): a local folder's name is never
 * sent to analytics. */
export function withoutLocalUrls(href: string): string {
  const u = new URL(href);

  for (const [key, value] of [...u.searchParams]) {
    let v = value;

    try {
      v = decodeURIComponent(value);
    } catch {
      // keep it as it is
    }
    if (isLocalUrl(value) || isLocalUrl(v)) u.searchParams.delete(key);
  }

  return u.toString();
}

/**
 * Whether the viewer keeps a dataset's address to itself — no registration
 * lookup, no claim banner. In the standalone viewer (no site, no account)
 * always. Otherwise for data on this computer, when the page came from the
 * internet: www.merfisheyes.com reading someone's own folder. A page served
 * from this computer too (a dev server, the end-to-end tests' stand-in
 * bucket) reads local data as it would a bucket's.
 */
export function keepPrivate(
  url: string | null | undefined,
  page: string | null = typeof window === "undefined"
    ? null
    : window.location.href,
): boolean {
  if (process.env.NEXT_PUBLIC_MERFISHEYES_STANDALONE === "1") return true;

  return isLocalUrl(url) && !isLocalUrl(page);
}
