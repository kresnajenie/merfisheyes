import { lookup } from "dns/promises";
import { isIP } from "net";

/**
 * A server-side GET of a URL someone else chose — an agent client's metadata
 * document, a store's sidecar JSON — bounded so it cannot be turned against
 * us: HTTPS only, no redirects, a time limit, a size limit, and never to a
 * private, loopback or link-local address (server-side request forgery). A
 * name is resolved and every address it resolves to is checked; the check
 * and the connection are separate lookups, so this narrows, not closes, a
 * DNS-rebinding window — enough for public metadata, not for secrets.
 */

export class LimitedFetchError extends Error {}

export interface LimitedFetchOptions {
  maxBytes: number;
  timeoutMs: number;
  /** Further restriction on the host, e.g. S3 only. */
  allowHost?: (hostname: string) => boolean;
}

function privateV4(ip: string): boolean {
  const [a, b] = ip.split(".").map(Number);

  return (
    a === 0 ||
    a === 10 ||
    a === 127 ||
    (a === 100 && b >= 64 && b <= 127) || // carrier-grade NAT
    (a === 169 && b === 254) || // link-local, cloud metadata
    (a === 172 && b >= 16 && b <= 31) ||
    (a === 192 && b === 168) ||
    (a === 198 && (b === 18 || b === 19)) ||
    a >= 224
  );
}

export function isPrivateAddress(ip: string): boolean {
  if (isIP(ip) === 4) return privateV4(ip);
  const v6 = ip.toLowerCase();

  if (v6.startsWith("::ffff:")) return privateV4(v6.slice(7));

  return (
    v6 === "::" ||
    v6 === "::1" ||
    v6.startsWith("fc") ||
    v6.startsWith("fd") || // unique local
    v6.startsWith("fe8") ||
    v6.startsWith("fe9") ||
    v6.startsWith("fea") ||
    v6.startsWith("feb") || // link-local
    v6.startsWith("ff") // multicast
  );
}

export async function assertPublicHttps(
  url: URL,
  opts: Pick<LimitedFetchOptions, "allowHost">,
) {
  if (url.protocol !== "https:")
    throw new LimitedFetchError("only https URLs are fetched");
  if (url.username || url.password)
    throw new LimitedFetchError("credentials in the URL");
  const host = url.hostname.replace(/^\[|\]$/g, "");

  if (opts.allowHost && !opts.allowHost(host))
    throw new LimitedFetchError(`host ${host} is not allowed`);
  const addresses = isIP(host)
    ? [host]
    : (await lookup(host, { all: true, verbatim: true }).catch(() => [])).map(
        (a) => a.address,
      );

  if (addresses.length === 0)
    throw new LimitedFetchError(`cannot resolve ${host}`);
  if (addresses.some(isPrivateAddress))
    throw new LimitedFetchError(`${host} is not a public address`);
}

export interface LimitedResponse {
  status: number;
  body: string;
  headers: Headers;
}

/** GET a URL within the limits. A non-2xx answer is returned, not thrown. */
export async function limitedFetch(
  raw: string,
  opts: LimitedFetchOptions,
): Promise<LimitedResponse> {
  let url: URL;

  try {
    url = new URL(raw);
  } catch {
    throw new LimitedFetchError("not a URL");
  }
  await assertPublicHttps(url, opts);

  const res = await fetch(url, {
    redirect: "manual",
    signal: AbortSignal.timeout(opts.timeoutMs),
    headers: { accept: "application/json" },
    cache: "no-store",
  });

  if (res.status >= 300 && res.status < 400)
    throw new LimitedFetchError("redirects are not followed");
  const declared = Number(res.headers.get("content-length") ?? "0");

  if (declared > opts.maxBytes)
    throw new LimitedFetchError(`larger than ${opts.maxBytes} bytes`);

  const reader = res.body?.getReader();
  const chunks: Uint8Array[] = [];
  let size = 0;

  while (reader) {
    const { done, value } = await reader.read();

    if (done) break;
    size += value.byteLength;
    if (size > opts.maxBytes) {
      await reader.cancel();
      throw new LimitedFetchError(`larger than ${opts.maxBytes} bytes`);
    }
    chunks.push(value);
  }

  return {
    status: res.status,
    body: Buffer.concat(chunks).toString("utf8"),
    headers: res.headers,
  };
}

/** Seconds a response may be cached, from Cache-Control max-age, clamped. */
export function maxAgeSeconds(
  headers: Headers,
  fallback: number,
  min: number,
  max: number,
): number {
  const m = /max-age=(\d+)/i.exec(headers.get("cache-control") ?? "");
  const noStore = /no-store|no-cache/i.test(headers.get("cache-control") ?? "");
  const s = noStore ? min : m ? Number(m[1]) : fallback;

  return Math.min(max, Math.max(min, s));
}
