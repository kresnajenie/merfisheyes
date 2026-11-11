import { createHash, createHmac, randomBytes, timingSafeEqual } from "crypto";

/** An opaque secret: 32 random bytes, base64url, with a readable prefix so a
 * leaked token is recognisable in logs and secret scanners. */
export function newSecret(prefix: "mea" | "mer" | "mec"): string {
  return `${prefix}_${randomBytes(32).toString("base64url")}`;
}

/** What is stored: codes and tokens are never kept in clear. */
export function sha256(value: string): string {
  return createHash("sha256").update(value).digest("hex");
}

/** PKCE S256 (RFC 7636 §4.6): BASE64URL(SHA256(verifier)) == challenge. */
export function pkceMatches(verifier: string, challenge: string): boolean {
  if (!/^[A-Za-z0-9\-._~]{43,128}$/.test(verifier)) return false;
  const computed = createHash("sha256").update(verifier).digest("base64url");

  return safeEqual(computed, challenge);
}

export function isValidChallenge(challenge: string): boolean {
  // a SHA-256 digest, base64url without padding
  return /^[A-Za-z0-9_-]{43}$/.test(challenge);
}

export function safeEqual(a: string, b: string): boolean {
  const x = Buffer.from(a);
  const y = Buffer.from(b);

  return x.length === y.length && timingSafeEqual(x, y);
}

function key(): string {
  const k = process.env.AUTH_SECRET || process.env.NEXTAUTH_SECRET;

  if (!k) throw new Error("AUTH_SECRET is not set");

  return k;
}

/** Sign a JSON value (the consent form's request), so the decision POST can
 * trust what the authorize page validated. */
export function sign(value: unknown): string {
  const body = Buffer.from(JSON.stringify(value)).toString("base64url");
  const mac = createHmac("sha256", key())
    .update(`oauth-consent.${body}`)
    .digest("base64url");

  return `${body}.${mac}`;
}

export function verifySigned<T>(signed: string): T | null {
  const [body, mac] = signed.split(".");

  if (!body || !mac) return null;
  const expected = createHmac("sha256", key())
    .update(`oauth-consent.${body}`)
    .digest("base64url");

  if (!safeEqual(mac, expected)) return null;
  try {
    return JSON.parse(Buffer.from(body, "base64url").toString()) as T;
  } catch {
    return null;
  }
}
