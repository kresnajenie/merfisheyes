import { createHash } from "crypto";
import { gunzipSync } from "zlib";

import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { normalizeS3Url } from "@/lib/utils/viewer-config";

const MANIFEST_MAX_BYTES = 5 * 1024 * 1024; // 5 MB cap on the manifest fetch

/** SSRF guard: only https, and reject internal / private / metadata hosts. */
function isSafePublicUrl(u: URL): boolean {
  if (u.protocol !== "https:") return false;
  const host = u.hostname.toLowerCase();

  if (host === "localhost" || host.endsWith(".local")) return false;
  // Literal private / link-local / loopback IPv4 ranges + IPv6 loopback.
  if (
    /^127\./.test(host) ||
    /^10\./.test(host) ||
    /^192\.168\./.test(host) ||
    /^169\.254\./.test(host) ||
    /^172\.(1[6-9]|2\d|3[01])\./.test(host) ||
    host === "0.0.0.0" ||
    host === "::1"
  ) {
    return false;
  }

  return true;
}

/** Fetch + parse the dataset manifest (manifest.json, else gzipped). Never throws. */
async function fetchManifest(baseUrl: string): Promise<any | null> {
  for (const suffix of ["/manifest.json", "/manifest.json.gz"]) {
    try {
      const target = new URL(baseUrl + suffix);

      if (!isSafePublicUrl(target)) return null;

      const controller = new AbortController();
      const timeout = setTimeout(() => controller.abort(), 10_000);
      const res = await fetch(target.toString(), {
        signal: controller.signal,
        redirect: "error", // don't follow redirects (SSRF)
      });

      clearTimeout(timeout);
      if (!res.ok) continue;

      const buf = Buffer.from(await res.arrayBuffer());

      if (buf.byteLength > MANIFEST_MAX_BYTES) return null;

      const text = suffix.endsWith(".gz")
        ? gunzipSync(buf).toString("utf-8")
        : buf.toString("utf-8");

      return JSON.parse(text);
    } catch {
      // try the next suffix
    }
  }

  return null;
}

/**
 * Register (or claim) a dataset that is loaded from a raw S3 URL so it can
 * carry ownership, view counts, and owner-saved viewer defaults. Dedups on the
 * normalized S3 base URL — re-registering the same URL never creates a
 * duplicate.
 */
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: { url?: string; asAdmin?: boolean };

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  if (!body.url) {
    return NextResponse.json({ error: "Missing url" }, { status: 400 });
  }

  const s3BaseUrl = normalizeS3Url(body.url);

  try {
    new URL(s3BaseUrl);
  } catch {
    return NextResponse.json({ error: "Invalid url" }, { status: 400 });
  }

  const userId = session.user.id;
  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
  // Admins may claim on behalf of "admin" (shared). Everyone else — and admins
  // who opt out — claim personally.
  const adminClaim = !!body.asAdmin && isAdmin;

  // Already registered? Claim it if unclaimed; otherwise return as-is.
  const existing = await prisma.dataset.findUnique({ where: { s3BaseUrl } });

  if (existing) {
    const unclaimed = existing.ownerId === null && !existing.adminOwned;

    if (unclaimed) {
      const claimed = await prisma.dataset.update({
        where: { id: existing.id },
        data: adminClaim ? { adminOwned: true } : { ownerId: userId },
      });

      return NextResponse.json({ dataset: claimed, claimed: true });
    }

    // Owned: "mine" if I'm the personal owner, or it's admin-owned and I'm an admin.
    const mine =
      existing.ownerId === userId || (existing.adminOwned && isAdmin);

    return NextResponse.json({
      dataset: existing,
      claimed: false,
      ownedByOther: !mine,
    });
  }

  // New registration — read counts/title from the manifest (best-effort).
  const manifest = await fetchManifest(s3BaseUrl);
  const stats = manifest?.statistics ?? {};
  // Single-molecule manifests report total_molecules; single-cell report
  // total_cells. Normalize to the pipeline's datasetType so the viewer routes
  // to /sm-viewer vs /viewer correctly.
  const isSingleMolecule = stats.total_molecules != null;
  const datasetType = isSingleMolecule
    ? "single_molecule"
    : (manifest?.type ?? null);
  const id = `s3_${nanoid(10)}`;
  const fingerprint = createHash("sha256").update(s3BaseUrl).digest("hex");

  const dataset = await prisma.dataset.create({
    data: {
      id,
      fingerprint,
      title: manifest?.name ?? "Untitled dataset",
      // total_cells/total_genes for single-cell; total_molecules/unique_genes
      // for single-molecule manifests.
      numCells: Number(stats.total_cells ?? stats.total_molecules) || 0,
      numGenes: Number(stats.total_genes ?? stats.unique_genes) || 0,
      datasetType,
      status: "COMPLETE",
      ownerId: adminClaim ? null : userId,
      adminOwned: adminClaim,
      ingestSource: "s3_registered",
      s3BaseUrl,
      manifestUrl: s3BaseUrl,
    },
  });

  return NextResponse.json({ dataset, claimed: true });
}
