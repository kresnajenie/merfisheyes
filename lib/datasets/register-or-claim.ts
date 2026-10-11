import { createHash } from "crypto";

import { nanoid } from "nanoid";
import type { Dataset } from "@prisma/client";

import { prisma } from "@/lib/prisma";
import { normalizeS3Url } from "@/lib/utils/viewer-config";
import { fetchJsonFromS3 } from "@/lib/catalog/derive-genes";
import { resolveDatasetOwner } from "@/lib/datasets/owner";

export interface RegisterOrClaimResult {
  dataset: Dataset;
  /** true when we created the row or claimed a previously-unowned one */
  claimed: boolean;
  /** true when the row was already owned by someone other than the actor */
  ownedByOther: boolean;
}

/**
 * Ensure a `Dataset` row exists for a raw S3 URL and is owned (never left
 * unowned). Dedups on the normalized S3 base URL (unique). If a row exists but
 * is unowned, it's claimed for the actor; if already owned, it's returned
 * as-is. New rows read counts/title from the manifest (best-effort).
 *
 * Shared by the register-s3 route and catalog import so ownership + dedup
 * behave identically. Throws on an invalid URL.
 */
export async function ensureDatasetForS3Url(
  url: string,
  {
    userId,
    isAdmin,
    asAdmin,
    kind,
  }: {
    userId: string;
    isAdmin: boolean;
    asAdmin: boolean;
    /**
     * What the caller opened the URL as. Only used when the server can't tell
     * for itself — it won't fetch from a non-public host, e.g. in local tests.
     */
    kind?: "image";
  },
): Promise<RegisterOrClaimResult> {
  const s3BaseUrl = normalizeS3Url(url);

  // Validate — throws for callers to handle.
  new URL(s3BaseUrl);

  const owner = resolveDatasetOwner({ isAdmin, asAdmin, userId });

  const existing = await prisma.dataset.findUnique({ where: { s3BaseUrl } });

  if (existing) {
    const unclaimed = existing.ownerId === null && !existing.adminOwned;

    if (unclaimed) {
      const dataset = await prisma.dataset.update({
        where: { id: existing.id },
        data: { ownerId: owner.ownerId, adminOwned: owner.adminOwned },
      });

      return { dataset, claimed: true, ownedByOther: false };
    }

    // Owned: "mine" if I'm the personal owner, or it's admin-owned and I'm an admin.
    const mine =
      existing.ownerId === userId || (existing.adminOwned && isAdmin);

    return { dataset: existing, claimed: false, ownedByOther: !mine };
  }

  // New registration — read counts/title from the manifest (best-effort).
  const base = s3BaseUrl.replace(/\/+$/, "");
  const manifest = await fetchJsonFromS3(`${base}/manifest.json`);

  // No manifest: an OME-Zarr image is a dataset too, known by the multiscale
  // metadata in its root zarr.json and named after its folder.
  if (!manifest) {
    const zarrRoot = await fetchJsonFromS3(`${base}/zarr.json`);

    if (
      zarrRoot?.attributes?.ome?.multiscales ||
      (!zarrRoot && kind === "image")
    ) {
      const dataset = await prisma.dataset.create({
        data: {
          id: `s3_${nanoid(10)}`,
          fingerprint: createHash("sha256").update(s3BaseUrl).digest("hex"),
          // "experiment/set1_dapi": images in different folders are often
          // named alike, so the folder is part of the name.
          title: decodeURIComponent(
            new URL(base).pathname
              .split("/")
              .filter(Boolean)
              .slice(-2)
              .join("/")
              .replace(/\.zarr$/, ""),
          ),
          numCells: 0,
          numGenes: 0,
          datasetType: "image",
          formatVersion: "zarr",
          status: "COMPLETE",
          ownerId: owner.ownerId,
          adminOwned: owner.adminOwned,
          ingestSource: "s3_registered",
          s3BaseUrl,
          manifestUrl: s3BaseUrl,
        },
      });

      return { dataset, claimed: true, ownedByOther: false };
    }
  }
  const stats = manifest?.statistics ?? {};
  // Single-molecule manifests report total_molecules; single-cell report
  // total_cells. Normalize to the pipeline's datasetType so the viewer routes
  // to /sm-viewer vs /viewer correctly.
  const isSingleMolecule = stats.total_molecules != null;
  const datasetType = isSingleMolecule
    ? "single_molecule"
    : (manifest?.type ?? null);

  const dataset = await prisma.dataset.create({
    data: {
      id: `s3_${nanoid(10)}`,
      fingerprint: createHash("sha256").update(s3BaseUrl).digest("hex"),
      title: manifest?.name ?? "Untitled dataset",
      // total_cells/total_genes for single-cell; total_molecules/unique_genes
      // for single-molecule manifests.
      numCells: Number(stats.total_cells ?? stats.total_molecules) || 0,
      numGenes: Number(stats.total_genes ?? stats.unique_genes) || 0,
      datasetType,
      status: "COMPLETE",
      ownerId: owner.ownerId,
      adminOwned: owner.adminOwned,
      ingestSource: "s3_registered",
      s3BaseUrl,
      manifestUrl: s3BaseUrl,
    },
  });

  return { dataset, claimed: true, ownedByOther: false };
}
