import { prisma } from "@/lib/prisma";
import {
  uploadBufferToS3,
  getObjectJson,
  deleteObject,
} from "@/lib/s3";

/**
 * A single-cell dataset can carry a single-molecule "overlay": a mapping.json
 * in its S3 prefix that tells the viewer to render the molecules of another
 * (app) dataset on top of the cells. This module is the one place that writes,
 * reads, and clears that link, plus the Project that groups the two datasets —
 * shared by the retroactive-overlay API and the duplicate-upload finalize path.
 */

const mappingKey = (scId: string) => `datasets/${scId}/mapping.json`;

/**
 * Build the overlay mapping.json for a single-molecule dataset. How the viewer
 * loads it depends on the dataset's storage:
 *   - app-uploaded (chunked, sm_ id)        → source "app", ref = dataset id
 *     (loaded via SingleMoleculeDataset.fromS3 by id / presigned API)
 *   - S3-registered (has s3BaseUrl, s3_ id) → source "s3", ref = the public URL
 *     (loaded via fromCustomS3) — fromS3(id) would 404 for these.
 * `datasetId` is always stored so the UI can identify the overlaid dataset
 * regardless of source.
 */
async function buildMapping(smId: string) {
  const sm = await prisma.dataset.findUnique({
    where: { id: smId },
    select: { id: true, s3BaseUrl: true, ingestSource: true },
  });
  const isS3Registered =
    sm?.ingestSource === "s3_registered" || (!!sm?.s3BaseUrl && smId.startsWith("s3_"));

  if (isS3Registered && sm?.s3BaseUrl) {
    return {
      linkColumn: "__all__",
      links: { __all__: sm.s3BaseUrl },
      source: "s3" as const,
      datasetId: smId,
    };
  }

  return {
    linkColumn: "__all__",
    links: { __all__: smId },
    source: "app" as const,
    datasetId: smId,
  };
}

/** The single-molecule dataset id currently overlaid on `scId`, or null. */
export async function getOverlaySmId(scId: string): Promise<string | null> {
  const m = await getObjectJson<{
    datasetId?: string;
    source?: string;
    links?: Record<string, string>;
  }>(mappingKey(scId)).catch(() => null);

  if (!m) return null;
  // Prefer the explicit datasetId; fall back to links.__all__ for app-source
  // mappings written before datasetId was stored.
  return m.datasetId ?? (m.source !== "s3" ? (m.links?.__all__ ?? null) : null);
}

/** The image overlay (`image` URL) in the dataset's mapping.json, if any. */
export async function getOverlayImage(
  scId: string,
): Promise<string | undefined> {
  const m = await getObjectJson<{ image?: string }>(mappingKey(scId)).catch(
    () => null,
  );

  return m?.image;
}

// The image overlay shares mapping.json with the molecule overlay but is
// independent of it, so writing or removing the molecules keeps the image.
export async function writeOverlay(scId: string, smId: string): Promise<void> {
  const mapping = await buildMapping(smId);
  const image = await getOverlayImage(scId);

  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify(image ? { ...mapping, image } : mapping)),
    "application/json",
  );
}

export async function removeOverlay(scId: string): Promise<void> {
  const image = await getOverlayImage(scId);

  if (image) {
    await uploadBufferToS3(
      mappingKey(scId),
      Buffer.from(JSON.stringify({ image })),
      "application/json",
    );

    return;
  }
  await deleteObject(mappingKey(scId));
}

/**
 * Set the image overlay (the public URL of an OME-Zarr image) or, with null,
 * clear it. The molecule overlay in the same mapping.json is left as it is.
 */
export async function writeOverlayImage(
  scId: string,
  image: string | null,
): Promise<void> {
  const mapping: Record<string, unknown> =
    (await getObjectJson<Record<string, unknown>>(mappingKey(scId)).catch(
      () => null,
    )) ?? {};

  if (image) mapping.image = image;
  else delete mapping.image;

  if (Object.keys(mapping).length === 0) {
    await deleteObject(mappingKey(scId));

    return;
  }
  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify(mapping)),
    "application/json",
  );
}

/**
 * Group two datasets in a Project owned by `ownerId`. Reuses an existing
 * project that already contains the single-cell dataset if one exists;
 * otherwise creates one titled `title`. Returns the project id.
 */
export async function groupInProject(
  ownerId: string,
  title: string,
  scId: string,
  smId: string,
): Promise<string> {
  const existing = await prisma.projectDataset.findFirst({
    where: { datasetId: scId, project: { ownerId } },
    select: { projectId: true },
  });

  const projectId =
    existing?.projectId ??
    (await prisma.project.create({ data: { title, ownerId } })).id;

  for (const [i, datasetId] of [scId, smId].entries()) {
    await prisma.projectDataset.upsert({
      where: { projectId_datasetId: { projectId, datasetId } },
      update: {},
      create: { projectId, datasetId, sortOrder: i },
    });
  }

  return projectId;
}
