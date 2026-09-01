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

export async function writeOverlay(scId: string, smId: string): Promise<void> {
  const mapping = await buildMapping(smId);

  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify(mapping)),
    "application/json",
  );
}

export async function removeOverlay(scId: string): Promise<void> {
  await deleteObject(mappingKey(scId));
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
