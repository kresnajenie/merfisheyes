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

/** The full overlay mapping.json the viewer expects for an app SM dataset. */
export function overlayMapping(smId: string) {
  return { linkColumn: "__all__", links: { __all__: smId }, source: "app" as const };
}

/** The sm_ dataset id currently overlaid on `scId`, or null. */
export async function getOverlaySmId(scId: string): Promise<string | null> {
  const m = await getObjectJson<{ links?: Record<string, string> }>(
    mappingKey(scId),
  ).catch(() => null);

  return m?.links?.__all__ ?? null;
}

export async function writeOverlay(scId: string, smId: string): Promise<void> {
  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify(overlayMapping(smId))),
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
