import { TerminateJobCommand } from "@aws-sdk/client-batch";
import { Prisma } from "@prisma/client";

import { batchClient } from "@/lib/batch";
import { getOverlaySmId, removeOverlay } from "@/lib/ingest/overlay";
import { prisma } from "@/lib/prisma";
import {
  abortMultipartUpload,
  deleteObject,
  deleteObjectsByPrefix,
  listInProgressMultipartUploads,
} from "@/lib/s3";

export interface DeletableDataset {
  id: string;
  title: string | null;
  ownerId: string | null;
  adminOwned: boolean;
  datasetType: string | null;
  status: string;
  batchJobId: string | null;
  processingParams: Prisma.JsonValue | null;
}

export const deletableDatasetSelect = {
  id: true,
  title: true,
  ownerId: true,
  adminOwned: true,
  datasetType: true,
  status: true,
  batchJobId: true,
  processingParams: true,
} as const;

const linkedSmId = (params: Prisma.JsonValue | null): string | null => {
  const id = (params as { linkedSmDatasetId?: unknown } | null)
    ?.linkedSmDatasetId;

  return typeof id === "string" && id ? id : null;
};

/** The single-molecule dataset overlaid on a single-cell dataset, or null. */
export async function overlayOf(
  dataset: DeletableDataset,
): Promise<string | null> {
  if (dataset.datasetType === "single_molecule") return null;

  return (
    (await getOverlaySmId(dataset.id)) ?? linkedSmId(dataset.processingParams)
  );
}

/**
 * The Explore card (curated or community, any review state) that shows this
 * dataset, if there is one. Deleting the dataset would leave that card
 * pointing at nothing, so it has to be withdrawn first.
 */
export function exploreCardFor(datasetId: string) {
  return prisma.catalogDataset.findFirst({
    where: {
      OR: [
        { sourceDatasetId: datasetId },
        { entries: { some: { datasetId } } },
      ],
    },
    select: { id: true, title: true },
  });
}

/**
 * Stop single-cell datasets from overlaying a molecule dataset that is about
 * to go: clear their mapping.json and the combined-upload link the mapping
 * route falls back to. Only the same owner's datasets can link to it.
 */
async function unlinkOverlay(sm: DeletableDataset): Promise<string[]> {
  const candidates = await prisma.dataset.findMany({
    where: {
      ...(sm.ownerId ? { ownerId: sm.ownerId } : { adminOwned: true }),
      NOT: { datasetType: "single_molecule" },
    },
    select: { id: true, processingParams: true },
  });
  const unlinked: string[] = [];

  await Promise.all(
    candidates.map(async (sc) => {
      const byParam = linkedSmId(sc.processingParams) === sm.id;
      const byMapping = (await getOverlaySmId(sc.id)) === sm.id;

      if (!byParam && !byMapping) return;
      if (byMapping) await removeOverlay(sc.id);
      if (byParam) {
        const { linkedSmDatasetId: _, ...rest } = sc.processingParams as Record<
          string,
          Prisma.JsonValue
        >;

        await prisma.dataset.update({
          where: { id: sc.id },
          data: { processingParams: rest as Prisma.InputJsonObject },
        });
      }
      unlinked.push(sc.id);
    }),
  );

  return unlinked;
}

/**
 * Permanently delete one dataset: its processing job, everything it has in
 * our bucket, and its row (which cascades to upload sessions, views and
 * project memberships). The row goes last, so a failure part-way leaves it
 * in place to retry. Data in someone else's bucket (S3-registered datasets)
 * is not ours to delete and is left alone.
 */
export async function deleteDatasetCompletely(
  dataset: DeletableDataset,
): Promise<{ unlinkedFrom: string[] }> {
  if (
    dataset.batchJobId &&
    (dataset.status === "QUEUED" || dataset.status === "PROCESSING")
  ) {
    await batchClient
      .send(
        new TerminateJobCommand({
          jobId: dataset.batchJobId,
          reason: "Dataset deleted by its owner",
        }),
      )
      .catch((e) =>
        console.warn(
          `Delete ${dataset.id}: could not stop its job:`,
          e.message,
        ),
      );
  }

  const unlinkedFrom =
    dataset.datasetType === "single_molecule"
      ? await unlinkOverlay(dataset)
      : [];

  // Incomplete multipart uploads hold staged parts that a prefix delete
  // doesn't see; the bucket's lifecycle rule is the backstop if this fails.
  const rawPrefix = `raw/${dataset.id}/`;

  await listInProgressMultipartUploads(rawPrefix)
    .then((uploads) =>
      Promise.all(uploads.map((u) => abortMultipartUpload(u.key, u.uploadId))),
    )
    .catch((e) =>
      console.warn(
        `Delete ${dataset.id}: could not abort multipart uploads:`,
        e.message,
      ),
    );

  await deleteObjectsByPrefix(rawPrefix);
  await deleteObjectsByPrefix(`datasets/${dataset.id}/`);
  await deleteObject(`thumbnails/${dataset.id}.jpg`);

  await prisma.dataset.delete({ where: { id: dataset.id } });

  return { unlinkedFrom };
}
