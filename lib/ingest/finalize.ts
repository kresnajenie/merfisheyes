import type { Prisma } from "@prisma/client";

import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { getObjectJson, listObjectKeys, uploadBufferToS3 } from "@/lib/s3";
import { sendDatasetReadyEmail } from "@/lib/ses";
import { imageOverlayFor } from "@/lib/ingest/overlay";
import { updateDatasetGenes } from "@/lib/datasets/read-genes";
import { usesPreparePipeline } from "@/lib/ingest/pipeline";

/**
 * Mark a dataset COMPLETE and register the worker's output.
 *
 * Shared by two callers that must not drift apart: the worker's callback, and
 * the /status reconciliation that recovers a job whose callback never arrived.
 * Registering the output is not optional — the viewer resolves files from the
 * most recent UploadSession's rows, and the raw upload's rows describe the
 * *input*, so without this /viewer cannot find manifest.json.
 */

export interface FinalizeStats {
  numCells?: number;
  numGenes?: number;
}

export interface FinalizeResult {
  registeredFiles: number;
  emailed: boolean;
}

export async function finalizeCompletedDataset(
  datasetId: string,
  opts: {
    stats?: FinalizeStats;
    fingerprint?: string;
    sendEmail?: boolean;
  } = {},
): Promise<FinalizeResult> {
  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    include: { owner: true },
  });

  if (!dataset) throw new Error(`Dataset ${datasetId} not found`);

  const outputFiles = await listObjectKeys(`datasets/${datasetId}/`);

  if (outputFiles.length === 0) {
    throw new Error(`No output objects under datasets/${datasetId}/`);
  }

  const uploadId = `up_out${nanoid(8)}`;
  // scripts/prepare.py wrote one SpatialData store: it opens in the
  // SpatialData viewer, whatever card it was uploaded from
  const prepared = usesPreparePipeline(dataset.processingParams);

  // Single-molecule: the worker's manifest (gene list, per-gene counts,
  // statistics) lives in S3; mirror it onto the row like a browser upload.
  const smManifest =
    dataset.datasetType === "single_molecule" && !prepared
      ? await getObjectJson<Record<string, unknown>>(
          `datasets/${datasetId}/manifest.json.gz`,
        ).catch(() => null)
      : null;

  await prisma.$transaction(
    async (tx) => {
      // Must be idempotent: the worker retries failed callbacks, Batch can
      // retry a job, and reconciliation may race a late-arriving callback.
      // Drop any previous worker-generated output registry (cascading its
      // UploadFile rows) so a re-delivery can't leave duplicate sessions.
      await tx.uploadSession.deleteMany({
        where: { datasetId, id: { startsWith: "up_out" } },
      });
      await tx.uploadSession.create({
        data: {
          id: uploadId,
          datasetId,
          totalFiles: outputFiles.length,
          completedFiles: outputFiles.length,
          // Output files are permanent; this session is a registry, not a
          // time-boxed upload window.
          expiresAt: new Date(Date.now() + 365 * 24 * 3600 * 1000),
        },
      });
      await tx.uploadFile.createMany({
        data: outputFiles.map((f) => ({
          uploadSessionId: uploadId,
          fileKey: f.key,
          fileSize: BigInt(f.size),
          status: "COMPLETE" as const,
          uploadedAt: new Date(),
        })),
      });
      await tx.dataset.update({
        where: { id: datasetId },
        data: {
          status: "COMPLETE",
          completedAt: new Date(),
          numCells: opts.stats?.numCells ?? dataset.numCells,
          numGenes: opts.stats?.numGenes ?? dataset.numGenes,
          ...(prepared
            ? { datasetType: "spatialdata", formatVersion: "zarr" }
            : {}),
          // Browser uploads store the SM manifest at initiate; server
          // processing writes it to S3, so copy it onto the row here — the
          // SM API and gene indexing read it from the row.
          ...(smManifest
            ? { manifestJson: smManifest as Prisma.InputJsonValue }
            : {}),
          errorMessage: null,
          processingProgress: {
            stage: "Done",
            percent: 100,
            updatedAt: new Date().toISOString(),
          },
        },
      });
    },
    { timeout: 60000 },
  );

  // If this SC dataset was uploaded alongside a single-molecule dataset and/or
  // an image (the combined upload flow), write a mapping.json so the SC viewer
  // draws the molecules over the cells and the image under them. Written
  // post-processing so the worker's output can't clobber it.
  const linked = dataset.processingParams as {
    linkedSmDatasetId?: unknown;
    linkedImageDatasetId?: unknown;
  } | null;
  const linkedSmId = linked?.linkedSmDatasetId;
  const linkedImageId = linked?.linkedImageDatasetId;

  try {
    const mapping = {
      ...(typeof linkedSmId === "string" && linkedSmId
        ? { linkColumn: "__all__", links: { __all__: linkedSmId }, source: "app" }
        : {}),
      ...(typeof linkedImageId === "string" && linkedImageId
        ? await imageOverlayFor(linkedImageId)
        : {}),
    };

    // prepare.py's own mapping.json names the store; it is not overwritten
    if (!prepared && Object.keys(mapping).length > 0) {
      await uploadBufferToS3(
        `datasets/${datasetId}/mapping.json`,
        Buffer.from(JSON.stringify(mapping)),
        "application/json",
      );
    }
  } catch (e) {
    // eslint-disable-next-line no-console
    console.error(
      "Ingest finalize: mapping.json write failed:",
      (e as Error).message,
    );
  }

  // The worker's hash of the raw upload replaces the raw_ placeholder. It is a
  // record of what was uploaded, never checked against other datasets: the
  // same files uploaded twice make two datasets (docs/plans/PREVIEW-SAVE.md,
  // D4). The column is indexed, not unique.
  if (opts.fingerprint && opts.fingerprint !== dataset.fingerprint) {
    await prisma.dataset.update({
      where: { id: datasetId },
      data: { fingerprint: opts.fingerprint },
    });
  }

  // Index the dataset's genes for gene search (best-effort — reads the
  // processed output from S3; it never throws, so awaiting can't fail the
  // finalize). Must be awaited: a fire-and-forget promise is killed when the
  // serverless route returns, leaving the row without genes.
  await updateDatasetGenes(prisma, {
    id: datasetId,
    datasetType: prepared ? "spatialdata" : dataset.datasetType,
    formatVersion: prepared ? "zarr" : dataset.formatVersion,
    manifestJson: dataset.manifestJson,
  });

  // Never fail a finalize because SES failed — the dataset is complete and
  // viewable either way.
  let emailed = false;

  if (opts.sendEmail !== false && dataset.owner?.email) {
    try {
      await sendDatasetReadyEmail({
        email: dataset.owner.email,
        datasetId,
        datasetType: prepared ? "spatialdata" : dataset.datasetType,
        // A molecule table prepared into a store still counts molecules
        molecules: prepared
          ? (dataset.processingParams as { kind?: unknown } | null)?.kind ===
            "single_molecule"
          : dataset.datasetType === "single_molecule",
        datasetName: dataset.title ?? undefined,
        metadata: {
          numCells: opts.stats?.numCells ?? undefined,
          numGenes: opts.stats?.numGenes ?? undefined,
          platform:
            (prepared ? "spatialdata" : dataset.datasetType) ?? undefined,
        },
      });
      emailed = true;
    } catch (e: any) {
      console.error("Ingest finalize: sending email failed:", e?.message);
    }
  }

  return { registeredFiles: outputFiles.length, emailed };
}
