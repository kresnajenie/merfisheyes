import type { Prisma } from "@prisma/client";

import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { getObjectJson, listObjectKeys, uploadBufferToS3 } from "@/lib/s3";
import { sendDatasetReadyEmail, sendDuplicateDatasetEmail } from "@/lib/ses";
import { updateDatasetGenes } from "@/lib/datasets/read-genes";

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

  // Dedup: the worker sends a SHA-256 hex digest of the raw upload (64 hex
  // chars; placeholders are "raw_"/"STUB_"-prefixed). If another dataset already
  // holds it, this is a re-upload of the identical file(s) — mark it FAILED
  // pointing at the original rather than minting a second copy. (The fingerprint
  // column is @unique, so at most one dataset can hold each value.)
  //
  // Admins are exempt: they're allowed to upload duplicates (e.g. re-processing,
  // fixtures), so their re-upload completes normally. Its fingerprint will
  // collide on the @unique column at the end, which the tail catch keeps as a
  // placeholder — harmless.
  const isAdmin =
    dataset.owner?.role === "ADMIN" || dataset.owner?.role === "SUPER_ADMIN";

  if (
    !isAdmin &&
    opts.fingerprint &&
    /^[0-9a-f]{64}$/.test(opts.fingerprint) &&
    opts.fingerprint !== dataset.fingerprint
  ) {
    const existing = await prisma.dataset.findFirst({
      where: { fingerprint: opts.fingerprint, id: { not: datasetId } },
      select: { id: true, title: true },
    });

    if (existing) {
      await prisma.dataset.update({
        where: { id: datasetId },
        data: {
          status: "FAILED",
          completedAt: new Date(),
          faultCategory: "USER",
          errorMessage:
            `Duplicate of already-ingested dataset ${existing.id}` +
            (existing.title ? ` ("${existing.title}")` : "") +
            ". Open that one instead of re-uploading.",
        },
      });
      console.warn(
        `Ingest finalize: ${datasetId} duplicates ${existing.id}; marked FAILED.`,
      );

      // Tell the owner it was a duplicate and where the original lives. Never
      // fail the finalize over a mail error.
      let emailed = false;

      if (dataset.owner?.email) {
        try {
          await sendDuplicateDatasetEmail({
            email: dataset.owner.email,
            uploadedName: dataset.title ?? undefined,
            existingId: existing.id,
            existingTitle: existing.title ?? undefined,
          });
          emailed = true;
        } catch (e: any) {
          console.error("Ingest finalize: duplicate email failed:", e?.message);
        }
      }

      return { registeredFiles: 0, emailed };
    }
  }

  const outputFiles = await listObjectKeys(`datasets/${datasetId}/`);

  if (outputFiles.length === 0) {
    throw new Error(`No output objects under datasets/${datasetId}/`);
  }

  const uploadId = `up_out${nanoid(8)}`;

  // Single-molecule: the worker's manifest (gene list, per-gene counts,
  // statistics) lives in S3; mirror it onto the row like a browser upload.
  const smManifest =
    dataset.datasetType === "single_molecule"
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

  // If this SC dataset was uploaded alongside a single-molecule dataset (the
  // combined SC+SM upload flow), write a mapping.json so the SC viewer overlays
  // the molecules. Written post-processing so the worker's output can't clobber it.
  const linkedSmId = (
    dataset.processingParams as { linkedSmDatasetId?: unknown } | null
  )?.linkedSmDatasetId;

  if (typeof linkedSmId === "string" && linkedSmId) {
    try {
      const mapping = {
        linkColumn: "__all__",
        links: { __all__: linkedSmId },
        source: "app",
      };

      await uploadBufferToS3(
        `datasets/${datasetId}/mapping.json`,
        Buffer.from(JSON.stringify(mapping)),
        "application/json",
      );
    } catch (e) {
      // eslint-disable-next-line no-console
      console.error(
        "Ingest finalize: mapping.json write failed:",
        (e as Error).message,
      );
    }
  }

  // The canonical fingerprint replaces the raw_ placeholder. It is @unique, so
  // a collision means an identical dataset already exists — keep the
  // placeholder rather than failing the whole finalize (dedup policy TBD).
  if (opts.fingerprint && opts.fingerprint !== dataset.fingerprint) {
    try {
      await prisma.dataset.update({
        where: { id: datasetId },
        data: { fingerprint: opts.fingerprint },
      });
    } catch (e: any) {
      if (e?.code === "P2002") {
        console.warn(
          `Ingest finalize: fingerprint ${opts.fingerprint} already exists; keeping placeholder for ${datasetId}`,
        );
      } else {
        throw e;
      }
    }
  }

  // Index the dataset's genes for gene search (best-effort — reads the
  // processed output from S3; it never throws, so awaiting can't fail the
  // finalize). Must be awaited: a fire-and-forget promise is killed when the
  // serverless route returns, leaving the row without genes.
  await updateDatasetGenes(prisma, {
    id: datasetId,
    datasetType: dataset.datasetType,
    formatVersion: dataset.formatVersion,
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
        datasetName: dataset.title ?? undefined,
        metadata: {
          numCells: opts.stats?.numCells ?? undefined,
          numGenes: opts.stats?.numGenes ?? undefined,
          platform: dataset.datasetType ?? undefined,
        },
      });
      emailed = true;
    } catch (e: any) {
      console.error("Ingest finalize: sending email failed:", e?.message);
    }
  }

  return { registeredFiles: outputFiles.length, emailed };
}
