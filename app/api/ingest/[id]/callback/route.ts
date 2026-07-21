import { createHmac, timingSafeEqual } from "crypto";

import { NextRequest, NextResponse } from "next/server";
import { nanoid } from "nanoid";

import { prisma } from "@/lib/prisma";
import { listObjectKeys } from "@/lib/s3";
import { sendDatasetReadyEmail } from "@/lib/ses";

// Worker → app callback. This is the ONLY writer of ingestion status, so DB and
// SES credentials never leave the app (design §5.4). Authenticated with an HMAC
// of the raw body using CALLBACK_SECRET — the worker only needs that secret and
// outbound HTTPS.

type CallbackStatus = "PROCESSING" | "COMPLETE" | "FAILED";

interface CallbackBody {
  status: CallbackStatus;
  /** Staged progress for the UI, e.g. { stage: "Expression chunks", percent: 40 } */
  progress?: { stage?: string; percent?: number; message?: string };
  stats?: { numCells?: number; numGenes?: number; spatialDimensions?: number };
  manifestKey?: string;
  fingerprint?: string;
  error?: string;
}

function verifySignature(rawBody: string, header: string | null): boolean {
  const secret = process.env.CALLBACK_SECRET;

  if (!secret || !header) return false;

  const provided = header.startsWith("sha256=") ? header.slice(7) : header;
  const expected = createHmac("sha256", secret).update(rawBody).digest("hex");

  const a = Buffer.from(provided, "hex");
  const b = Buffer.from(expected, "hex");

  // timingSafeEqual throws on length mismatch, so guard first.
  return a.length === b.length && timingSafeEqual(a, b);
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  try {
    const { id: datasetId } = await params;

    if (!process.env.CALLBACK_SECRET) {
      console.error("Ingest callback: CALLBACK_SECRET is not configured");

      return NextResponse.json(
        { error: "Callback not configured" },
        { status: 503 },
      );
    }

    // Read the raw body for signature verification, then parse.
    const rawBody = await request.text();

    if (!verifySignature(rawBody, request.headers.get("x-ingest-signature"))) {
      return NextResponse.json({ error: "Invalid signature" }, { status: 401 });
    }

    let body: CallbackBody;

    try {
      body = JSON.parse(rawBody);
    } catch {
      return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
    }

    const dataset = await prisma.dataset.findUnique({
      where: { id: datasetId },
      include: { owner: true },
    });

    if (!dataset) {
      return NextResponse.json({ error: "Dataset not found" }, { status: 404 });
    }

    // ── PROCESSING / progress ──────────────────────────────────────────
    if (body.status === "PROCESSING") {
      await prisma.dataset.update({
        where: { id: datasetId },
        data: {
          status: "PROCESSING",
          processingProgress: {
            ...(body.progress ?? {}),
            updatedAt: new Date().toISOString(),
          },
        },
      });

      return NextResponse.json({ success: true, status: "PROCESSING" });
    }

    // ── FAILED ─────────────────────────────────────────────────────────
    if (body.status === "FAILED") {
      await prisma.dataset.update({
        where: { id: datasetId },
        data: {
          status: "FAILED",
          errorMessage: body.error?.slice(0, 2000) || "Processing failed",
          processingProgress: {
            ...(body.progress ?? {}),
            failedAt: new Date().toISOString(),
          },
        },
      });

      return NextResponse.json({ success: true, status: "FAILED" });
    }

    // ── COMPLETE ───────────────────────────────────────────────────────
    if (body.status !== "COMPLETE") {
      return NextResponse.json(
        { error: `Unknown status: ${body.status}` },
        { status: 400 },
      );
    }

    // The viewer resolves files from the most recent UploadSession's rows, so
    // the worker's OUTPUT must be registered here — otherwise /viewer cannot
    // find manifest.json (the raw upload's rows describe the input file).
    const outputFiles = await listObjectKeys(`datasets/${datasetId}/`);

    if (outputFiles.length === 0) {
      return NextResponse.json(
        { error: `No output objects under datasets/${datasetId}/` },
        { status: 400 },
      );
    }

    const uploadId = `up_out${nanoid(8)}`;

    await prisma.$transaction(
      async (tx) => {
        // COMPLETE must be idempotent: the worker retries failed callbacks, and
        // Batch can retry a job. Drop any previous worker-generated output
        // registry (cascades its UploadFile rows) before writing the new one,
        // so a re-delivery can't leave duplicate sessions behind.
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
            numCells: body.stats?.numCells ?? dataset.numCells,
            numGenes: body.stats?.numGenes ?? dataset.numGenes,
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

    // The canonical fingerprint replaces the raw_ placeholder. It is @unique, so
    // a collision means an identical dataset already exists — keep the
    // placeholder rather than failing the whole callback (dedup policy TBD).
    if (body.fingerprint && body.fingerprint !== dataset.fingerprint) {
      try {
        await prisma.dataset.update({
          where: { id: datasetId },
          data: { fingerprint: body.fingerprint },
        });
      } catch (e: any) {
        if (e?.code === "P2002") {
          console.warn(
            `Ingest callback: fingerprint ${body.fingerprint} already exists; keeping placeholder for ${datasetId}`,
          );
        } else {
          throw e;
        }
      }
    }

    // Email the uploader. Never fail the callback because SES failed —
    // the dataset is complete and viewable either way.
    let emailed = false;

    if (dataset.owner?.email) {
      try {
        await sendDatasetReadyEmail({
          email: dataset.owner.email,
          datasetId,
          datasetName: dataset.title ?? undefined,
          metadata: {
            numCells: body.stats?.numCells ?? undefined,
            numGenes: body.stats?.numGenes ?? undefined,
            platform: dataset.datasetType ?? undefined,
          },
        });
        emailed = true;
      } catch (e: any) {
        console.error("Ingest callback: sending email failed:", e?.message);
      }
    }

    return NextResponse.json({
      success: true,
      status: "COMPLETE",
      registeredFiles: outputFiles.length,
      emailed,
    });
  } catch (err: any) {
    console.error("Ingest callback error:", err);

    return NextResponse.json(
      { error: "Internal server error", message: err.message },
      { status: 500 },
    );
  }
}
