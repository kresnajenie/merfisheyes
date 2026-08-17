import { createHmac, timingSafeEqual } from "crypto";

import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { finalizeCompletedDataset } from "@/lib/ingest/finalize";
import { markDatasetFailed } from "@/lib/ingest/mark-failed";

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

    // Existence check only — finalizeCompletedDataset re-reads what it needs,
    // including the owner for the notification email.
    const exists = await prisma.dataset.findUnique({
      where: { id: datasetId },
      select: { id: true },
    });

    if (!exists) {
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
      // Classify the fault + email the owner, via the shared helper the
      // reconcile path also uses.
      await markDatasetFailed(
        datasetId,
        body.error?.slice(0, 2000) || "Processing failed",
        { progress: body.progress },
      );

      return NextResponse.json({ success: true, status: "FAILED" });
    }

    // ── COMPLETE ───────────────────────────────────────────────────────
    if (body.status !== "COMPLETE") {
      return NextResponse.json(
        { error: `Unknown status: ${body.status}` },
        { status: 400 },
      );
    }

    // Registering output, updating the row, the fingerprint swap and the email
    // all live in finalizeCompletedDataset — /status reconciliation performs the
    // exact same completion when a callback is lost, and the two must not drift.
    const { registeredFiles, emailed } = await finalizeCompletedDataset(
      datasetId,
      { stats: body.stats, fingerprint: body.fingerprint },
    );

    return NextResponse.json({
      success: true,
      status: "COMPLETE",
      registeredFiles,
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
