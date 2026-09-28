// app/api/datasets/[datasetId]/upload-policy/route.ts
// Re-signs the prefix-scoped POST policy for an in-progress chunked upload,
// so uploads that outlive the policy's 1-hour expiry can keep going.
import { NextRequest, NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";
import { auth } from "@/lib/auth";
import { generatePrefixUploadPolicy } from "@/lib/s3";

const corsHeaders = {
  "Access-Control-Allow-Origin": process.env.CORS_ORIGIN || "*",
  "Access-Control-Allow-Methods": "POST, OPTIONS",
  "Access-Control-Allow-Headers": "Content-Type",
};

export async function OPTIONS() {
  return NextResponse.json({}, { headers: corsHeaders });
}

export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  try {
    const { datasetId } = await params;
    const { uploadId } = await request.json();

    if (!uploadId) {
      return NextResponse.json(
        { error: "uploadId is required" },
        { status: 400, headers: corsHeaders },
      );
    }

    const session = await auth();

    if (!session?.user) {
      return NextResponse.json(
        { error: "Sign in to upload" },
        { status: 401, headers: corsHeaders },
      );
    }

    const uploadSession = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      include: { dataset: true },
    });

    if (!uploadSession || uploadSession.datasetId !== datasetId) {
      return NextResponse.json(
        { error: "Upload session not found" },
        { status: 404, headers: corsHeaders },
      );
    }

    const { dataset } = uploadSession;
    const isAdmin =
      session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
    const owned =
      dataset.ownerId === session.user.id || (isAdmin && dataset.adminOwned);

    if (!owned) {
      return NextResponse.json(
        { error: "Forbidden" },
        { status: 403, headers: corsHeaders },
      );
    }

    if (dataset.status !== "UPLOADING") {
      return NextResponse.json(
        { error: "Dataset is no longer accepting uploads" },
        { status: 409, headers: corsHeaders },
      );
    }

    // Same cap as initiate: the largest file registered for this session.
    const { _max } = await prisma.uploadFile.aggregate({
      where: { uploadSessionId: uploadId },
      _max: { fileSize: true },
    });
    const maxFileSize = Math.max(1, Number(_max.fileSize ?? 0));

    const policy = await generatePrefixUploadPolicy(
      `datasets/${datasetId}/`,
      maxFileSize,
    );

    await prisma.uploadSession.update({
      where: { id: uploadId },
      data: { expiresAt: new Date(policy.expiresAt) },
    });

    return NextResponse.json(
      { success: true, policy },
      { headers: corsHeaders },
    );
  } catch (error: any) {
    console.error("Refresh upload policy error:", error);

    return NextResponse.json(
      { error: "Internal server error", message: error.message },
      { status: 500, headers: corsHeaders },
    );
  }
}
