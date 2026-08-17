import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { uploadBufferToS3, getObjectBytes } from "@/lib/s3";

const MAX_BYTES = 6 * 1024 * 1024; // 6 MB cap on an uploaded thumbnail

function thumbKey(projectId: string) {
  return `thumbnails/projects/${projectId}.jpg`;
}

function isAdmin(session: { user: { role?: string | null } }): boolean {
  return session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";
}

/**
 * GET — serve the project's uploaded thumbnail image (public; used in <img src>).
 * Proxies the object from the app's S3 bucket so it works regardless of bucket
 * privacy. Only relevant when the owner uploaded a file; a project that instead
 * reuses a dataset's thumbnail stores that dataset's URL directly.
 */
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { id } = await params;
  const obj = await getObjectBytes(thumbKey(id));

  if (!obj) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return new Response(new Uint8Array(obj.body), {
    headers: {
      "Content-Type": obj.contentType || "image/jpeg",
      "Cache-Control": "public, max-age=300",
    },
  });
}

/**
 * POST — owner/admin uploads an image as the project's thumbnail. Body is the
 * raw image bytes. Stores it in S3 and sets Project.thumbnailUrl to this route
 * (with a cache-busting version).
 */
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { id } = await params;
  const project = await prisma.project.findUnique({
    where: { id },
    select: { id: true, ownerId: true },
  });

  if (!project) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }
  if (project.ownerId !== session.user.id && !isAdmin(session)) {
    return NextResponse.json({ error: "Forbidden" }, { status: 403 });
  }

  const buf = Buffer.from(await request.arrayBuffer());

  if (buf.byteLength === 0) {
    return NextResponse.json({ error: "Empty image" }, { status: 400 });
  }
  if (buf.byteLength > MAX_BYTES) {
    return NextResponse.json({ error: "Image too large" }, { status: 413 });
  }

  const contentType = request.headers.get("content-type") || "image/jpeg";

  await uploadBufferToS3(thumbKey(id), buf, contentType);

  const thumbnailUrl = `/api/projects/${id}/thumbnail?v=${Date.now()}`;

  await prisma.project.update({ where: { id }, data: { thumbnailUrl } });

  // Keep a live community submission's card in sync.
  await prisma.catalogDataset.updateMany({
    where: { isCommunity: true, sourceProjectId: id },
    data: { thumbnailUrl },
  });

  return NextResponse.json({ thumbnailUrl });
}
