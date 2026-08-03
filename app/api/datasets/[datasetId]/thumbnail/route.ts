import { NextRequest, NextResponse } from "next/server";
import { revalidatePath } from "next/cache";

import { requireUser, ownerOrAdminError } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { uploadBufferToS3, getObjectBytes } from "@/lib/s3";

const MAX_BYTES = 6 * 1024 * 1024; // 6 MB cap on an uploaded thumbnail

function thumbKey(datasetId: string) {
  return `thumbnails/${datasetId}.jpg`;
}

/**
 * GET — serve the dataset's thumbnail image (public; used in <img src>). Proxies
 * the object from the app's S3 bucket so it works regardless of bucket privacy.
 */
export async function GET(
  _req: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { datasetId } = await params;
  const obj = await getObjectBytes(thumbKey(datasetId));

  if (!obj) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  return new Response(new Uint8Array(obj.body), {
    headers: {
      "Content-Type": obj.contentType || "image/jpeg",
      // Short cache; the stored thumbnailUrl carries a ?v= cache-buster on change.
      "Cache-Control": "public, max-age=300",
    },
  });
}

/**
 * POST — owner/admin uploads a screenshot as the dataset's thumbnail. Body is
 * the raw image bytes (image/jpeg). Stores it in S3 and sets Dataset.thumbnailUrl
 * to this route (with a cache-busting version).
 */
export async function POST(
  request: NextRequest,
  { params }: { params: Promise<{ datasetId: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { datasetId } = await params;
  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: { id: true, ownerId: true, adminOwned: true, s3BaseUrl: true },
  });

  if (!dataset) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  const forbidden = ownerOrAdminError(dataset, session);

  if (forbidden) return forbidden;

  const buf = Buffer.from(await request.arrayBuffer());

  if (buf.byteLength === 0) {
    return NextResponse.json({ error: "Empty image" }, { status: 400 });
  }
  if (buf.byteLength > MAX_BYTES) {
    return NextResponse.json({ error: "Image too large" }, { status: 413 });
  }

  await uploadBufferToS3(thumbKey(datasetId), buf, "image/jpeg");

  const thumbnailUrl = `/api/datasets/${datasetId}/thumbnail?v=${Date.now()}`;

  await prisma.dataset.update({
    where: { id: datasetId },
    data: { thumbnailUrl },
  });

  // Propagate into the explore catalog so cards + detail pages reflect it.
  // The catalog entry's s3BaseUrl may carry a trailing slash the Dataset's
  // normalized one lacks, so match both variants.
  if (dataset.s3BaseUrl) {
    const variants = [dataset.s3BaseUrl, `${dataset.s3BaseUrl}/`];
    const entries = await prisma.catalogDatasetEntry.findMany({
      where: { s3BaseUrl: { in: variants } },
      select: { id: true, catalogId: true },
    });

    if (entries.length > 0) {
      const catalogIds = [...new Set(entries.map((e) => e.catalogId))];

      await prisma.catalogDatasetEntry.updateMany({
        where: { id: { in: entries.map((e) => e.id) } },
        data: { thumbnailUrl },
      });
      await prisma.catalogDataset.updateMany({
        where: { id: { in: catalogIds } },
        data: { thumbnailUrl },
      });

      // Bust the ISR cache so explore reflects the new thumbnail immediately.
      const cats = await prisma.catalogDataset.findMany({
        where: { id: { in: catalogIds } },
        select: { id: true, bilCode: true },
      });

      revalidatePath("/explore");
      for (const c of cats) {
        revalidatePath(`/explore/${c.id}`);
        if (c.bilCode) revalidatePath(`/explore/bil/${c.bilCode}`);
      }
    }
  }

  return NextResponse.json({ thumbnailUrl });
}
