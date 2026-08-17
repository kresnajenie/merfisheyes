import { NextRequest, NextResponse } from "next/server";

import { requireUser, ownerOrAdminError } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";
import { validateViewerConfig } from "@/lib/utils/viewer-config";
import {
  sanitizeMetadataPatch,
  sanitizeMetadataBag,
} from "@/lib/datasets/metadata";

/**
 * Owner (or admin) edit of a dataset's metadata and viewer defaults. Single
 * save path for the account page (title) and the in-viewer "Save current view
 * as default" (viewerConfig). Gated: only the owner or an admin may write.
 *
 * NOTE: for datasets in our own S3 bucket we could additionally bake the
 * config into the S3 palette/coord files (permanent). That write-back is not
 * yet implemented; the DB `viewerConfig` override applies on load for all
 * datasets in the meantime.
 */
export async function PATCH(
  request: NextRequest,
  { params }: { params: Promise<{ id: string }> },
) {
  const { error, session } = await requireUser();

  if (error) return error;

  const { id } = await params;
  const dataset = await prisma.dataset.findUnique({
    where: { id },
    select: { id: true, ownerId: true, adminOwned: true },
  });

  if (!dataset) {
    return NextResponse.json({ error: "Not found" }, { status: 404 });
  }

  const forbidden = ownerOrAdminError(dataset, session);

  if (forbidden) return forbidden;

  let body: Record<string, unknown>;

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  const data: Record<string, unknown> = {};

  if (typeof body.title === "string") {
    data.title = body.title.trim().slice(0, 500);
  }

  if (typeof body.thumbnailUrl === "string" || body.thumbnailUrl === null) {
    data.thumbnailUrl =
      typeof body.thumbnailUrl === "string"
        ? body.thumbnailUrl.trim().slice(0, 1000) || null
        : null;
  }

  if (body.viewerConfig !== undefined) {
    // null clears the saved defaults; otherwise sanitize.
    data.viewerConfig =
      body.viewerConfig === null
        ? null
        : (validateViewerConfig(body.viewerConfig) ?? undefined);
  }

  // Catalog-style metadata (species, tissue, description, tags, links, …).
  Object.assign(data, sanitizeMetadataPatch(body));

  // Flexible display-metadata bag (investigator, authors, citation, …).
  if ("metadata" in body) {
    data.metadata = sanitizeMetadataBag(body.metadata);
  }

  if (Object.keys(data).length === 0) {
    return NextResponse.json({ error: "Nothing to update" }, { status: 400 });
  }

  const updated = await prisma.dataset.update({
    where: { id },
    data,
    select: { id: true, title: true, viewerConfig: true },
  });

  return NextResponse.json({ dataset: updated });
}
