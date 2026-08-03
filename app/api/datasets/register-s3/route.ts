import { NextRequest, NextResponse } from "next/server";

import { requireUser } from "@/lib/admin-auth";
import { ensureDatasetForS3Url } from "@/lib/datasets/register-or-claim";

/**
 * Register (or claim) a dataset that is loaded from a raw S3 URL so it can
 * carry ownership, view counts, and owner-saved viewer defaults. Dedups on the
 * normalized S3 base URL — re-registering the same URL never creates a
 * duplicate. Ownership logic lives in ensureDatasetForS3Url (shared with
 * catalog import).
 */
export async function POST(request: NextRequest) {
  const { error, session } = await requireUser();

  if (error) return error;

  let body: { url?: string; asAdmin?: boolean };

  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ error: "Invalid JSON" }, { status: 400 });
  }

  if (!body.url) {
    return NextResponse.json({ error: "Missing url" }, { status: 400 });
  }

  const isAdmin =
    session.user.role === "ADMIN" || session.user.role === "SUPER_ADMIN";

  try {
    const { dataset, claimed, ownedByOther } = await ensureDatasetForS3Url(
      body.url,
      { userId: session.user.id, isAdmin, asAdmin: !!body.asAdmin },
    );

    return NextResponse.json({ dataset, claimed, ownedByOther });
  } catch {
    return NextResponse.json({ error: "Invalid url" }, { status: 400 });
  }
}
