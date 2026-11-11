import { auth } from "@/lib/auth";
import { prisma } from "@/lib/prisma";

/**
 * Verify the current request is from any authenticated user.
 * Returns the session if valid, or a 401 Response to send back.
 */
export async function requireUser() {
  const session = await auth();

  if (!session?.user) {
    return { error: new Response("Unauthorized", { status: 401 }), session: null };
  }

  return { error: null, session };
}

/**
 * Verify the current request is from an authenticated ADMIN or SUPER_ADMIN user.
 * Returns the session if valid, or a Response to send back.
 */
export async function requireAdmin() {
  const session = await auth();

  if (!session?.user) {
    return { error: new Response("Unauthorized", { status: 401 }), session: null };
  }

  if (session.user.role !== "ADMIN" && session.user.role !== "SUPER_ADMIN") {
    return { error: new Response("Forbidden", { status: 403 }), session: null };
  }

  return { error: null, session };
}

/**
 * Verify the current request is from an authenticated SUPER_ADMIN user.
 * Returns the session if valid, or a Response to send back.
 */
export async function requireSuperAdmin() {
  const session = await auth();

  if (!session?.user) {
    return { error: new Response("Unauthorized", { status: 401 }), session: null };
  }

  if (session.user.role !== "SUPER_ADMIN") {
    return { error: new Response("Forbidden", { status: 403 }), session: null };
  }

  return { error: null, session };
}

export function isAdminRole(role: string | null | undefined): boolean {
  return role === "ADMIN" || role === "SUPER_ADMIN";
}

/**
 * After requireUser(), check the current user may edit a dataset:
 *   - the personal owner (ownerId === me), any role, OR
 *   - any admin when the dataset is admin-owned (shared).
 * A personal dataset is editable only by its owner — an admin's personal
 * dataset is NOT shared with other admins. Unclaimed datasets
 * (no owner, not admin-owned) are editable by nobody until claimed.
 * Returns a 403 Response to send back, or null if OK.
 */
export function ownerOrAdminError(
  dataset: { ownerId?: string | null; adminOwned?: boolean | null },
  session: { user?: { id?: string | null; role?: string | null } } | null,
): Response | null {
  const uid = session?.user?.id;

  // Personal owner edits their own dataset.
  if (dataset.ownerId && uid && dataset.ownerId === uid) return null;

  // Admin-owned (shared) datasets: any admin may edit.
  if (dataset.adminOwned && isAdminRole(session?.user?.role)) return null;

  return new Response("Forbidden", { status: 403 });
}

/**
 * For the calls an upload makes as it goes — marking files done, progress,
 * completing it, the "ready" email: the signed-in owner of the dataset (as
 * ownerOrAdminError), and, when an uploadId is given, an upload session that
 * belongs to this dataset. These routes used to trust whoever knew the ids.
 * Returns the session, or a Response to send back.
 */
export async function requireUploadOwner(
  datasetId: string,
  uploadId?: string | null,
) {
  const { error, session } = await requireUser();

  if (error) return { error, session: null };

  const dataset = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: { ownerId: true, adminOwned: true },
  });

  if (!dataset) {
    return { error: new Response("Not found", { status: 404 }), session: null };
  }
  const denied = ownerOrAdminError(dataset, session);

  if (denied) return { error: denied, session: null };

  if (uploadId) {
    const upload = await prisma.uploadSession.findUnique({
      where: { id: uploadId },
      select: { datasetId: true },
    });

    if (upload?.datasetId !== datasetId) {
      return {
        error: new Response("Upload not found", { status: 404 }),
        session: null,
      };
    }
  }

  return { error: null, session };
}
