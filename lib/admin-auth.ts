import { auth } from "@/lib/auth";

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

function isAdminRole(role: string | null | undefined): boolean {
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
