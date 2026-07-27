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

/**
 * After requireUser(), check the current user may edit `ownerId`'s dataset:
 * the owner themselves, or any ADMIN/SUPER_ADMIN. Ownerless (null) rows are
 * editable by admins only. Returns a 403 Response to send back, or null if OK.
 */
export function ownerOrAdminError(
  ownerId: string | null | undefined,
  session: { user?: { id?: string | null; role?: string | null } } | null,
): Response | null {
  const role = session?.user?.role;

  if (role === "ADMIN" || role === "SUPER_ADMIN") return null;

  const uid = session?.user?.id;

  if (ownerId && uid && ownerId === uid) return null;

  return new Response("Forbidden", { status: 403 });
}
