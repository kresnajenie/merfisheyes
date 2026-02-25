import { auth } from "@/lib/auth";

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
