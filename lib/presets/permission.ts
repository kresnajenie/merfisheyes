import { auth } from "@/lib/auth";
import { isAdminRole } from "@/lib/admin-auth";
import { prisma } from "@/lib/prisma";

/** The permission to save and publish technology presets, and nothing else. */
export const PRESETS_WRITE = "presets:write";

/** Whether this user may edit presets: an admin, or granted the permission. */
export async function canEditPresets(
  user: { id?: string | null; role?: string | null } | null | undefined,
): Promise<boolean> {
  if (!user?.id) return false;
  if (isAdminRole(user.role)) return true;

  const row = await prisma.user.findUnique({
    where: { id: user.id },
    select: { permissions: true },
  });

  return !!row?.permissions.includes(PRESETS_WRITE);
}

/**
 * Verify the current request is from a user who may edit presets.
 * Returns the session if valid, or a Response to send back.
 */
export async function requirePresetEditor() {
  const session = await auth();

  if (!session?.user) {
    return {
      error: new Response("Unauthorized", { status: 401 }),
      session: null,
    };
  }
  if (!(await canEditPresets(session.user))) {
    return { error: new Response("Forbidden", { status: 403 }), session: null };
  }

  return { error: null, session };
}
