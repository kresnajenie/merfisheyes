import type { Prisma } from "@prisma/client";

/**
 * Resolve a dataset's owner tuple from the actor + their choice.
 *
 * Two owned states, never empty (callers must pass a real `userId`, i.e. the
 * flow requires sign-in):
 *   - admin choosing "admin" → collective: { ownerId: null, adminOwned: true }
 *   - otherwise             → personal:   { ownerId: userId, adminOwned: false }
 *
 * `asAdmin` is only honored for actual admins, mirroring register-s3.
 */
export function resolveDatasetOwner({
  isAdmin,
  asAdmin,
  userId,
}: {
  isAdmin: boolean;
  asAdmin: boolean;
  userId: string;
}): { ownerId: string | null; adminOwned: boolean } {
  const adminClaim = asAdmin && isAdmin;

  return adminClaim
    ? { ownerId: null, adminOwned: true }
    : { ownerId: userId, adminOwned: false };
}

/**
 * The datasets that are a person's own: their personal ones, plus the
 * admin-owned (shared) ones when they are an admin. "My datasets" on the
 * account page and in the agent API (lib/agent/visibility.ts).
 */
export function ownedDatasetsWhere(
  userId: string,
  isAdmin: boolean,
): Prisma.DatasetWhereInput {
  return isAdmin
    ? { OR: [{ ownerId: userId }, { adminOwned: true }] }
    : { ownerId: userId };
}
