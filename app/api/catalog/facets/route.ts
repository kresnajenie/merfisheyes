import { NextResponse } from "next/server";

import { auth } from "@/lib/auth";
import { catalogFacets } from "@/lib/explore/search";
import { isAdminRole } from "@/lib/admin-auth";

/**
 * Distinct species / tissue / platform values already used across the published
 * catalog. Powers the "pick an existing value before adding a new one"
 * autocompletes in the account editor, so people converge on canonical terms
 * (no "Mouse" vs "Mice" split) instead of inventing a new spelling each time.
 * The values of the rows Explore shows the caller (lib/explore/visibility.ts):
 * internal rows' values only for admins.
 */
export async function GET() {
  const session = await auth();

  return NextResponse.json(
    await catalogFacets(isAdminRole(session?.user?.role)),
  );
}
