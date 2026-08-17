import { NextResponse } from "next/server";

import { prisma } from "@/lib/prisma";

/**
 * Distinct species / tissue / platform values already used across the published
 * catalog. Powers the "pick an existing value before adding a new one"
 * autocompletes in the account editor, so people converge on canonical terms
 * (no "Mouse" vs "Mice" split) instead of inventing a new spelling each time.
 * Public — these are already visible on Explore.
 */
export async function GET() {
  const where = { isPublished: true } as const;

  const [speciesRaw, tissueRaw, platformRaw] = await Promise.all([
    prisma.catalogDataset.findMany({
      where: { ...where, species: { not: null } },
      select: { species: true },
      distinct: ["species"],
      orderBy: { species: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...where, tissue: { not: null } },
      select: { tissue: true },
      distinct: ["tissue"],
      orderBy: { tissue: "asc" },
    }),
    prisma.catalogDataset.findMany({
      where: { ...where, platform: { not: null } },
      select: { platform: true },
      distinct: ["platform"],
      orderBy: { platform: "asc" },
    }),
  ]);

  return NextResponse.json({
    species: speciesRaw.map((r) => r.species).filter(Boolean),
    tissues: tissueRaw.map((r) => r.tissue).filter(Boolean),
    platforms: platformRaw.map((r) => r.platform).filter(Boolean),
  });
}
