import type { TechnologyPreset } from "@/lib/config/technology-presets";

import { prisma } from "@/lib/prisma";

export interface PublishedPreset {
  hash: string;
  config: TechnologyPreset;
}

/**
 * The preset in force per technology: its latest publish row. A technology
 * whose latest row has no version is turned off and maps to null (which is
 * different from never published: the built-in preset does not apply either).
 */
export async function publishedPresets(): Promise<
  Record<string, PublishedPreset | null>
> {
  const rows = await prisma.technologyPresetPublish.findMany({
    orderBy: { createdAt: "desc" },
    distinct: ["technology"],
    select: {
      technology: true,
      version: { select: { hash: true, config: true } },
    },
  });

  return Object.fromEntries(
    rows.map((r) => [
      r.technology,
      r.version && {
        hash: r.version.hash,
        config: r.version.config as TechnologyPreset,
      },
    ]),
  );
}
