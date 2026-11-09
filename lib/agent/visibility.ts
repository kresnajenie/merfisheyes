import type { Prisma } from "@prisma/client";
import type { Caller } from "@/lib/agent/caller";

import { prisma } from "@/lib/prisma";
import { catalogRowWhere } from "@/lib/explore/visibility";
import { ownedDatasetsWhere } from "@/lib/datasets/owner";

/**
 * What an agent may see (docs/specs/MCP-PROTOCOL.md §2). Catalog rows follow
 * Explore's rule (lib/explore/visibility.ts). Datasets have no visibility
 * flag on the site — any id opens — so the agent has its own rule, here and
 * only here: a dataset is visible when a visible catalog row names it, or
 * when it is the caller's own and inside the grant's reach.
 */

type User = Extract<Caller, { kind: "user" }>;

/** The caller's projects within the grant's reach. Admins get no extra
 * reach: an agent sees the person's own projects only. */
export function projectsWhere(caller: User): Prisma.ProjectWhereInput {
  return caller.projects === "all"
    ? { ownerId: caller.userId }
    : { ownerId: caller.userId, id: { in: caller.projects } };
}

/** The caller's own datasets within the grant's reach. */
export function ownDatasetsWhere(caller: User): Prisma.DatasetWhereInput {
  const own = ownedDatasetsWhere(caller.userId, caller.isAdmin);

  return caller.projects === "all"
    ? own
    : {
        AND: [own, { projects: { some: { project: projectsWhere(caller) } } }],
      };
}

/** A catalog row by id or BIL code, if visible. */
export function findCatalogRow(caller: Caller, idOrCode: string) {
  const isAdmin = caller.kind === "user" && caller.isAdmin;

  return prisma.catalogDataset.findFirst({
    where: {
      AND: [
        { OR: [{ id: idOrCode }, { bilCode: idOrCode }] },
        catalogRowWhere(isAdmin),
      ],
    },
    include: { entries: { orderBy: { sortOrder: "asc" } } },
  });
}

/** A dataset by id, if visible to the caller. */
export async function findDataset(caller: Caller, id: string) {
  const dataset = await prisma.dataset.findUnique({ where: { id } });

  if (!dataset) return null;
  if (caller.kind === "user") {
    const own = await prisma.dataset.count({
      where: { AND: [{ id }, ownDatasetsWhere(caller)] },
    });

    if (own) return dataset;
  }
  const isAdmin = caller.kind === "user" && caller.isAdmin;
  const named = await prisma.catalogDataset.count({
    where: {
      AND: [
        catalogRowWhere(isAdmin),
        {
          entries: {
            some: {
              OR: [
                { datasetId: id },
                ...(dataset.s3BaseUrl
                  ? [{ s3BaseUrl: dataset.s3BaseUrl }]
                  : []),
              ],
            },
          },
        },
      ],
    },
  });

  return named ? dataset : null;
}
