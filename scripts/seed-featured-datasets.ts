/**
 * Seed the featured-datasets S3 prefix into the catalog.
 *
 * For each atlas in scripts/featured-datasets.json this:
 *   1. registers every S3 URL as a Dataset row (admin-owned) via the same
 *      ensureDatasetForS3Url() the "Add this dataset" button uses,
 *   2. writes each dataset's viewerConfig so 3D stacks open with the 3D camera
 *      instead of the viewer's hard-coded 2D default,
 *   3. puts the atlas's datasets in a Project,
 *   4. creates one CatalogDataset per atlas with its datasets as entries, and
 *      publishes it (featuring the ones flagged in the manifest).
 *
 * Idempotent: re-running reuses rows by their natural keys (Dataset by
 * s3BaseUrl, Project/CatalogDataset by title) and rewrites the child rows.
 *
 * DRY RUN BY DEFAULT — it prints the plan and writes nothing until you pass
 * --apply.
 *
 *   npx tsx scripts/seed-featured-datasets.ts --owner you@example.com
 *   npx tsx scripts/seed-featured-datasets.ts --owner you@example.com --apply
 *
 * DATABASE_URL must point at the database you intend to change.
 */

import { readFileSync } from "fs";
import { join } from "path";

// Supabase's direct host (db.<ref>.supabase.co) is IPv6-only, so DATABASE_URL
// is unreachable from an IPv4-only machine even though it works on Vercel.
// --use-direct-url runs against DIRECT_URL (the IPv4 pooler) instead. The swap
// must happen before @/lib/prisma is loaded, hence the dynamic imports below.
if (process.argv.includes("--use-direct-url")) {
  if (!process.env.DIRECT_URL) {
    throw new Error("--use-direct-url given but DIRECT_URL is not set");
  }
  process.env.DATABASE_URL = process.env.DIRECT_URL;
}

type Entry = {
  slug: string;
  s3BaseUrl: string;
  label: string;
  sortOrder: number;
  datasetType: string;
  viewerConfig: Record<string, unknown>;
};

type Atlas = {
  key: string;
  featured: boolean;
  project: string;
  title: string;
  description: string;
  species: string | null;
  disease: string | null;
  tissue: string | null;
  platform: string | null;
  institute: string | null;
  externalLink: string | null;
  publicationLink: string | null;
  tags: string[];
  numCells: number;
  numGenes: number;
  datasets: Entry[];
};

const APPLY = process.argv.includes("--apply");
const ownerFlag = process.argv.indexOf("--owner");
const OWNER_EMAIL =
  ownerFlag > -1 ? process.argv[ownerFlag + 1] : process.env.SUPER_ADMIN_EMAIL;

const atlases: Atlas[] = JSON.parse(
  readFileSync(join(__dirname, "featured-datasets.json"), "utf8"),
);

function log(...a: unknown[]) {
  // eslint-disable-next-line no-console
  console.log(...a);
}

async function main() {
  const { prisma } = await import("@/lib/prisma");
  const { ensureDatasetForS3Url } = await import(
    "@/lib/datasets/register-or-claim"
  );

  // Identify the target by host only — never echo the connection string.
  let where = "(DATABASE_URL unset)";

  try {
    const u = new URL(process.env.DATABASE_URL ?? "");

    where = `${u.hostname}:${u.port || "5432"}`;
  } catch {
    /* leave as unset */
  }

  log(`\ndatabase : ${where}`);
  log(`mode     : ${APPLY ? "APPLY — this writes to the database" : "DRY RUN — nothing will be written"}`);

  if (!OWNER_EMAIL) {
    throw new Error("Pass --owner <email> or set SUPER_ADMIN_EMAIL — projects need an owner.");
  }

  const owner = await prisma.user.findUnique({
    where: { email: OWNER_EMAIL },
    select: { id: true, email: true, role: true },
  });

  if (!owner) throw new Error(`No user with email ${OWNER_EMAIL}`);
  if (owner.role !== "ADMIN" && owner.role !== "SUPER_ADMIN") {
    throw new Error(`${owner.email} has role ${owner.role}; needs ADMIN or SUPER_ADMIN`);
  }
  log(`owner    : ${owner.email} (${owner.role})`);

  const total = atlases.reduce((n, a) => n + a.datasets.length, 0);

  log(`plan     : ${atlases.length} atlases, ${total} datasets, ` +
      `${atlases.filter((a) => a.featured).length} featured\n`);

  for (const a of atlases) {
    const n3d = a.datasets.filter((d) => d.viewerConfig.viewMode === "3D").length;

    log(`── ${a.title}`);
    log(`   ${a.datasets.length} datasets (${n3d} open in 3D) · ` +
        `featured=${a.featured} · project "${a.project}"`);

    // 1. datasets ------------------------------------------------------
    const ids: { entry: Entry; datasetId: string }[] = [];

    for (const e of a.datasets) {
      if (!APPLY) {
        log(`   would register  ${e.slug}  viewMode=${e.viewerConfig.viewMode}`);
        continue;
      }
      const { dataset, claimed } = await ensureDatasetForS3Url(e.s3BaseUrl, {
        userId: owner.id,
        isAdmin: true,
        asAdmin: true,
      });

      await prisma.dataset.update({
        where: { id: dataset.id },
        data: { title: e.label, viewerConfig: e.viewerConfig },
      });
      ids.push({ entry: e, datasetId: dataset.id });
      log(`   ${claimed ? "registered" : "reused    "}  ${dataset.id}  ${e.slug}  viewMode=${e.viewerConfig.viewMode}`);
    }

    if (!APPLY) {
      log(`   would create project + catalog card, publish${a.featured ? " + feature" : ""}\n`);
      continue;
    }

    // 2. project -------------------------------------------------------
    const existingProject = await prisma.project.findFirst({
      where: { title: a.project, ownerId: owner.id },
      select: { id: true },
    });
    const project =
      existingProject ??
      (await prisma.project.create({
        data: {
          title: a.project,
          description: a.description,
          species: a.species,
          disease: a.disease,
          tissue: a.tissue,
          platform: a.platform,
          institute: a.institute,
          tags: a.tags,
          externalLink: a.externalLink,
          publicationLink: a.publicationLink,
          ownerId: owner.id,
        },
      }));

    for (const { entry, datasetId } of ids) {
      await prisma.projectDataset.upsert({
        where: { projectId_datasetId: { projectId: project.id, datasetId } },
        update: { sortOrder: entry.sortOrder },
        create: { projectId: project.id, datasetId, sortOrder: entry.sortOrder },
      });
    }
    log(`   project ${project.id} · ${ids.length} datasets attached`);

    // 3. catalog card --------------------------------------------------
    const existingCatalog = await prisma.catalogDataset.findFirst({
      where: { title: a.title },
      select: { id: true },
    });
    const data = {
      title: a.title,
      description: a.description,
      species: a.species,
      disease: a.disease,
      tissue: a.tissue,
      platform: a.platform,
      institute: a.institute,
      tags: a.tags,
      externalLink: a.externalLink,
      publicationLink: a.publicationLink,
      numCells: a.numCells,
      numGenes: a.numGenes,
      isPublished: true,
      isFeatured: a.featured,
      isBil: false,
      isInternal: false,
      isCommunity: false,
      sourceProjectId: project.id,
      createdBy: owner.id,
    };
    const catalog = existingCatalog
      ? await prisma.catalogDataset.update({ where: { id: existingCatalog.id }, data })
      : await prisma.catalogDataset.create({ data });

    // Entries are derived, so rewrite them wholesale rather than diffing.
    await prisma.catalogDatasetEntry.deleteMany({ where: { catalogId: catalog.id } });
    await prisma.catalogDatasetEntry.createMany({
      data: ids.map(({ entry, datasetId }) => ({
        catalogId: catalog.id,
        label: entry.label,
        datasetType: entry.datasetType,
        s3BaseUrl: entry.s3BaseUrl,
        datasetId,
        sortOrder: entry.sortOrder,
      })),
    });
    log(`   catalog ${catalog.id} · published${a.featured ? " · FEATURED" : ""} · ${ids.length} entries\n`);
  }

  if (!APPLY) {
    log("Nothing was written. Re-run with --apply to commit.\n");

    return;
  }

  const featured = await prisma.catalogDataset.count({ where: { isFeatured: true } });

  log(`done. ${featured} featured catalog rows now exist (including any pre-existing).\n`);
}

main()
  .catch((e) => {
    // Connection strings can appear in Prisma errors — strip them.
    const safe = String(e?.message ?? e)
      .replace(/postgres(ql)?:\/\/\S+/gi, "<connection-string>")
      .replace(/[\w.%+-]+:[^@\s]+@/g, "<credentials>@");

    // eslint-disable-next-line no-console
    console.error("\nFAILED:", safe, "\n");
    process.exitCode = 1;
  })
  .finally(async () => {
    const { prisma } = await import("@/lib/prisma");

    await prisma.$disconnect();
  });
