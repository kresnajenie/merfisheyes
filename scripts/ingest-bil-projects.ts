/**
 * Organize the BIL catalog datasets into Projects, fixing the messy
 * auto-registered Dataset titles along the way.
 *
 * Background: viewing a BIL dataset via /viewer/from-s3 auto-registers a
 * `Dataset` row titled from the S3 manifest's folder name — so the account
 * pickers fill up with dozens of rows all named "combined_output" or
 * "cellpose-detected_transcripts". This script walks a catalog-import JSON
 * (scripts/new-datasets/bil-completed.json) and, per BIL item:
 *
 *   1. find-or-creates the Dataset row for each entry's s3BaseUrl (via the
 *      same ensureDatasetForS3Url the app uses; rows become admin-owned)
 *   2. renames it to "{BIL title} — {entry label}" and fills empty metadata
 *      (description/species/tissue/platform/thumbnail) from the BIL item
 *   3. creates one Project per BIL item (owner = --owner account) carrying
 *      the BIL metadata, and attaches the datasets in entry order
 *
 * Idempotent: projects are matched by metadata.bilCode, memberships upserted.
 * Rows personally owned by someone else are attached but never renamed.
 *
 *   npx tsx scripts/ingest-bil-projects.ts                     # dry run
 *   npx tsx scripts/ingest-bil-projects.ts --apply             # write
 *   npx tsx scripts/ingest-bil-projects.ts --owner=me@x.com    # other owner
 *
 * Requires DATABASE_URL (+ AWS creds only for non-public buckets) in the env.
 */
import { readFileSync } from "fs";

import { prisma } from "../lib/prisma";
import { ensureDatasetForS3Url } from "../lib/datasets/register-or-claim";
import { updateDatasetGenes } from "../lib/datasets/read-genes";
import { normalizeS3Url } from "../lib/utils/viewer-config";

interface BilEntry {
  label: string;
  datasetType: string;
  s3BaseUrl: string;
  thumbnailUrl?: string | null;
}

interface BilItem {
  title: string;
  description?: string | null;
  species?: string | null;
  tissue?: string | null;
  platform?: string | null;
  thumbnailUrl?: string | null;
  bilCode: string;
  externalLink?: string | null;
  publicationLink?: string | null;
  metadata?: Record<string, unknown>;
  entries: BilEntry[];
}

const apply = process.argv.includes("--apply");
const ownerEmail =
  process.argv.find((a) => a.startsWith("--owner="))?.slice("--owner=".length) ??
  "kresna.jenie@gmail.com";
const file =
  process.argv.find((a) => a.startsWith("--file="))?.slice("--file=".length) ??
  "scripts/new-datasets/bil-completed.json";

const cap = (s: string, n: number) => (s.length > n ? s.slice(0, n) : s);

async function main() {
  const items = JSON.parse(readFileSync(file, "utf-8")) as BilItem[];
  const owner = await prisma.user.findUnique({
    where: { email: ownerEmail },
    select: { id: true, email: true, role: true },
  });

  if (!owner) throw new Error(`No user with email ${ownerEmail}`);
  console.log(
    `${items.length} BIL item(s) from ${file} → projects owned by ${owner.email}` +
      `${apply ? "" : "  (dry run — pass --apply to write)"}\n`,
  );

  let projectsCreated = 0;
  let projectsUpdated = 0;
  let renamed = 0;
  let attached = 0;

  for (const item of items) {
    if (!item.bilCode || !Array.isArray(item.entries)) {
      console.log(`  !! skipping malformed item: ${item.title?.slice(0, 50)}`);
      continue;
    }

    const datasetIds: string[] = [];
    const lines: string[] = [];

    for (const entry of item.entries) {
      if (!entry.s3BaseUrl) continue;
      const url = normalizeS3Url(entry.s3BaseUrl);
      const newTitle = cap(`${item.title} — ${entry.label}`, 500);

      if (!apply) {
        const existing = await prisma.dataset.findUnique({
          where: { s3BaseUrl: url },
          select: { id: true, title: true, ownerId: true, adminOwned: true },
        });

        lines.push(
          existing
            ? `    ${existing.id} "${existing.title?.slice(0, 40)}" → "${newTitle.slice(0, 60)}"`
            : `    (new registration) → "${newTitle.slice(0, 60)}"`,
        );
        continue;
      }

      // Register/claim exactly like the app would, as collective admin rows.
      const { dataset, ownedByOther } = await ensureDatasetForS3Url(url, {
        userId: owner.id,
        isAdmin: true,
        asAdmin: true,
      });

      datasetIds.push(dataset.id);

      if (ownedByOther) {
        lines.push(
          `    ${dataset.id} owned by another user — attached, not renamed`,
        );
        continue;
      }

      // Rename + fill empty display fields from the BIL item (never clobber
      // values someone already set).
      await prisma.dataset.update({
        where: { id: dataset.id },
        data: {
          title: newTitle,
          description: dataset.description ?? item.description ?? undefined,
          species: dataset.species ?? item.species ?? undefined,
          tissue: dataset.tissue ?? item.tissue ?? undefined,
          platform: dataset.platform ?? item.platform ?? undefined,
          externalLink: dataset.externalLink ?? item.externalLink ?? undefined,
          thumbnailUrl:
            dataset.thumbnailUrl ??
            entry.thumbnailUrl ??
            item.thumbnailUrl ??
            undefined,
          metadata:
            dataset.metadata && Object.keys(dataset.metadata as object).length
              ? undefined
              : ((item.metadata ?? {}) as object),
        },
      });
      renamed++;

      // Newly created rows have no genes yet — index them (best-effort).
      if ((dataset.genes ?? []).length === 0) {
        await updateDatasetGenes(prisma, {
          id: dataset.id,
          datasetType: dataset.datasetType,
          formatVersion: dataset.formatVersion,
          manifestJson: dataset.manifestJson,
          ingestSource: dataset.ingestSource,
          s3BaseUrl: dataset.s3BaseUrl,
        });
      }
    }

    if (!apply) {
      console.log(`  ▸ ${item.bilCode}  "${item.title.slice(0, 60)}" (${item.entries.length} entries)`);
      lines.forEach((l) => console.log(l));
      continue;
    }

    // One project per BIL item, matched by metadata.bilCode (idempotent).
    const projectData = {
      title: cap(item.title, 500),
      description: item.description ?? null,
      species: item.species ?? null,
      tissue: item.tissue ?? null,
      platform: item.platform ?? null,
      thumbnailUrl: item.thumbnailUrl ?? null,
      externalLink: item.externalLink ?? null,
      publicationLink: item.publicationLink ?? null,
      metadata: { ...(item.metadata ?? {}), bilCode: item.bilCode } as object,
    };

    const existingProject = await prisma.project.findFirst({
      where: {
        ownerId: owner.id,
        metadata: { path: ["bilCode"], equals: item.bilCode },
      },
      select: { id: true },
    });

    const project = existingProject
      ? await prisma.project.update({
          where: { id: existingProject.id },
          data: projectData,
        })
      : await prisma.project.create({
          data: { ...projectData, ownerId: owner.id },
        });

    if (existingProject) projectsUpdated++;
    else projectsCreated++;

    for (let i = 0; i < datasetIds.length; i++) {
      await prisma.projectDataset.upsert({
        where: {
          projectId_datasetId: { projectId: project.id, datasetId: datasetIds[i] },
        },
        update: { sortOrder: i },
        create: { projectId: project.id, datasetId: datasetIds[i], sortOrder: i },
      });
      attached++;
    }

    console.log(
      `  ✓ ${item.bilCode}  "${item.title.slice(0, 60)}" — ${datasetIds.length} dataset(s)`,
    );
  }

  console.log(
    `\nDone: ${projectsCreated} project(s) created, ${projectsUpdated} updated, ` +
      `${renamed} dataset(s) renamed, ${attached} membership(s) upserted.`,
  );
}

main()
  .catch((e) => {
    console.error(e);
    process.exitCode = 1;
  })
  .finally(() => prisma.$disconnect());
