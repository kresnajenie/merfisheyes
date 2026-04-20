/**
 * Seed script for BIL (Brain Image Library) datasets.
 * Reads from bil-examples.json and upserts into CatalogDataset + CatalogDatasetEntry.
 *
 * Run with: npx tsx prisma/seed-bil-data.ts
 *
 * Idempotent: uses bilCode as unique key. Re-running updates existing records.
 */
import { PrismaClient, Prisma } from "@prisma/client";
import * as fs from "fs";
import * as path from "path";

const prisma = new PrismaClient();

interface SeedEntry {
  label: string;
  datasetType: "single_cell" | "single_molecule";
  s3BaseUrl?: string;
  datasetId?: string;
  thumbnailUrl?: string;
}

interface SeedDataset {
  title: string;
  description: string;
  species: string;
  tissue: string;
  platform: string;
  bilCode: string;
  externalLink: string;
  genes: string[];
  metadata: Record<string, unknown>;
  entries: SeedEntry[];
}

async function main() {
  const jsonPath = path.join(__dirname, "bil-examples.json");
  const raw = fs.readFileSync(jsonPath, "utf-8");
  const datasets: SeedDataset[] = JSON.parse(raw);

  console.log(`Seeding ${datasets.length} BIL datasets from bil-examples.json...`);

  let created = 0;
  let updated = 0;
  let skipped = 0;

  for (const ds of datasets) {
    if (!ds.bilCode) {
      console.log(`  SKIP — no bilCode: ${ds.title.slice(0, 60)}`);
      skipped++;
      continue;
    }

    // Upsert by bilCode
    const existing = await prisma.catalogDataset.findUnique({
      where: { bilCode: ds.bilCode },
      include: { entries: true },
    });

    // Filter out invalid entries (e.g., placeholder "..." strings)
    const validEntries = ds.entries.filter(
      (e): e is SeedEntry => typeof e === "object" && e !== null && !!e.label && !!e.datasetType,
    );

    const catalogData: Prisma.CatalogDatasetUpdateInput = {
      title: ds.title,
      description: ds.description,
      species: ds.species,
      tissue: ds.tissue,
      platform: ds.platform,
      externalLink: ds.externalLink,
      genes: ds.genes,
      metadata: ds.metadata as Prisma.InputJsonValue,
      numGenes: ds.genes.length,
      isBil: true,
      isPublished: true,
    };

    if (existing) {
      // Delete old entries and recreate
      await prisma.catalogDatasetEntry.deleteMany({
        where: { catalogId: existing.id },
      });

      await prisma.catalogDataset.update({
        where: { id: existing.id },
        data: {
          ...catalogData,
          entries: {
            create: validEntries.map((e, i) => ({
              label: e.label,
              datasetType: e.datasetType,
              s3BaseUrl: e.s3BaseUrl ?? null,
              datasetId: e.datasetId ?? null,
              thumbnailUrl: e.thumbnailUrl ?? null,
              sortOrder: i,
            })),
          },
        },
      });

      updated++;
      process.stdout.write(`  Updated: ${ds.bilCode} (${validEntries.length} entries)\n`);
    } else {
      await prisma.catalogDataset.create({
        data: {
          title: ds.title,
          description: ds.description,
          species: ds.species,
          tissue: ds.tissue,
          platform: ds.platform,
          externalLink: ds.externalLink,
          bilCode: ds.bilCode,
          genes: ds.genes,
          metadata: ds.metadata as Prisma.InputJsonValue,
          numGenes: ds.genes.length,
          isBil: true,
          isPublished: true,
          entries: {
            create: validEntries.map((e, i) => ({
              label: e.label,
              datasetType: e.datasetType,
              s3BaseUrl: e.s3BaseUrl ?? null,
              datasetId: e.datasetId ?? null,
              thumbnailUrl: e.thumbnailUrl ?? null,
              sortOrder: i,
            })),
          },
        },
      });

      created++;
      process.stdout.write(`  Created: ${ds.bilCode} (${validEntries.length} entries)\n`);
    }
  }

  const total = await prisma.catalogDataset.count({ where: { isBil: true } });

  console.log(`\nDone. Created: ${created}, Updated: ${updated}, Skipped: ${skipped}`);
  console.log(`Total BIL datasets in database: ${total}`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exit(1);
  })
  .finally(() => prisma.$disconnect());
