/**
 * Register the built spiralia datasets and group them into one Project.
 *
 *   npx tsx scripts/spiralia/register_datasets.ts --dry-run
 *   npx tsx scripts/spiralia/register_datasets.ts
 *
 * Reads scripts/spiralia/build_report.csv (what was actually built) and the
 * embryo metadata for stages, then upserts one admin-owned Dataset per embryo
 * keyed on its S3 URL, and a Project holding them ordered by developmental
 * stage — that order is what a stage filmstrip scrubs through.
 *
 * Idempotent: `s3BaseUrl` is unique, so re-running updates rather than
 * duplicates. Safe to re-run after rebuilding a dataset.
 */

import { readFileSync } from "node:fs";

import { PrismaClient } from "@prisma/client";

const prisma = new PrismaClient();

const S3 = "https://merfisheyes-bil.s3.us-west-2.amazonaws.com/yiqun-spiralia";
const PROJECT_TITLE = "Spiralia embryo atlas";
const META = "/home/data/yiqun-spiralia/Sep2026/embryo_metadata_deDup_May12_2026.csv";
const REPORT = "/home/kjenie/merfisheyes/scripts/spiralia/build_report.csv";

/** Developmental order; anything unlisted sorts after, alphabetically. */
const STAGE_ORDER = [
  "1-cell", "2-cells", "4-cells", "4-8-cells", "8-cells", "12-cells",
  "16-cells", "16-20-cells", "20-24-cells", "24-cells", "Organizer",
];

function csv(path: string): Record<string, string>[] {
  const [head, ...rows] = readFileSync(path, "utf8").trim().split("\n");
  const cols = head.split(",");

  return rows.map((r) => {
    const cells = r.split(",");

    return Object.fromEntries(cols.map((c, i) => [c, cells[i]]));
  });
}

async function main() {
  const dryRun = process.argv.includes("--dry-run");
  const meta = new Map(csv(META).map((r) => [r.embryo, r]));
  const built = csv(REPORT).filter((r) => r.status === "ok");

  const admin = await prisma.user.findFirst({
    where: { role: "ADMIN" },
    select: { id: true, email: true },
  });

  if (!admin) throw new Error("no ADMIN user to own these datasets");

  const rows = built
    .map((r) => {
      const m = meta.get(r.embryo);
      const stage = m?.stage ?? "";
      const rank = STAGE_ORDER.indexOf(stage);

      return {
        embryo: r.embryo,
        stage,
        // Unlisted stages sort last but stay grouped.
        rank: rank === -1 ? STAGE_ORDER.length : rank,
        molecules: Number(r.molecules),
        genes: Number(r.genes),
        cells: Number(r.cells),
      };
    })
    .sort((a, b) => a.rank - b.rank || a.embryo.localeCompare(b.embryo));

  console.log(`${rows.length} datasets, owner ${admin.email} (adminOwned)`);
  if (dryRun) {
    for (const r of rows.slice(0, 5)) {
      console.log(`  ${r.embryo.padEnd(15)} ${r.stage.padEnd(13)} ${r.molecules.toLocaleString()} molecules`);
    }
    console.log(`  … and ${rows.length - 5} more`);
    console.log("dry run — nothing written");

    return;
  }

  const ids: string[] = [];

  for (const r of rows) {
    const s3BaseUrl = `${S3}/${r.embryo}_lm`;
    // Deterministic id so a re-run targets the same row even if the unique
    // s3BaseUrl lookup is ever changed.
    const id = `lm-spiralia-${r.embryo.toLowerCase().replace(/[^a-z0-9]+/g, "-")}`;
    const common = {
      title: r.embryo,
      numCells: r.molecules,
      numGenes: r.genes,
      datasetType: "labelled_single_molecule",
      status: "COMPLETE" as const,
      s3BaseUrl,
      adminOwned: true,
      ownerId: admin.id,
      species: "Spiralia",
      description: `${r.stage} — ${r.molecules.toLocaleString()} molecules, ${r.genes} genes, ${r.cells} cells.`,
      tags: [r.stage].filter(Boolean),
      metadata: { stage: r.stage, embryo: r.embryo },
    };

    const existing = await prisma.dataset.findUnique({
      where: { s3BaseUrl },
      select: { id: true },
    });

    if (existing) {
      await prisma.dataset.update({ where: { id: existing.id }, data: common });
      ids.push(existing.id);
    } else {
      await prisma.dataset.create({
        data: { id, fingerprint: `lm-spiralia-${r.embryo}`, ...common },
      });
      ids.push(id);
    }
  }

  let project = await prisma.project.findFirst({
    where: { title: PROJECT_TITLE, ownerId: admin.id },
    select: { id: true },
  });

  if (!project) {
    project = await prisma.project.create({
      data: {
        title: PROJECT_TITLE,
        ownerId: admin.id,
        species: "Spiralia",
        description:
          "MERFISH single molecules with per-molecule gene, RNA domain and cell " +
          "labels, across development from the 1-cell stage to 24 cells.",
      },
      select: { id: true },
    });
  }

  // sortOrder is the stage order — what the filmstrip scrubs through.
  for (const [i, datasetId] of ids.entries()) {
    await prisma.projectDataset.upsert({
      where: { projectId_datasetId: { projectId: project.id, datasetId } },
      create: { projectId: project.id, datasetId, sortOrder: i },
      update: { sortOrder: i },
    });
  }

  console.log(`registered ${ids.length} datasets`);
  console.log(`project ${project.id} — ${PROJECT_TITLE}`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exit(1);
  })
  .finally(() => prisma.$disconnect());
