/**
 * Backfill Dataset.genes for existing COMPLETE datasets, so pre-existing
 * uploads become gene-searchable (new uploads are indexed at completion; see
 * lib/datasets/read-genes.ts, which this script reuses).
 *
 * Reads each dataset's genes from what already exists — the S3 expr index
 * (single cell) or the manifest (single molecule) — and writes them to the
 * row. Skips rows that already have genes. Safe to re-run.
 *
 *   npx tsx scripts/backfill-dataset-genes.ts            # dry run
 *   npx tsx scripts/backfill-dataset-genes.ts --apply    # write to the DB
 *
 * Requires DATABASE_URL + AWS creds in the environment (e.g. run with
 * `node --env-file=.env` semantics via: npx dotenv -e .env -- npx tsx ...`
 * or export them first).
 */
import { prisma } from "../lib/prisma";
import { readDatasetGenes } from "../lib/datasets/read-genes";

const apply = process.argv.includes("--apply");

async function main() {
  const datasets = await prisma.dataset.findMany({
    where: { status: "COMPLETE", genes: { isEmpty: true } },
    select: {
      id: true,
      title: true,
      datasetType: true,
      formatVersion: true,
      manifestJson: true,
      ingestSource: true,
      s3BaseUrl: true,
    },
    orderBy: { createdAt: "asc" },
  });

  console.log(
    `${datasets.length} COMPLETE dataset(s) without genes${apply ? "" : " (dry run — pass --apply to write)"}`,
  );

  let filled = 0;
  let skipped = 0;

  for (const d of datasets) {
    const genes = await readDatasetGenes(d);

    if (!genes) {
      skipped++;
      console.log(
        `  – ${d.id} (${d.datasetType ?? "?"} / ${d.formatVersion ?? "?"}) "${(d.title ?? "").slice(0, 40)}": no genes found (zarr / missing objects)`,
      );
      continue;
    }

    if (apply) {
      await prisma.dataset.update({ where: { id: d.id }, data: { genes } });
    }
    filled++;
    console.log(
      `  ✓ ${d.id} "${(d.title ?? "").slice(0, 40)}": ${genes.length} genes${apply ? "" : " (would write)"}`,
    );
  }

  console.log(`\nDone: ${filled} filled, ${skipped} skipped.`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exitCode = 1;
  })
  .finally(() => prisma.$disconnect());
