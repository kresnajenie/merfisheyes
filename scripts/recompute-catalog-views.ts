/**
 * Recompute CatalogDataset.viewCount from the per-dataset counts.
 *
 * The column is normally kept current by /api/datasets/[datasetId]/view, which
 * bumps it alongside Dataset.viewCount. It drifts when a catalog row's entries
 * change after the fact — the admin editor deletes and recreates entries, so
 * adding an entry for a dataset that already had views leaves the catalog
 * total behind. Run this to re-derive the whole column. Safe to re-run.
 *
 *   npx tsx scripts/recompute-catalog-views.ts            # dry run
 *   npx tsx scripts/recompute-catalog-views.ts --apply    # write to the DB
 *
 * Requires DATABASE_URL in the environment.
 */
import { prisma } from "../lib/prisma";

const apply = process.argv.includes("--apply");

/** Sum of Dataset.viewCount over every dataset a row's entries resolve to —
 * by datasetId, or by s3BaseUrl for curated entries. EXISTS rather than a join
 * so a dataset reachable through several entries counts once, matching how the
 * live increment behaves. */
const EXPECTED_SQL = `
  SELECT c."id",
         c."title",
         c."view_count" AS current,
         COALESCE((
           SELECT SUM(d."view_count")
             FROM "datasets" d
            WHERE EXISTS (
                  SELECT 1
                    FROM "catalog_dataset_entries" e
                   WHERE e."catalog_id" = c."id"
                     AND (e."dataset_id" = d."id" OR e."s3_base_url" = d."s3_base_url")
            )
         ), 0)::int AS expected
    FROM "catalog_datasets" c
   ORDER BY expected DESC`;

async function main() {
  const rows =
    await prisma.$queryRawUnsafe<
      { id: string; title: string; current: number; expected: number }[]
    >(EXPECTED_SQL);
  const drifted = rows.filter((r) => r.current !== r.expected);

  console.log(
    `${rows.length} catalog row(s), ${drifted.length} with a stale count${
      apply ? "" : " (dry run — pass --apply to write)"
    }`,
  );

  for (const r of drifted) {
    if (apply) {
      await prisma.catalogDataset.update({
        where: { id: r.id },
        data: { viewCount: r.expected },
      });
    }
    console.log(
      `  ${apply ? "✓" : "–"} ${r.id} "${r.title.slice(0, 40)}": ${r.current} → ${r.expected}`,
    );
  }

  console.log(`\nDone: ${apply ? drifted.length : 0} updated.`);
}

main()
  .catch((e) => {
    console.error(e);
    process.exitCode = 1;
  })
  .finally(() => prisma.$disconnect());
