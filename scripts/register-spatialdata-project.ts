/**
 * Register SpatialData stores on S3 as admin-owned datasets and put them, in
 * order, in one project — the shape the viewers' film strip reads.
 *
 * Each URL goes through the same ensureDatasetForS3Url() as "Add this
 * dataset", so the rows are what the app itself would create. Idempotent:
 * datasets are reused by URL, the project by title, and its member list is
 * rewritten to match.
 *
 * DRY RUN BY DEFAULT — pass --apply to write.
 *
 *   npx tsx scripts/register-spatialdata-project.ts --title "12_16_2025 sets" \
 *     --owner admin@example.com <store URL>... [--apply] [--use-direct-url]
 *
 * DATABASE_URL must point at the database you intend to change.
 */

if (process.argv.includes("--use-direct-url")) {
  if (!process.env.DIRECT_URL) {
    throw new Error("--use-direct-url given but DIRECT_URL is not set");
  }
  process.env.DATABASE_URL = process.env.DIRECT_URL;
}

const args = process.argv.slice(2);
const flag = (name: string) => {
  const i = args.indexOf(name);

  return i > -1 ? args[i + 1] : undefined;
};
const APPLY = args.includes("--apply");
const TITLE = flag("--title");
const OWNER_EMAIL = flag("--owner") ?? process.env.SUPER_ADMIN_EMAIL;
const DESCRIPTION = flag("--description") ?? null;
const URLS = args.filter(
  (a, i) => a.startsWith("http") && !["--title", "--owner", "--description"].includes(args[i - 1]),
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
  let where = "(DATABASE_URL unset)";

  try {
    const u = new URL(process.env.DATABASE_URL ?? "");

    where = `${u.hostname}:${u.port || "5432"}`;
  } catch {
    /* leave as unset */
  }
  log(`database : ${where}`);
  log(`mode     : ${APPLY ? "APPLY — this writes to the database" : "DRY RUN — nothing will be written"}`);
  if (!TITLE) throw new Error("Pass --title <project title>");
  if (!OWNER_EMAIL) throw new Error("Pass --owner <email> or set SUPER_ADMIN_EMAIL");
  if (URLS.length === 0) throw new Error("Pass the stores' URLs, in order");

  const owner = await prisma.user.findUnique({
    where: { email: OWNER_EMAIL },
    select: { id: true, email: true, role: true },
  });

  if (!owner) throw new Error(`No user with email ${OWNER_EMAIL}`);
  if (owner.role !== "ADMIN" && owner.role !== "SUPER_ADMIN") {
    throw new Error(`${owner.email} has role ${owner.role}; needs ADMIN or SUPER_ADMIN`);
  }
  log(`owner    : ${owner.email} (${owner.role})`);
  log(`project  : "${TITLE}", ${URLS.length} datasets\n`);

  const ids: string[] = [];

  for (const url of URLS) {
    if (!APPLY) {
      log(`  would register  ${url}`);
      continue;
    }
    const { dataset, claimed } = await ensureDatasetForS3Url(url, {
      userId: owner.id,
      isAdmin: true,
      asAdmin: true,
      kind: "spatialdata",
    });

    ids.push(dataset.id);
    log(`  ${claimed ? "registered" : "existing  "}  ${dataset.id}  ${dataset.title}  ` +
        `${dataset.datasetType}  ${dataset.numCells} cells`);
  }
  if (!APPLY) return;

  const project =
    (await prisma.project.findFirst({ where: { title: TITLE, ownerId: owner.id } })) ??
    (await prisma.project.create({
      data: { title: TITLE, description: DESCRIPTION, ownerId: owner.id },
    }));

  await prisma.projectDataset.deleteMany({ where: { projectId: project.id } });
  await prisma.projectDataset.createMany({
    data: ids.map((datasetId, sortOrder) => ({ projectId: project.id, datasetId, sortOrder })),
  });
  log(`\nproject ${project.id}: ${ids.length} datasets in order`);
}

main()
  .then(() => process.exit(0))
  .catch((e) => {
    // eslint-disable-next-line no-console
    console.error(e);
    process.exit(1);
  });
