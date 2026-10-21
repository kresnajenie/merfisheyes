/**
 * Sync the Explore catalog with a SpatialData-DB catalog release.
 *
 * SpatialData-DB publishes its catalog as `catalog.parquet` on a GitHub
 * release of theislab/spatialdata-db-client (near-nightly, `catalog-<date>`).
 * Each row is one store on their public bucket. This script turns every row
 * into one Explore card (source "sddb", one `spatialdata` entry opened by the
 * store's public URL), registers the store as an admin-owned Dataset (the
 * same path as the viewer's "Add this dataset" banner, so views are counted
 * and the banner never shows), upserts by the catalog uid, and unpublishes
 * rows whose uid left the catalog. Spec: docs/specs/sddb-catalog.md.
 *
 * Idempotent: a second run on the same release changes nothing. Run by hand,
 * on the dev database first, then production (there is no schedule):
 *
 *   npx tsx scripts/sync-sddb-catalog.ts --env=.env                 # dry run, newest release
 *   npx tsx scripts/sync-sddb-catalog.ts --env=.env --limit=30 --apply   # a sample, spread over technologies
 *   npx tsx scripts/sync-sddb-catalog.ts --env=.env --apply
 *   npx tsx scripts/sync-sddb-catalog.ts --env=.env.prod --direct --apply
 *   npx tsx scripts/sync-sddb-catalog.ts --env=.env --file=path/to/catalog.parquet --tag=local
 *
 * Options: --tag=catalog-YYYY-MM-DD (default: the newest release carrying the
 * asset), --owner=<email> (the admin whose account registers the stores;
 * default kresna.jenie@gmail.com), --no-register (cards only),
 * --concurrency=8 (store reads in flight while registering).
 * Needs DATABASE_URL in the environment, or --env=<file> (with --direct to
 * use that file's DIRECT_URL, the pooler port Supabase accepts migrations on).
 * Releases are cached under .cache/sddb/<tag>/ (tags are immutable).
 */
import { createWriteStream, existsSync, mkdirSync, readFileSync, writeFileSync } from "fs";
import { join } from "path";
import { Readable } from "stream";
import { pipeline } from "stream/promises";

import {
  CATALOG_ASSET,
  CATALOG_REPO,
  PINNED_RELEASE,
} from "../lib/sddb/catalog-schema";
import { readCatalogFile } from "../lib/sddb/catalog-file";
import {
  mapCatalogRow,
  sampleRows,
  staleUids,
  type MappedCatalogRow,
} from "../lib/sddb/catalog-row";

const arg = (name: string) =>
  process.argv.find((a) => a.startsWith(`--${name}=`))?.slice(name.length + 3);
const apply = process.argv.includes("--apply");
const register = !process.argv.includes("--no-register");
const limit = Number(arg("limit") ?? 0) || 0;
const ownerEmail = arg("owner") ?? "kresna.jenie@gmail.com";
const concurrency = Number(arg("concurrency") ?? 8) || 8;

function loadEnv(file: string, direct: boolean) {
  for (const line of readFileSync(file, "utf8").split("\n")) {
    const m = /^([A-Z_]+)=(.*)$/.exec(line);

    if (m) process.env[m[1]] = m[2].replace(/^"|"$/g, "");
  }
  if (direct && process.env.DIRECT_URL) {
    process.env.DATABASE_URL = process.env.DIRECT_URL;
  }
}

interface ReleaseAsset {
  tag: string;
  url: string;
  size: number;
}

/** The newest release that carries the catalog asset (assets land a few hours after the tag). */
async function findRelease(tag: string | undefined): Promise<ReleaseAsset> {
  const url = tag
    ? `https://api.github.com/repos/${CATALOG_REPO}/releases/tags/${tag}`
    : `https://api.github.com/repos/${CATALOG_REPO}/releases?per_page=10`;
  const res = await fetch(url, {
    headers: { accept: "application/vnd.github+json", "user-agent": "merfisheyes-sync" },
  });

  if (!res.ok) throw new Error(`GitHub ${res.status} for ${url}`);
  const body = await res.json();
  const releases = (Array.isArray(body) ? body : [body]) as {
    tag_name: string;
    assets: { name: string; browser_download_url: string; size: number }[];
  }[];

  for (const r of releases) {
    if (r.tag_name === "catalog-latest") continue; // rolling; not reproducible
    const asset = r.assets.find((a) => a.name === CATALOG_ASSET);

    if (asset) {
      return { tag: r.tag_name, url: asset.browser_download_url, size: asset.size };
    }
  }
  throw new Error(`no release with ${CATALOG_ASSET}${tag ? ` at ${tag}` : ""}`);
}

async function download(asset: ReleaseAsset, cacheDir: string): Promise<string> {
  const dir = join(cacheDir, asset.tag);
  const path = join(dir, CATALOG_ASSET);

  if (existsSync(path)) return path;
  mkdirSync(dir, { recursive: true });
  const res = await fetch(asset.url, { headers: { "user-agent": "merfisheyes-sync" } });

  if (!res.ok || !res.body) throw new Error(`download ${res.status} for ${asset.url}`);
  await pipeline(Readable.fromWeb(res.body as never), createWriteStream(path));
  writeFileSync(join(dir, "source.json"), JSON.stringify(asset, null, 2));

  return path;
}

/** Run `fn` over `items` with at most `n` in flight. */
async function pool<T, R>(items: T[], n: number, fn: (item: T) => Promise<R>): Promise<R[]> {
  const out: R[] = new Array(items.length);
  let next = 0;
  const worker = async () => {
    while (next < items.length) {
      const i = next++;

      out[i] = await fn(items[i]);
    }
  };

  await Promise.all(Array.from({ length: Math.min(n, items.length) }, worker));

  return out;
}

/**
 * Register every store as an admin-owned Dataset and name it from the
 * catalog. Returns uid → dataset id; a store whose root cannot be read is
 * logged and its card still synced, without a dataset id.
 */
async function registerStores(
  mapped: MappedCatalogRow[],
  prisma: typeof import("../lib/prisma").prisma,
): Promise<Map<string, string>> {
  const { ensureDatasetForS3Url } = await import("../lib/datasets/register-or-claim");
  const owner = await prisma.user.findUnique({
    where: { email: ownerEmail },
    select: { id: true, role: true },
  });

  if (!owner) throw new Error(`No user with email ${ownerEmail}`);
  const isAdmin = owner.role === "ADMIN" || owner.role === "SUPER_ADMIN";

  if (!isAdmin) throw new Error(`${ownerEmail} is not an admin; stores would be registered to that account personally`);

  const ids = new Map<string, string>();
  let created = 0;
  let failed = 0;

  await pool(mapped, concurrency, async (m) => {
    try {
      const { dataset, claimed, ownedByOther } = await ensureDatasetForS3Url(m.entry.s3BaseUrl, {
        userId: owner.id,
        isAdmin,
        asAdmin: true,
        kind: "spatialdata",
      });

      ids.set(m.uid, dataset.id);
      if (claimed) created++;
      // Name and describe the dataset from the catalog, unless someone owns
      // it personally or has already given it a title of their own.
      const bareTitle = dataset.title === m.entry.s3BaseUrl.split("/").pop()!.replace(/\.zarr$/, "");

      if (!ownedByOther && (claimed || bareTitle || dataset.title === m.data.title)) {
        await prisma.dataset.update({
          where: { id: dataset.id },
          data: {
            title: m.data.title,
            description: m.data.description,
            species: m.data.species,
            tissue: m.data.tissue,
            platform: m.data.platform,
            externalLink: m.data.externalLink,
            publicationLink: m.data.publicationLink,
            tags: m.data.tags,
            metadata: m.data.metadata,
            // The catalog counts the default table; the registration read
            // whichever table it found first (a Visium HD store's 2 µm bins).
            ...(m.data.numCells != null ? { numCells: m.data.numCells } : {}),
            ...(m.data.numGenes != null ? { numGenes: m.data.numGenes } : {}),
          },
        });
      }
    } catch (e) {
      failed++;
      console.warn(`  !! ${m.uid}: ${(e as Error).message}`);
    }
  });
  console.log(`  registered: ${ids.size} datasets (${created} new, ${failed} failed)`);

  return ids;
}

async function main() {
  const envFile = arg("env");

  if (envFile) loadEnv(envFile, process.argv.includes("--direct"));
  if (!process.env.DATABASE_URL) throw new Error("DATABASE_URL is not set");

  const file = arg("file");
  let tag = arg("tag");
  let path: string;

  if (file) {
    path = file;
    tag ??= "local";
  } else {
    const asset = await findRelease(tag);

    tag = asset.tag;
    path = await download(asset, arg("cache") ?? ".cache/sddb");
  }

  const catalog = await readCatalogFile(path);

  console.log(
    `${catalog.rows.length} rows from ${path} (release ${tag}, generated ${catalog.generatedAt ?? "?"})` +
      (tag !== PINNED_RELEASE ? `; schema pinned at ${PINNED_RELEASE}` : "") +
      (apply ? "" : "  (dry run — pass --apply to write)"),
  );
  if (catalog.missingColumns.length) {
    console.warn(`  known columns missing from this release: ${catalog.missingColumns.join(", ")}`);
  }

  const rows = limit ? sampleRows(catalog.rows, limit) : catalog.rows;
  const mapped = rows.map((row) =>
    mapCatalogRow(row, { release: tag!, generatedAt: catalog.generatedAt }),
  );
  const uids = new Set(mapped.map((m) => m.uid));

  if (uids.size !== mapped.length) throw new Error("duplicate uid in the catalog file");
  if (limit) {
    console.log(`  --limit=${limit}: a sample spread over technologies; nothing is unpublished and the uid-set gate is skipped`);
  }

  // The database only now: the env is loaded, so the client sees DATABASE_URL.
  const { prisma } = await import("../lib/prisma");
  const existing = await prisma.catalogDataset.findMany({
    where: { source: "sddb" },
    select: { id: true, sddbUid: true, isPublished: true, sourceRelease: true },
  });
  const byUid = new Map(existing.map((e) => [e.sddbUid!, e]));
  const stale = limit ? [] : staleUids(byUid.keys(), uids);
  const toCreate = mapped.filter((m) => !byUid.has(m.uid));
  const toUpdate = mapped.filter((m) => byUid.has(m.uid));
  const registered = await prisma.dataset.count({
    where: { s3BaseUrl: { in: mapped.map((m) => m.entry.s3BaseUrl) } },
  });

  console.log(
    `  database: ${existing.length} sddb rows (${existing.filter((e) => e.isPublished).length} published)` +
      `\n  create ${toCreate.length}, update ${toUpdate.length}, unpublish ${stale.length}` +
      (register ? `; register ${mapped.length - registered} stores (${registered} already)` : "; no registration"),
  );
  for (const m of mapped.slice(0, 2)) {
    console.log(`  e.g. ${m.uid}: "${m.data.title}" · ${m.data.species} · ${m.data.tissue} · ${m.data.platform} · ${m.data.license} → ${m.entry.s3BaseUrl}`);
  }
  if (!apply) {
    await prisma.$disconnect();

    return;
  }

  const datasetIds = register ? await registerStores(mapped, prisma) : new Map<string, string>();

  const BATCH = 50;
  let done = 0;

  for (let i = 0; i < mapped.length; i += BATCH) {
    const slice = mapped.slice(i, i + BATCH);

    await prisma.$transaction(
      async (tx) => {
        for (const m of slice) {
          const row = await tx.catalogDataset.upsert({
            where: { sddbUid: m.uid },
            create: { ...m.data },
            update: { ...m.data },
            select: { id: true, entries: { select: { id: true, datasetType: true } } },
          });
          const entry = row.entries.find((e) => e.datasetType === "spatialdata");
          const data = { ...m.entry, datasetId: datasetIds.get(m.uid) ?? null };

          if (entry) {
            await tx.catalogDatasetEntry.update({ where: { id: entry.id }, data });
          } else {
            await tx.catalogDatasetEntry.create({ data: { ...data, catalogId: row.id, sortOrder: 0 } });
          }
        }
      },
      { timeout: 120_000, maxWait: 20_000 },
    );
    done += slice.length;
    if (done % 500 === 0 || done === mapped.length) console.log(`  ${done}/${mapped.length}`);
  }

  if (stale.length) {
    await prisma.catalogDataset.updateMany({
      where: { source: "sddb", sddbUid: { in: stale } },
      data: { isPublished: false },
    });
  }

  // The gate: what we publish is exactly the release's uid set.
  const published = await prisma.catalogDataset.findMany({
    where: { source: "sddb", isPublished: true },
    select: { sddbUid: true },
  });
  const ours = new Set(published.map((p) => p.sddbUid));
  const equal = ours.size === uids.size && [...uids].every((u) => ours.has(u));

  console.log(
    `  published sddb rows: ${ours.size}` +
      (limit ? " (partial sync)" : `; equals release uid set: ${equal}`),
  );
  await prisma.$disconnect();
  if (!limit && !equal) process.exit(1);
}

main().catch((e) => {
  console.error(e);
  process.exit(1);
});
