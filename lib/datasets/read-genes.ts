import { gunzipSync } from "zlib";

import { getObjectBytes } from "@/lib/s3";

/**
 * Extract a dataset's gene symbols server-side, from artifacts that already
 * exist when processing completes — no client involvement and no worker/image
 * change:
 *
 *   - single molecule → the manifest's gene list. Browser-written manifests
 *     carry `uniqueGenes`; Python/worker-written ones carry
 *     `genes.unique_gene_names`. Read from the row's stored manifestJson when
 *     present, else the manifest object itself.
 *   - single cell (chunked) → expr/index.json, whose `genes` entries carry
 *     every gene name.
 *
 * Objects are read from the app bucket at datasets/{id}/… for uploaded
 * datasets, or over plain HTTPS from the dataset's public s3BaseUrl for
 * s3-registered rows (their data lives in an external bucket).
 *
 * Returns null when genes can't be determined (e.g. zarr-format datasets,
 * missing objects) — callers treat that as "leave the row as is". Used by the
 * upload complete routes, the ingest finalize, and the backfill script.
 */
export interface GenesSource {
  id: string;
  datasetType: string | null;
  formatVersion?: string | null;
  manifestJson?: unknown;
  ingestSource?: string | null;
  s3BaseUrl?: string | null;
}

export async function readDatasetGenes(
  dataset: GenesSource,
): Promise<string[] | null> {
  try {
    if (dataset.datasetType === "single_molecule") {
      return await readSingleMoleculeGenes(dataset);
    }

    // Zarr-format single cell stores no expr/index.json; skip for now.
    if (dataset.formatVersion === "zarr") return null;

    return await readChunkedGenes(dataset);
  } catch {
    return null;
  }
}

function sanitizeGenes(raw: unknown): string[] | null {
  if (!Array.isArray(raw)) return null;
  const genes = Array.from(
    new Set(
      raw.filter((g): g is string => typeof g === "string" && g.length > 0),
    ),
  );

  return genes.length > 0 ? genes : null;
}

/** Gene list from a manifest of either lineage (browser or Python). */
function genesFromManifest(manifest: unknown): string[] | null {
  const m = manifest as {
    uniqueGenes?: unknown;
    genes?: { unique_gene_names?: unknown } | unknown;
  } | null;

  return (
    sanitizeGenes(m?.uniqueGenes) ??
    sanitizeGenes(
      (m?.genes as { unique_gene_names?: unknown } | undefined)
        ?.unique_gene_names,
    )
  );
}

/**
 * Fetch a dataset-relative object: from the app bucket for uploads, or the
 * public s3BaseUrl for s3-registered rows. Returns raw bytes or null.
 */
async function fetchDatasetObject(
  dataset: GenesSource,
  relKey: string,
): Promise<Buffer | null> {
  if (dataset.ingestSource === "s3_registered" && dataset.s3BaseUrl) {
    try {
      const res = await fetch(
        `${dataset.s3BaseUrl.replace(/\/+$/, "")}/${relKey}`,
      );

      if (!res.ok) return null;

      return Buffer.from(await res.arrayBuffer());
    } catch {
      return null;
    }
  }

  const obj = await getObjectBytes(`datasets/${dataset.id}/${relKey}`);

  return obj?.body ?? null;
}

/** Transparently gunzip when the payload is gzipped. */
function toText(buf: Buffer): string {
  return buf[0] === 0x1f && buf[1] === 0x8b
    ? gunzipSync(buf).toString("utf-8")
    : buf.toString("utf-8");
}

async function readSingleMoleculeGenes(
  dataset: GenesSource,
): Promise<string[] | null> {
  // Preferred: the manifest stored on the row at initiate time.
  const fromRow = genesFromManifest(dataset.manifestJson);

  if (fromRow) return fromRow;

  // Fallback: the manifest object itself.
  const buf = await fetchDatasetObject(dataset, "manifest.json.gz");

  if (!buf) return null;

  return genesFromManifest(JSON.parse(toText(buf)));
}

async function readChunkedGenes(
  dataset: GenesSource,
): Promise<string[] | null> {
  const buf = await fetchDatasetObject(dataset, "expr/index.json");

  if (!buf) return null;
  const index = JSON.parse(toText(buf)) as { genes?: { name?: unknown }[] };

  return sanitizeGenes(index.genes?.map((g) => g?.name));
}

/**
 * Best-effort: read the genes and store them on the dataset row. Never throws
 * — completion/callback flows must not fail because gene indexing did.
 */
export async function updateDatasetGenes(
  prisma: {
    dataset: {
      update: (args: {
        where: { id: string };
        data: { genes: string[] };
      }) => Promise<unknown>;
    };
  },
  dataset: GenesSource,
): Promise<boolean> {
  try {
    const genes = await readDatasetGenes(dataset);

    if (!genes) return false;
    await prisma.dataset.update({
      where: { id: dataset.id },
      data: { genes },
    });

    return true;
  } catch {
    return false;
  }
}
