import { prisma } from "@/lib/prisma";
import {
  uploadBufferToS3,
  getObjectJson,
  deleteObject,
  listObjectKeys,
} from "@/lib/s3";
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";
import {
  rootHasPoints,
  type MoleculeLinkFormat,
} from "@/lib/overlays/moleculeSource";

export { rootHasPoints };

/**
 * A single-cell dataset can carry a single-molecule "overlay": a mapping.json
 * in its S3 prefix that tells the viewer to render the molecules of another
 * (app) dataset on top of the cells. This module is the one place that writes,
 * reads, and clears that link, plus the Project that groups the two datasets —
 * shared by the retroactive-overlay API and the duplicate-upload finalize path.
 */

const mappingKey = (scId: string) => `datasets/${scId}/mapping.json`;

/**
 * Build the overlay mapping.json for a molecule dataset. How the viewer
 * loads it depends on the dataset's storage:
 *   - app-uploaded (chunked, sm_ id)        → source "app", ref = dataset id
 *     (loaded via SingleMoleculeDataset.fromS3 by id / presigned API)
 *   - S3-registered (has s3BaseUrl, s3_ id) → source "s3", ref = the public URL
 *     (loaded via fromCustomS3) — fromS3(id) would 404 for these.
 *   - a SpatialData store (what scripts/prepare.py makes of a molecule
 *     table) → `linkFormat: "store"`: its tiled points element is drawn
 * `datasetId` is always stored so the UI can identify the overlaid dataset
 * regardless of source.
 */
/** A dataset registered by its public URL (s3_), not uploaded to the app. */
function isS3Registered(
  id: string,
  d: { ingestSource: string | null; s3BaseUrl: string | null } | null,
): boolean {
  return (
    d?.ingestSource === "s3_registered" ||
    (!!d?.s3BaseUrl && id.startsWith("s3_"))
  );
}

export async function moleculeLinkFor(smId: string) {
  const sm = await prisma.dataset.findUnique({
    where: { id: smId },
    select: {
      id: true,
      s3BaseUrl: true,
      ingestSource: true,
      datasetType: true,
    },
  });
  const store =
    sm?.datasetType === "spatialdata" ? { linkFormat: "store" as const } : {};

  if (isS3Registered(smId, sm) && sm?.s3BaseUrl) {
    return {
      linkColumn: "__all__",
      links: { __all__: sm.s3BaseUrl },
      source: "s3" as const,
      datasetId: smId,
      ...store,
    };
  }

  return {
    linkColumn: "__all__",
    links: { __all__: smId },
    source: "app" as const,
    datasetId: smId,
    ...store,
  };
}

/** The molecule link's fields of a mapping.json; the rest is someone else's. */
const MOLECULE_FIELDS = [
  "linkColumn",
  "links",
  "source",
  "datasetId",
  "linkFormat",
];

/**
 * Whether a dataset is a SpatialData store with molecules, whatever else it
 * holds: a store with cells too lends just its molecules. An uploaded store
 * is looked for in the app's bucket; one registered by its public URL is read
 * there, a folder resolved to the store it names as the viewer resolves it.
 */
export async function storeHasPoints(datasetId: string): Promise<boolean> {
  const ds = await prisma.dataset.findUnique({
    where: { id: datasetId },
    select: { s3BaseUrl: true, ingestSource: true },
  });

  if (isS3Registered(datasetId, ds) && ds?.s3BaseUrl) {
    const target = await resolveZarrUrl(ds.s3BaseUrl);

    return target?.kind === "spatialdata" && rootHasPoints(target.root);
  }
  const keys = await listObjectKeys(`datasets/${datasetId}/`);

  return keys.some((k) => /\.zarr\/points\/[^/]+\/zarr\.json$/.test(k.key));
}

/** A dataset the user owns, or, for an admin, one the admins own; else null. */
export async function ownedDataset(
  id: string,
  userId: string,
  isAdmin: boolean,
) {
  const d = await prisma.dataset.findUnique({
    where: { id },
    select: {
      id: true,
      ownerId: true,
      adminOwned: true,
      datasetType: true,
      title: true,
    },
  });

  if (!d) return null;
  const owned = d.ownerId === userId || (isAdmin && d.adminOwned);

  return owned ? d : null;
}

/**
 * Whether a dataset has molecules an overlay can draw: a single molecule
 * dataset, or a SpatialData store with them (what scripts/prepare.py makes of
 * a molecule table, or a store with cells too).
 */
export async function hasMolecules(d: {
  id: string;
  datasetType: string | null;
}): Promise<boolean> {
  return (
    d.datasetType === "single_molecule" ||
    (d.datasetType === "spatialdata" && (await storeHasPoints(d.id)))
  );
}

/** The single-molecule dataset id currently overlaid on `scId`, or null. */
export async function getOverlaySmId(scId: string): Promise<string | null> {
  const m = await getObjectJson<{
    datasetId?: string;
    source?: string;
    links?: Record<string, string>;
  }>(mappingKey(scId)).catch(() => null);

  if (!m) return null;

  // Prefer the explicit datasetId; fall back to links.__all__ for app-source
  // mappings written before datasetId was stored.
  return m.datasetId ?? (m.source !== "s3" ? (m.links?.__all__ ?? null) : null);
}

/**
 * The public URL of the molecules overlaid on `scId` when they were added by
 * their link rather than as a dataset (writeOverlayUrl), or null.
 */
export async function getOverlayUrl(scId: string): Promise<string | null> {
  const m = await getObjectJson<{
    datasetId?: string;
    source?: string;
    links?: Record<string, string>;
  }>(mappingKey(scId)).catch(() => null);

  return m && !m.datasetId && m.source === "s3"
    ? (m.links?.__all__ ?? null)
    : null;
}

/**
 * The image drawn under a dataset's cells, as its mapping.json names it:
 *   - an image dataset in the app's bucket → `image` is its id, `imageSource`
 *     "app" (opened by id through the presigned API);
 *   - an image at a public URL → `image` is the URL (`imageSource` "s3" or
 *     absent, which is how mappings written before image datasets look).
 * `imageDatasetId` is set whenever the image is a registered dataset, so the
 * UI can say which one it is.
 */
export interface ImageOverlay {
  image: string;
  imageSource?: "s3" | "app";
  imageDatasetId?: string;
}

const imageFields = (m: Partial<ImageOverlay> | null): ImageOverlay | null =>
  m?.image
    ? {
        image: m.image,
        ...(m.imageSource ? { imageSource: m.imageSource } : {}),
        ...(m.imageDatasetId ? { imageDatasetId: m.imageDatasetId } : {}),
      }
    : null;

/** The image overlay in the dataset's mapping.json, if any. */
export async function getOverlayImage(
  scId: string,
): Promise<ImageOverlay | null> {
  return imageFields(
    await getObjectJson<Partial<ImageOverlay>>(mappingKey(scId)).catch(
      () => null,
    ),
  );
}

/** How a mapping.json names an image dataset; null when there is no such image. */
export async function imageOverlayFor(
  imageDatasetId: string,
): Promise<ImageOverlay | null> {
  const img = await prisma.dataset.findUnique({
    where: { id: imageDatasetId },
    select: { datasetType: true, s3BaseUrl: true, ingestSource: true },
  });

  if (img?.datasetType !== "image") return null;

  return img.ingestSource === "s3_registered" && img.s3BaseUrl
    ? { image: img.s3BaseUrl, imageSource: "s3", imageDatasetId }
    : { image: imageDatasetId, imageSource: "app", imageDatasetId };
}

/** mapping.json without its molecule link: the image, a store's `cells`… */
async function mappingWithoutMolecules(
  scId: string,
): Promise<Record<string, unknown>> {
  const mapping: Record<string, unknown> =
    (await getObjectJson<Record<string, unknown>>(mappingKey(scId)).catch(
      () => null,
    )) ?? {};

  for (const field of MOLECULE_FIELDS) delete mapping[field];

  return mapping;
}

// The image overlay shares mapping.json with the molecule overlay but is
// independent of it, so writing or removing the molecules keeps the image —
// and a store's own `cells` entry
/**
 * Overlay the molecules at a public URL, by their link alone: no dataset is
 * made of them. `linkFormat` is what moleculeUrlFormat found there.
 */
export async function writeOverlayUrl(
  scId: string,
  url: string,
  linkFormat: MoleculeLinkFormat,
): Promise<void> {
  const rest = await mappingWithoutMolecules(scId);

  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(
      JSON.stringify({
        ...rest,
        linkColumn: "__all__",
        links: { __all__: url },
        source: "s3",
        linkFormat,
      }),
    ),
    "application/json",
  );
}

export async function writeOverlay(scId: string, smId: string): Promise<void> {
  const mapping = await moleculeLinkFor(smId);
  const rest = await mappingWithoutMolecules(scId);

  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify({ ...rest, ...mapping })),
    "application/json",
  );
}

export async function removeOverlay(scId: string): Promise<void> {
  const rest = await mappingWithoutMolecules(scId);

  if (Object.keys(rest).length > 0) {
    await uploadBufferToS3(
      mappingKey(scId),
      Buffer.from(JSON.stringify(rest)),
      "application/json",
    );

    return;
  }
  await deleteObject(mappingKey(scId));
}

/**
 * Set the image overlay or, with null, clear it. The molecule overlay in the
 * same mapping.json is left as it is.
 */
export async function writeOverlayImage(
  scId: string,
  image: ImageOverlay | null,
): Promise<void> {
  const mapping: Record<string, unknown> =
    (await getObjectJson<Record<string, unknown>>(mappingKey(scId)).catch(
      () => null,
    )) ?? {};

  delete mapping.image;
  delete mapping.imageSource;
  delete mapping.imageDatasetId;
  Object.assign(mapping, imageFields(image));

  if (Object.keys(mapping).length === 0) {
    await deleteObject(mappingKey(scId));

    return;
  }
  await uploadBufferToS3(
    mappingKey(scId),
    Buffer.from(JSON.stringify(mapping)),
    "application/json",
  );
}

/**
 * Group two datasets in a Project owned by `ownerId`. Reuses an existing
 * project that already contains the single-cell dataset if one exists;
 * otherwise creates one titled `title`. Returns the project id.
 */
export async function groupInProject(
  ownerId: string,
  title: string,
  scId: string,
  smId: string,
): Promise<string> {
  const existing = await prisma.projectDataset.findFirst({
    where: { datasetId: scId, project: { ownerId } },
    select: { projectId: true },
  });

  const projectId =
    existing?.projectId ??
    (await prisma.project.create({ data: { title, ownerId } })).id;

  for (const [i, datasetId] of [scId, smId].entries()) {
    await prisma.projectDataset.upsert({
      where: { projectId_datasetId: { projectId, datasetId } },
      update: {},
      create: { projectId, datasetId, sortOrder: i },
    });
  }

  return projectId;
}
