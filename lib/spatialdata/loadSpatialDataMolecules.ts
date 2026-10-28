import type { SpatialDataStore } from "./store";
import type { DatasetLinkConfig } from "@/lib/config/dataset-links";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

/**
 * The molecules a link points at (`mapping.json`'s `links`), however they are
 * kept: a molecule folder of one file per gene, or — `linkFormat: "store"` —
 * a SpatialData store whose tiled points element stands for them (drawn by
 * the scene's PointsLayer). `source` says whether `ref` is a public URL
 * ("s3") or an app dataset id ("app").
 */
export async function openLinkedMolecules(
  ref: string,
  source: "s3" | "app",
  linkFormat: DatasetLinkConfig["linkFormat"] = "molecules",
): Promise<SingleMoleculeDataset> {
  if (linkFormat === "store") {
    const { openSpatialDataStore, openSpatialDataStoreFromUrl } = await import(
      "./store"
    );
    const { findPointsElements, openPointsElement } = await import("./points");
    const { pointsMoleculeDataset } = await import("./pointsMolecules");
    const { resolveZarrUrl } = await import("./resolveUrl");
    // A store registered by its public URL may be a folder whose mapping.json
    // names the store, as /spatialdata-viewer/from-s3 accepts; only its points
    // are drawn, whatever cells it has too
    const store =
      source === "app"
        ? await openSpatialDataStore(ref)
        : await openSpatialDataStoreFromUrl(
            (await resolveZarrUrl(ref))?.url ?? ref,
          );
    const [name] = await findPointsElements(store);

    if (!name) throw new Error(`${ref} holds no molecules`);

    return pointsMoleculeDataset(
      await openPointsElement(store, name),
      `${ref}_${name}`,
      ref,
    );
  }
  const { SingleMoleculeDataset } = await import("@/lib/SingleMoleculeDataset");

  return source === "app"
    ? SingleMoleculeDataset.fromS3(ref)
    : SingleMoleculeDataset.fromCustomS3(ref);
}

/**
 * The molecules linked by the mapping.json beside a store: drawn over its
 * cells, or on their own when it has none. Null when it links none.
 */
export async function loadSpatialDataMolecules({
  sidecar,
}: SpatialDataStore): Promise<SingleMoleculeDataset | null> {
  const bytes = await sidecar.get("/mapping.json");

  if (!bytes) return null;
  const mapping: Partial<DatasetLinkConfig> = JSON.parse(
    new TextDecoder().decode(bytes),
  );
  const molecules = mapping.links?.__all__;

  if (!molecules) return null;

  return openLinkedMolecules(
    molecules,
    mapping.source ?? "s3",
    mapping.linkFormat,
  );
}
