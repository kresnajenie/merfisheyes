/**
 * What is drawn with a cell dataset besides its cells, found the way the
 * dataset was found: the parts dropped beside it for a preview, the
 * hardcoded link registry, the `mapping.json` beside a public folder, or
 * the app's own mapping route for an uploaded dataset. Each part is handed
 * over as soon as it is open (an image need not wait for the molecules).
 *
 * Pure apart from what `io` does, so it can be tested with fakes; the
 * hook `useOverlays` wires it to the stores.
 */
import type { OpenedElements } from "@/lib/scene/registry";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type {
  DatasetLinkConfig,
  ImageOnlyLinkConfig,
} from "@/lib/config/dataset-links";
import type { SpatialDataImage } from "@/lib/spatialdata/image";

import { moleculeElements } from "@/lib/scene/registry";

export interface OverlayIO {
  /** The mapping.json beside a public folder, or null. */
  fetchMapping(
    dataset: StandardizedDataset,
  ): Promise<DatasetLinkConfig | ImageOnlyLinkConfig | null>;
  /** The app's mapping route for an uploaded dataset, or null. */
  fetchAppMapping(datasetId: string): Promise<DatasetLinkConfig | null>;
  openImage(ref: {
    image: string;
    source?: "s3" | "app";
  }): Promise<SpatialDataImage>;
  loadMolecules(
    ref: string,
    source: "s3" | "app",
    linkFormat?: DatasetLinkConfig["linkFormat"],
  ): Promise<SingleMoleculeDataset>;
  /** Make a cluster column available on the dataset; true when fetched. */
  loadColumn(dataset: StandardizedDataset, column: string): Promise<boolean>;
}

export interface OverlayEvents {
  /** A part is open; drawn from now on. */
  element(patch: Partial<OpenedElements>): void;
  /** Which molecules the cells link to (the right-click split). */
  link(config: DatasetLinkConfig | null): void;
  /** The molecules of an `__all__` link are being read. */
  loading(on: boolean): void;
  /** A public URL the molecules came from (their saved defaults). */
  source(url: string): void;
  /** The link column was fetched onto the dataset. */
  column(): void;
  /** Something to tell the user. */
  warn(message: string): void;
  /** True once the caller has moved on: nothing more is handed over. */
  stale(): boolean;
}

/** The parts dropped beside the cells for a preview, already open. */
export interface LocalOverlays {
  image: SpatialDataImage | null;
  shapes: (() => Promise<unknown>) | null;
  molecules: SingleMoleculeDataset | null;
}

const imageOnly = (
  mapping: DatasetLinkConfig | ImageOnlyLinkConfig,
): mapping is ImageOnlyLinkConfig => !("links" in mapping) || !mapping.links;

export async function resolveOverlays(
  dataset: StandardizedDataset,
  local: LocalOverlays | null,
  registry: DatasetLinkConfig | null,
  io: OverlayIO,
  on: OverlayEvents,
): Promise<void> {
  // An image is placed in raw coordinates; a normalized scene has none
  const raw = dataset.normalized === false;

  if (local) {
    if (local.image && !raw) {
      on.warn(
        "The image is not drawn: these cells were saved with normalized coordinates, and an image is placed in the original ones.",
      );
    }
    on.element({
      image: raw ? (local.image ?? undefined) : undefined,
      shapes: (local.shapes as OpenedElements["shapes"]) ?? undefined,
      ...moleculeElements(local.molecules),
    });
    on.link(null);

    return;
  }

  on.link(registry);
  if (registry) return;

  const image = async (
    mapping: { image?: string; imageSource?: "s3" | "app" } | null,
  ) => {
    if (!mapping?.image) return;
    if (!raw) {
      console.warn("mapping.json image skipped: dataset is normalized");

      return;
    }
    try {
      const opened = await io.openImage({
        image: mapping.image,
        source: mapping.imageSource,
      });

      if (!on.stale()) on.element({ image: opened });
    } catch (error) {
      console.warn("Failed to load overlay image:", error);
    }
  };

  const everyCell = async (
    ref: string,
    source: "s3" | "app",
    linkFormat?: DatasetLinkConfig["linkFormat"],
  ) => {
    if (!ref) return;
    on.loading(true);
    try {
      const molecules = await io.loadMolecules(ref, source, linkFormat);

      if (on.stale()) return;
      on.element(moleculeElements(molecules));
      if (source === "s3") on.source(ref);
    } catch (error) {
      console.warn("Failed to load SM overlay dataset:", error);
    } finally {
      on.loading(false);
    }
  };

  // Only a public folder has a mapping.json beside it
  if (dataset.metadata?.customS3BaseUrl) {
    const mapping = await io.fetchMapping(dataset);

    if (!mapping || on.stale()) return;
    void image(mapping);
    if (imageOnly(mapping)) return;
    on.link(mapping);
    // "__all__" means every cell links to the same molecule set: an overlay
    if (mapping.linkColumn === "__all__") {
      await everyCell(
        mapping.links.__all__,
        mapping.source ?? "s3",
        mapping.linkFormat,
      );

      return;
    }
    // Otherwise the link column decides per cell (the right-click split);
    // fetched now so the split opens at once
    if (dataset.clusters?.some((c) => c.column === mapping.linkColumn)) return;
    if (
      dataset.allClusterColumnNames &&
      !dataset.allClusterColumnNames.includes(mapping.linkColumn)
    ) {
      console.warn(
        `mapping.json linkColumn "${mapping.linkColumn}" not found in dataset columns`,
      );

      return;
    }
    try {
      if ((await io.loadColumn(dataset, mapping.linkColumn)) && !on.stale()) {
        on.column();
      }
    } catch (error) {
      console.warn(
        `Failed to preload link column "${mapping.linkColumn}":`,
        error,
      );
    }

    return;
  }

  // An uploaded dataset opens by id: the app's own mapping route, where
  // only the "__all__" overlay is written (the combined SC+SM upload)
  if (dataset.id) {
    const mapping = await io.fetchAppMapping(dataset.id);

    if (!mapping || on.stale()) return;
    void image(mapping);
    if (mapping.linkColumn === "__all__" && mapping.links) {
      on.link(mapping);
      await everyCell(
        mapping.links.__all__,
        mapping.source ?? "s3",
        mapping.linkFormat,
      );
    }
  }
}
