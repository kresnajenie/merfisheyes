/**
 * The ONE way a dataset is opened for viewing, whatever surface asks: a full
 * page by id or by public URL, a split panel, a SpatialData store. Every
 * caller gets the same steps in the same order — look the dataset up,
 * reuse it when it is already open, read it, apply its open state (see
 * `open-dataset.ts`) — and the same result shape, so a viewer only decides
 * what to render.
 */
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { VisualizationState } from "@/lib/stores/createVisualizationStore";
import type { SingleMoleculeVisualizationState } from "@/lib/stores/createSingleMoleculeVisualizationStore";
import type { OpenedElements } from "@/lib/scene/registry";
import type { DatasetLinkConfig } from "@/lib/config/dataset-links";

import { applyCellOpenState, applySMOpenState } from "./open-dataset";

import { tryReadCellVizFromUrl } from "@/lib/hooks/useUrlVizSync";
import { openElements } from "@/lib/scene/registry";
import { loadSpatialDataCells } from "@/lib/spatialdata/loadSpatialDataCells";
import { resolveZarrUrl } from "@/lib/spatialdata/resolveUrl";
import {
  elementNames,
  openSpatialDataStore,
  openSpatialDataStoreFromUrl,
} from "@/lib/spatialdata/store";
import { normalizeS3Url, type ViewerConfig } from "@/lib/utils/viewer-config";

export type DatasetKind = "cell" | "sm" | "spatialdata";

/** Where a dataset comes from: the app's bucket by id, or a public URL. */
export interface DatasetSource {
  kind: DatasetKind;
  datasetId?: string | null;
  s3Url?: string | null;
}

/** What the database knows of the dataset (the registration store's shape). */
export interface DatasetRegistration {
  dbId: string | null;
  ownerId: string | null;
  adminOwned: boolean;
  registered: boolean;
  viewerConfig: ViewerConfig | null;
  s3Url: string | null;
  kind: "image" | "spatialdata" | null;
}

export interface LoadedDataset {
  /** The cells, when the dataset has any: they are the scene. */
  cells: StandardizedDataset | null;
  /** Molecules shown on their own, when there are no cells. */
  molecules: SingleMoleculeDataset | null;
  /** What else the scene draws (an image, polygons, molecules over cells). */
  elements: OpenedElements;
  /** Which molecules the cells link to per cell (the right-click split). */
  linkConfig?: DatasetLinkConfig | null;
  registration: DatasetRegistration;
  /** A SpatialData store's own URL, once a folder URL has been resolved. */
  storeUrl: string | null;
}

export interface LoadDatasetOptions {
  panel: "left" | "right";
  viz: VisualizationState;
  smViz: SingleMoleculeVisualizationState;
  /** Datasets already open in this tab, to reuse instead of reading again. */
  findCells: (id: string) => StandardizedDataset | undefined;
  findMolecules: (id: string) => SingleMoleculeDataset | undefined;
  onProgress: (progress: number, message: string) => void;
  /** True once the caller has moved on; nothing is applied after that. */
  isCancelled: () => boolean;
}

const unregistered = (s3Url: string | null): DatasetRegistration => ({
  dbId: null,
  ownerId: null,
  adminOwned: false,
  registered: false,
  viewerConfig: null,
  s3Url,
  kind: null,
});

/** The registration of a public URL; unregistered when nobody has added it. */
async function registrationByUrl(url: string): Promise<DatasetRegistration> {
  const s3Url = normalizeS3Url(url);

  try {
    const r = await fetch(
      `/api/datasets/by-url?url=${encodeURIComponent(s3Url)}`,
    );

    if (r.ok) {
      const j = await r.json();

      if (j.registered) {
        return {
          dbId: j.id,
          ownerId: j.ownerId ?? null,
          adminOwned: !!j.adminOwned,
          registered: true,
          viewerConfig: (j.viewerConfig as ViewerConfig | null) ?? null,
          s3Url,
          kind: null,
        };
      }
    }
  } catch {
    // best-effort; the viewer still loads unregistered
  }

  return unregistered(s3Url);
}

/**
 * The metadata route of an app dataset. 202 means it is still being
 * uploaded or processed; that is said plainly instead of failing later.
 */
async function datasetMeta(
  datasetId: string,
  route: "datasets" | "single-molecule",
): Promise<Record<string, any>> {
  const r = await fetch(`/api/${route}/${datasetId}`);
  const body = await r.json().catch(() => ({}));

  if (r.status === 202) {
    throw new Error(
      body.message || "Dataset is not ready yet. Try again shortly.",
    );
  }
  if (!r.ok) {
    throw new Error(body.message || `Dataset fetch failed: ${r.status}`);
  }

  return body;
}

const registrationOf = (
  datasetId: string,
  meta: Record<string, any>,
): DatasetRegistration => ({
  dbId: datasetId,
  ownerId: meta.ownerId ?? null,
  adminOwned: !!meta.adminOwned,
  registered: true,
  viewerConfig: (meta.viewerConfig as ViewerConfig | null) ?? null,
  s3Url: null,
  kind: null,
});

export async function loadDataset(
  source: DatasetSource,
  options: LoadDatasetOptions,
): Promise<LoadedDataset | null> {
  switch (source.kind) {
    case "cell":
      return loadCells(source, options);
    case "sm":
      return loadMolecules(source, options);
    case "spatialdata":
      return loadStore(source, options);
  }
}

async function loadCells(
  { datasetId, s3Url }: DatasetSource,
  options: LoadDatasetOptions,
): Promise<LoadedDataset | null> {
  const { panel, viz, findCells, onProgress, isCancelled } = options;

  if (!datasetId && !s3Url) throw new Error("No dataset given");
  let registration: DatasetRegistration;
  let meta: Record<string, any> | null = null;
  let dataset = s3Url ? undefined : findCells(datasetId!);

  if (s3Url) {
    registration = await registrationByUrl(s3Url);
  } else {
    try {
      meta = await datasetMeta(datasetId!, "datasets");
      registration = registrationOf(datasetId!, meta);
    } catch (error) {
      // One that is already open (just uploaded) shows without its record
      if (!dataset) throw error;
      registration = unregistered(null);
    }
  }
  // A SpatialData store reached by a cell-viewer link (the upload bar gave
  // those out) opens as the store it is, in the same scene
  if (!dataset && meta?.datasetType === "spatialdata") {
    return loadStore({ kind: "spatialdata", datasetId }, options);
  }
  const config = registration.viewerConfig;

  if (!dataset) {
    const { StandardizedDataset } = await import("@/lib/StandardizedDataset");
    // Priority column precedence: URL state > owner-saved config > heuristic
    const priorityColumn =
      tryReadCellVizFromUrl(panel)?.c || config?.priorityColumn || undefined;

    dataset = s3Url
      ? await StandardizedDataset.fromCustomS3(
          normalizeS3Url(s3Url),
          onProgress,
          priorityColumn,
        )
      : meta?.formatVersion === "zarr"
        ? await StandardizedDataset.fromS3Zarr(
            datasetId!,
            onProgress,
            priorityColumn,
          )
        : await StandardizedDataset.fromS3(
            datasetId!,
            onProgress,
            priorityColumn,
          );
  }
  if (isCancelled()) return null;
  // Awaited so the first frame paints already-aligned coordinates
  onProgress(100, "Applying alignment...");
  await applyCellOpenState({ dataset, config, store: viz });

  return {
    cells: dataset,
    molecules: null,
    elements: {},
    registration,
    storeUrl: null,
  };
}

async function loadMolecules(
  { datasetId, s3Url }: DatasetSource,
  { panel, smViz, findMolecules, onProgress, isCancelled }: LoadDatasetOptions,
): Promise<LoadedDataset | null> {
  if (!datasetId && !s3Url) throw new Error("No dataset given");
  let registration: DatasetRegistration;

  if (s3Url) {
    registration = await registrationByUrl(s3Url);
  } else {
    // Best-effort: one just uploaded is in the store before it has a record
    registration = await datasetMeta(datasetId!, "single-molecule")
      .then((meta) => registrationOf(datasetId!, meta))
      .catch(() => unregistered(null));
  }
  let dataset = s3Url ? undefined : findMolecules(datasetId!);

  if (!dataset) {
    const { SingleMoleculeDataset } = await import(
      "@/lib/SingleMoleculeDataset"
    );

    dataset = s3Url
      ? await SingleMoleculeDataset.fromCustomS3(
          normalizeS3Url(s3Url),
          onProgress,
        )
      : await SingleMoleculeDataset.fromS3(datasetId!, onProgress);
  }
  if (isCancelled()) return null;
  applySMOpenState({
    dataset,
    config: registration.viewerConfig,
    store: smViz,
    panel,
  });

  return {
    cells: null,
    molecules: dataset,
    elements: {},
    registration,
    storeUrl: null,
  };
}

/**
 * A SpatialData store: its table makes the cells' scene, with the store's
 * other elements drawn in it; without a table, its points element or a
 * linked molecule folder makes the molecules' scene, the image under them;
 * an image alone is shown on its own.
 */
async function loadStore(
  { datasetId, s3Url }: DatasetSource,
  { viz, onProgress, isCancelled }: LoadDatasetOptions,
): Promise<LoadedDataset | null> {
  onProgress(0, "Opening store...");
  let store;
  let id: string;
  let name: string;
  let folder: string | undefined;
  let storeUrl: string | null = null;
  let registration: DatasetRegistration;

  if (s3Url) {
    const resolved = await resolveZarrUrl(s3Url);
    const base = resolved?.url ?? s3Url.replace(/\/+$/, "");

    storeUrl = base;
    folder = base.slice(0, base.lastIndexOf("/"));
    name = base.slice(folder.length + 1).replace(/\.zarr$/, "");
    id = `spatialdata_${name}`;
    store = await openSpatialDataStoreFromUrl(base, resolved?.root);
    // One nobody has registered still opens; the banner offers to add it
    registration = { ...(await registrationByUrl(base)), kind: "spatialdata" };
  } else if (datasetId) {
    const meta = await datasetMeta(datasetId, "datasets");

    id = datasetId;
    name = meta.title || datasetId;
    store = await openSpatialDataStore(datasetId);
    registration = { ...registrationOf(datasetId, meta), kind: "spatialdata" };
  } else {
    throw new Error("No dataset given");
  }
  if (isCancelled()) return null;

  let cells: StandardizedDataset | null = null;

  if (elementNames(store.keys, "tables").length > 0) {
    cells = await loadSpatialDataCells(store, id, name, onProgress, folder);
    if (isCancelled()) return null;
    await applyCellOpenState({ dataset: cells, config: null, store: viz });
  }
  const elements = await openElements(store, {
    id,
    name,
    cells: cells ?? undefined,
    onProgress,
  });

  if (isCancelled()) return null;
  let molecules: SingleMoleculeDataset | null = null;

  if (!cells) {
    molecules = elements.points ?? elements.molecules ?? null;
    delete elements.molecules;
    delete elements.points;
  }
  if (!cells && !molecules && !elements.image) {
    throw new Error(
      "This SpatialData store has no cells, molecules or image to show",
    );
  }

  return { cells, molecules, elements, registration, storeUrl };
}
