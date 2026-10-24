"use client";

/**
 * The element registry: every kind of element a scene can draw, as one
 * entry that names how it is found in a store, opened, drawn, kept in step
 * with the panel's stores, and controlled. Viewers iterate `ELEMENTS`
 * instead of naming kinds; adding a kind is adding an entry (and its
 * conformance case): see docs/specs/LAYER-CONTRACT.md.
 */
import type { ComponentType } from "react";
import type { ElementKind, Layer, SceneHandle } from "./contract";
import type { SpatialDataStore } from "@/lib/spatialdata/store";
import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { SpatialDataImage } from "@/lib/spatialdata/image";
import type { CellPolygons } from "@/lib/spatialdata/shapes";
import type { MoleculeLayerOptions } from "@/lib/webgl/MoleculeLayer";
import type { DatasetLinkConfig } from "@/lib/config/dataset-links";

import { PointsStatusPanel } from "@/components/points-status-panel";
import { SpatialDataImageControls } from "@/components/spatialdata-image-controls";
import { SpatialDataShapesControls } from "@/components/spatialdata-shapes-controls";
import { useMoleculeLayerSync } from "@/lib/hooks/useMoleculeLayerSync";
import { usePointsLayerSync } from "@/lib/hooks/usePointsLayerSync";
import { ImageLayer } from "@/lib/spatialdata/ImageLayer";
import { openSpatialDataImage } from "@/lib/spatialdata/image";
import { loadSpatialDataMolecules } from "@/lib/spatialdata/loadSpatialDataMolecules";
import {
  findPointsElements,
  openPointsElement,
} from "@/lib/spatialdata/points";
import { PointsLayer } from "@/lib/spatialdata/PointsLayer";
import {
  pointsElementOf,
  pointsMoleculeDataset,
} from "@/lib/spatialdata/pointsMolecules";
import {
  placesAsWritten,
  readCoordinateSystems,
} from "@/lib/spatialdata/coordinates";
import { cellShapesName, loadCellPolygons } from "@/lib/spatialdata/shapes";
import { ShapesLayer } from "@/lib/spatialdata/ShapesLayer";
import { elementNames } from "@/lib/spatialdata/store";
import { MoleculeLayer } from "@/lib/webgl/MoleculeLayer";

/**
 * What was opened of a store or folder, by kind: the metadata a layer needs
 * to start drawing, never the data itself (tiles, molecules and polygons are
 * read a view at a time by the layer).
 */
export interface OpenedElements {
  /** A multiscale image, drawn under everything else. */
  image?: SpatialDataImage;
  /** Cell polygons, read the first time they are drawn. */
  shapes?: () => Promise<CellPolygons>;
  /** Molecules in the Zarr points format, drawn a view at a time. */
  points?: SingleMoleculeDataset;
  /** Molecules as per-gene files (the legacy folder), drawn whole per gene. */
  molecules?: SingleMoleculeDataset;
  /** Elements the store has that are not drawn, and why. */
  problems?: string[];
}

export interface OpenContext {
  /** Dataset id and display name the opened parts are named after. */
  id: string;
  name: string;
  /** The store's cells, when it has a table: shapes belong to them. */
  cells?: StandardizedDataset;
  onProgress?: (progress: number, message: string) => void;
}

/** What a viewer hands the layers it creates, beyond the element. */
export interface LayerContext {
  molecules?: MoleculeLayerOptions & {
    /** After each reconciliation of the selection with the layer. */
    onSynced?: (layer: MoleculeLayer) => void;
  };
}

export interface ElementDefinition<K extends ElementKind = ElementKind> {
  kind: K;
  /** Names of this kind's elements in a store; empty when it has none. */
  detect(store: SpatialDataStore): Promise<string[]>;
  /** Open one element's metadata; null when it cannot be shown here. */
  open(
    store: SpatialDataStore,
    name: string,
    ctx: OpenContext,
  ): Promise<NonNullable<OpenedElements[K]> | null>;
  /** The layer that draws it. */
  createLayer(opened: NonNullable<OpenedElements[K]>, ctx: LayerContext): Layer;
  /** What the scene must offer; the layer stays out of one that does not. */
  canAttach?(scene: SceneHandle): boolean;
  /** Keeps the layer in step with the panel's stores (renders nothing). */
  Bind?: ComponentType<{ layer: Layer; ctx: LayerContext }>;
  /** Its controls, mounted by `LayerControls` while the layer exists. */
  Controls?: ComponentType<{ layer: Layer }>;
  /** Where the controls go: the viewer's button rail, or a floating panel. */
  slot?: "rail" | "panel";
}

/**
 * Points and polygons are drawn as written, in the cells' coordinates; one
 * whose transforms put it elsewhere is left out rather than drawn wrong
 * (every SpatialData-DB store's sit in the cells' frame: sddb-profile.md).
 */
async function assertInCellsFrame(store: SpatialDataStore, path: string) {
  const placed = (await readCoordinateSystems(store)).place(path);

  if (!placesAsWritten(placed)) {
    throw new Error(
      placed.ok
        ? `${path} declares a different coordinate frame from the cells (scale ${placed.map.scale.join(" × ")}); the store's transforms need fixing`
        : `${path}: ${placed.reason}`,
    );
  }
}

const image: ElementDefinition<"image"> = {
  kind: "image",
  detect: async (store) => elementNames(store.keys, "images"),
  open: (store) => openSpatialDataImage(store),
  createLayer: (opened) => new ImageLayer(opened),
  Controls: SpatialDataImageControls as ComponentType<{ layer: Layer }>,
  slot: "panel",
};

const shapes: ElementDefinition<"shapes"> = {
  kind: "shapes",
  detect: async (store) => {
    const name = cellShapesName(store);

    return name ? [name] : [];
  },
  // Read the first time they are drawn, and kept; a failure can be retried
  open: async (store, name, { cells }) => {
    if (!cells) return null;
    await assertInCellsFrame(store, `shapes/${name}`);
    let polygons: Promise<CellPolygons> | null = null;

    return () => {
      polygons ??= (async () => {
        const ids = await cells.adapter?.readInstanceIds?.();

        if (!ids) throw new Error("The table names no cell ids");

        return loadCellPolygons(store, name, ids);
      })();
      polygons.catch(() => {
        polygons = null;
      });

      return polygons;
    };
  },
  createLayer: (load) => new ShapesLayer(load),
  // Polygons take the colour of their cell
  canAttach: (scene) => !!scene.cellColors,
  Controls: SpatialDataShapesControls as ComponentType<{ layer: Layer }>,
  slot: "rail",
};

function PointsBinding({ layer }: { layer: Layer }) {
  usePointsLayerSync(layer as PointsLayer);

  return null;
}

const points: ElementDefinition<"points"> = {
  kind: "points",
  detect: (store) => findPointsElements(store),
  open: async (store, name, { id, name: datasetName }) => {
    await assertInCellsFrame(store, `points/${name}`);

    return pointsMoleculeDataset(
      await openPointsElement(store, name),
      `${id}_${name}`,
      `${datasetName} — ${name}`,
    );
  },
  createLayer: (opened) => new PointsLayer(pointsElementOf(opened)!),
  Bind: PointsBinding,
  Controls: PointsStatusPanel as ComponentType<{ layer: Layer }>,
  slot: "panel",
};

function MoleculesBinding({ layer, ctx }: { layer: Layer; ctx: LayerContext }) {
  useMoleculeLayerSync(layer as MoleculeLayer, ctx.molecules?.onSynced);

  return null;
}

const molecules: ElementDefinition<"molecules"> = {
  kind: "molecules",
  // Linked by the mapping.json beside the store
  detect: async (store) => {
    const bytes = await store.sidecar.get("/mapping.json");

    if (!bytes) return [];
    const mapping: Partial<DatasetLinkConfig> = JSON.parse(
      new TextDecoder().decode(bytes),
    );

    return mapping.links?.__all__ ? ["__all__"] : [];
  },
  open: (store) => loadSpatialDataMolecules(store),
  createLayer: (opened, ctx) => new MoleculeLayer(opened, ctx.molecules),
  Bind: MoleculesBinding,
};

/** Every kind, in draw order (what is attached first is drawn under). */
export const ELEMENTS: ElementDefinition[] = [
  image as ElementDefinition,
  shapes as ElementDefinition,
  points as ElementDefinition,
  molecules as ElementDefinition,
];

export function elementDefinition<K extends ElementKind>(
  kind: K,
): ElementDefinition<K> {
  return ELEMENTS.find((d) => d.kind === kind) as ElementDefinition<K>;
}

/**
 * Open every element of a store that can be shown with the cells given
 * (or without any). One element that fails to open is logged and left out,
 * so a broken image never hides the cells. A store's own points stand in
 * for a linked molecule folder, so at most one molecule kind comes back.
 */
export async function openElements(
  store: SpatialDataStore,
  ctx: OpenContext,
): Promise<OpenedElements> {
  const opened: OpenedElements = {};

  for (const def of ELEMENTS) {
    if (def.kind === "molecules" && opened.points) continue;
    try {
      const [name] = await def.detect(store);

      if (!name) continue;
      ctx.onProgress?.(90, `Opening ${def.kind}...`);
      const element = await def.open(store, name, ctx);

      if (!element) continue;
      // A linked store's molecules are its points element: drawn as points
      if (
        def.kind === "molecules" &&
        pointsElementOf(element as SingleMoleculeDataset)
      ) {
        opened.points = element as SingleMoleculeDataset;
        continue;
      }
      (opened as Record<string, unknown>)[def.kind] = element;
    } catch (error) {
      console.error(`Error opening ${def.kind}:`, error);
      (opened.problems ??= []).push(
        `${def.kind}: ${error instanceof Error ? error.message : String(error)}`,
      );
    }
  }

  return opened;
}

/**
 * A molecule dataset as an element: a points element when it stands for
 * one, else per-gene files. Nothing for null.
 */
export function moleculeElements(
  ds: SingleMoleculeDataset | null | undefined,
): Pick<OpenedElements, "points" | "molecules"> {
  if (!ds) return {};

  return pointsElementOf(ds) ? { points: ds } : { molecules: ds };
}
