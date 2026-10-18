"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";
import type { SceneHandle } from "@/lib/scene/contract";
import type { LayerContext, OpenedElements } from "@/lib/scene/registry";
import type { MoleculeLayer } from "@/lib/webgl/MoleculeLayer";

import { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";
import { toast } from "react-toastify";
import gsap from "gsap";

import { getClusterValue } from "@/lib/StandardizedDataset";
import { initializeScene } from "@/lib/webgl/scene-manager";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";
import { useLocalBundleStore } from "@/lib/stores/localBundleStore";
import {
  createPointCloudFromBuffers,
  updatePointCloudAttributes,
  updateDotSize,
} from "@/lib/webgl/point-cloud";
import { LayerHost } from "@/lib/scene/LayerHost";
import { moleculeElements } from "@/lib/scene/registry";
import { useLayers } from "@/lib/scene/useLayers";
import { LayerBindings, LayerControls } from "@/components/layer-controls";
import {
  updateGeneVisualization,
  updateCelltypeVisualization,
  updateNumericalCelltypeVisualization,
  updateCombinedVisualization,
  updateCoexpressionVisualization,
  type AdvancedVizSettings,
} from "@/lib/webgl/visualization-utils";
import {
  usePanelVisualizationStore,
  usePanelId,
  usePanelSingleMoleculeApi,
  usePanelSingleMoleculeVisualizationApi,
} from "@/lib/hooks/usePanelStores";
import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import {
  getDatasetLinkConfig,
  fetchMappingConfig,
  type DatasetLinkConfig,
} from "@/lib/config/dataset-links";
import { VisualizationLegends } from "@/components/visualization-legends";
import { SingleMoleculeLegends } from "@/components/single-molecule-legends";
import { getEffectiveColumnType } from "@/lib/utils/column-type-utils";
import {
  detectExpressionKind,
  nativeModeFor,
  transformExpression,
} from "@/lib/utils/expression-scale";
import { ExportBoxOverlay } from "@/components/export-box-overlay";
import { WebGLErrorOverlay } from "@/components/webgl-error-overlay";
import {
  markSceneReady,
  markRenderComplete,
  markGeneRendered,
} from "@/lib/utils/test-hooks";
import { SpatialScaleBar } from "@/components/spatial-scale-bar";
import { VISUALIZATION_CONFIG } from "@/lib/config/visualization.config";
import {
  applyTransformsToPositions,
  computeSampleCentroids,
} from "@/lib/utils/sample-transforms";
import { useTourPlayer } from "@/lib/hooks/useTourPlayer";
import { useTourStore } from "@/lib/stores/tourStore";
import { openImageRef, type SpatialDataImage } from "@/lib/spatialdata/image";

export type { SceneHandle } from "@/lib/scene/contract";

// On-demand rendering: a frame is drawn when the controls move, a layer
// changes or an effect here mutates what is drawn (invalidateRef); this is
// the safety-net rate for anything that forgets to say so.
const IDLE_FPS = 10;

interface ThreeSceneProps {
  dataset?: StandardizedDataset | null;
  /**
   * What is drawn with the cells, opened by the viewer (a SpatialData
   * store's image, polygons and points). Molecules given here take the
   * place of any linked by a mapping.json.
   */
  elements?: OpenedElements;
}

export function ThreeScene({ dataset, elements }: ThreeSceneProps) {
  const molecules = elements?.points ?? elements?.molecules ?? null;
  const moleculesRef = useRef(molecules);

  moleculesRef.current = molecules;

  const sceneHandleRef = useRef<SceneHandle | null>(null);
  // The scene's side of the layer contract: layers join and leave through
  // it, and it ticks them from the scene loop (lib/scene/LayerHost.ts).
  const hostRef = useRef<LayerHost | null>(null);
  // Ask the scene loop for a frame (on-demand rendering): called after every
  // mutation of what is drawn. Layers do it themselves through `changed()`.
  const invalidateRef = useRef<() => void>(() => {});
  // Image named by the dataset's mapping.json (`image`) or dropped beside
  // it, drawn under the cells; one given in `elements` takes precedence.
  const [overlayImage, setOverlayImage] = useState<SpatialDataImage | null>(
    null,
  );
  // Cell polygons dropped with the cells for a preview
  const [overlayShapes, setOverlayShapes] = useState<
    OpenedElements["shapes"] | null
  >(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const pointCloudRef = useRef<THREE.Points | null>(null);
  const [pointCloudVersion, setPointCloudVersion] = useState(0);
  const [webglFailed, setWebglFailed] = useState(false);
  const sceneRef = useRef<THREE.Scene | null>(null);
  // Outer group sits at the data center and carries rotation/flip so the
  // dataset spins in place; inner group is offset by -center and holds the
  // point clouds (SC + SM overlay).
  const sceneGroupRef = useRef<THREE.Group | null>(null);
  const innerGroupRef = useRef<THREE.Group | null>(null);
  // Untransformed scene-space positions (xyz interleaved). Used as the base
  // when applying per-sample transforms so edits compose from identity.
  const basePositionsRef = useRef<Float32Array | null>(null);
  const sceneScaleRef = useRef<number>(1);
  const geneToastIdRef = useRef<string | number | null>(null);

  // Raycaster and interaction state
  const raycasterRef = useRef<THREE.Raycaster>(new THREE.Raycaster());
  // World-space lookAt captured at scene init — restored by the camera reset.
  const initialLookAtRef = useRef<THREE.Vector3 | null>(null);
  // gsap.ticker callback for the continuous O-orbit. Set while orbiting.
  const orbitTickerRef = useRef<
    ((time: number, deltaTime: number) => void) | null
  >(null);
  // Mirror orbit settings into refs so the ticker reads current values
  // without having to restart the orbit on every slider drag.
  const orbitSpeedRef = useRef(1.0);
  const orbitYBobRef = useRef(0.3);
  // Distance-from-target filter state — pushed to shader uniforms each frame.
  // The radius is in absolute world units (microns for raw-coord datasets).
  const targetFilterEnabledRef = useRef(false);
  const targetFilterRadiusRef = useRef(5000);
  const mouseRef = useRef<THREE.Vector2>(new THREE.Vector2());
  const cameraRef = useRef<THREE.Camera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const controlsRef = useRef<any>(null);
  const tooltipRef = useRef<HTMLDivElement | null>(null);
  const hoveredPointRef = useRef<number | null>(null);
  const lastCameraPositionRef = useRef<THREE.Vector3>(new THREE.Vector3());

  // Store current visualization data for tooltips
  const geneExpressionRef = useRef<number[] | null>(null);
  // Tracks the last applied raw/log display mode so a toggle re-auto-scales the
  // gene range (like a gene change), since counts and log1p have different spans.
  const previousModeRef = useRef<string | null>(null);
  const colorPaletteRef = useRef<Record<string, string>>({});
  const clusterRef = useRef<any>(null);
  const isNumericalClusterRef = useRef<boolean>(false);
  const pinnedClustersRef = useRef<
    Array<{
      column: string;
      cluster: any;
      palette: Record<string, string>;
      isNumerical: boolean;
    }>
  >([]);
  const baseDotSizeRef = useRef<number>(5);
  const cameraDistanceRef = useRef<number>(1);

  // SM overlay state (for __all__ mapping). Gene selection lives in the
  // global SM viz store so the picker + legends drive what's rendered.
  const smDatasetRef = useRef<SingleMoleculeDataset | null>(null);
  const [smDataset, setSmDataset] = useState<SingleMoleculeDataset | null>(
    null,
  );
  const [smLoading, setSmLoading] = useState(false);
  // Last cmd/ctrl+click world coordinates; rendered as a small overlay chip
  // so the user can read them directly without opening devtools.
  const [lastClickCoords, setLastClickCoords] = useState<{
    x: number;
    y: number;
    z: number;
  } | null>(null);
  // The overlay's stores are the panel's own, like the gene picker and the
  // legends that drive them — so each side of a split has its own overlay.
  const smStoreApi = usePanelSingleMoleculeApi();
  const smVizApi = usePanelSingleMoleculeVisualizationApi();

  // What the scene draws besides the cells, as layers of the contract: the
  // elements the viewer opened, and the overlays this scene resolves itself
  // (a mapping.json, a folder dropped for a preview). The layers outlive
  // scene rebuilds; `onSceneReady` attaches them to each new scene.
  const handedImage = elements?.image;
  const handedShapes = elements?.shapes;
  const layerElements = useMemo<OpenedElements>(
    () => ({
      image: handedImage ?? overlayImage ?? undefined,
      shapes: handedShapes ?? overlayShapes ?? undefined,
      ...moleculeElements(smDataset),
    }),
    [handedImage, handedShapes, overlayImage, overlayShapes, smDataset],
  );
  const layerCtx: LayerContext = {
    molecules: {
      // Closing a gene's loading toast cancels its download; take it out of
      // the selection so nothing restarts it.
      onCancelGene: (gene: string) => smVizApi.getState().removeGene(gene),
    },
  };
  const { onSceneReady, layers } = useLayers(layerElements, layerCtx);
  // The overlay's molecule layer, for the per-frame uniform push
  const smLayerRef = useRef<MoleculeLayer | null>(null);

  smLayerRef.current = (layers.molecules as MoleculeLayer | undefined) ?? null;

  // Get visualization settings from store
  const {
    mode,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    geneEverywhere,
    hiddenCelltypes,
    colorPalette,
    alphaScale,
    sizeScale,
    geneScaleMin,
    geneScaleMax,
    setGeneScaleMin,
    setGeneScaleMax,
    geneScaleMode,
    geneScaleModeUserSet,
    setGeneScaleModeAuto,
    setDetectedGeneKind,
    numericalScaleMin,
    numericalScaleMax,
    setNumericalScaleMin,
    setNumericalScaleMax,
    toggleCelltype,
    clusterVersion,
    incrementClusterVersion,
    cameraResetSignal,
    orbitSpeed,
    orbitYBob,
    targetFilterEnabled,
    targetFilterRadius,
    scLayerVisible,
    columnTypeOverrides,
    viewMode,
    targetPx,
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    sceneRotation,
    flipX,
    flipY,
    pinnedTooltipColumns,
    colormap,
    transformColumn,
    sampleTransforms,
    transformVersion,
    coexpressEnabled,
    selectedGene2,
    coexpressSwapped,
    gene2ScaleMin,
    gene2ScaleMax,
    setGene2ScaleMin,
    setGene2ScaleMax,
    exportBoxEnabled,
    exportBoxWidthMm,
    exportBoxHeightMm,
    exportBoxCenterPx,
    setExportBoxCenterPx,
  } = usePanelVisualizationStore();

  // Mirror orbit settings into refs (read by gsap.ticker each frame so the
  // sliders update a running orbit live without restarting it).
  useEffect(() => {
    orbitSpeedRef.current = orbitSpeed;
  }, [orbitSpeed]);
  useEffect(() => {
    orbitYBobRef.current = orbitYBob;
  }, [orbitYBob]);
  useEffect(() => {
    targetFilterEnabledRef.current = targetFilterEnabled;
  }, [targetFilterEnabled]);
  useEffect(() => {
    targetFilterRadiusRef.current = targetFilterRadius;
  }, [targetFilterRadius]);

  // Master show/hide for the single-cell point cloud. Mirror into a ref so the
  // point-cloud build effect can stamp visibility onto a freshly created cloud,
  // and apply live whenever the flag toggles.
  const scLayerVisibleRef = useRef(scLayerVisible);

  useEffect(() => {
    scLayerVisibleRef.current = scLayerVisible;
    if (pointCloudRef.current) {
      pointCloudRef.current.visible = scLayerVisible;
      invalidateRef.current();
    }
  }, [scLayerVisible]);

  // Tour player — drives camera + SM gene selection from the loaded tour JSON.
  // No-op until the user uploads a tour file and presses Play.
  const setTargetFilterEnabledFn = usePanelVisualizationStore(
    (s) => s.setTargetFilterEnabled,
  );
  const setTargetFilterRadiusFn = usePanelVisualizationStore(
    (s) => s.setTargetFilterRadius,
  );

  useTourPlayer({
    cameraRef,
    controlsRef,
    orbitTickerRef,
    isNormalized: dataset?.normalized !== false,
    targetFilterEnabledRef,
    targetFilterRadiusRef,
    setTargetFilterEnabled: setTargetFilterEnabledFn,
    setTargetFilterRadius: setTargetFilterRadiusFn,
    scPointCloudRef: pointCloudRef,
  });

  // Split screen support
  const panelId = usePanelId();
  const { enableSplit, setRightPanelS3 } = useSplitScreenStore();
  const linkConfigRef = useRef(dataset ? getDatasetLinkConfig(dataset) : null);

  // Update link config when dataset changes:
  // 1. Check hardcoded registry first (instant)
  // 2. Try fetching mapping.json from custom S3 (async, background)
  // 3. If mapping found, lazy-load the link column so right-click works immediately
  useEffect(() => {
    // Drop any SM overlay carried over from a previously-viewed dataset.
    // Only an "__all__" mapping (handled below) re-establishes one, so a
    // non-"__all__" (multi-slice) or unmapped dataset never inherits a
    // stale overlay. The right-click split view for non-"__all__" mappings
    // is unaffected — it reads linkConfigRef, not this overlay state.
    const clearSmOverlay = () => {
      const prev = smDatasetRef.current;

      if (!prev) return;
      smDatasetRef.current = null;
      setSmDataset(null);
      // Also wipe the SM viz selection so the gene picker + legends (gated on
      // selectedGenesLegend, not on the dataset) disappear with the overlay.
      smVizApi.getState().clearGenes();
      smStoreApi.getState().removeDataset(prev.id);
    };

    // Load the single shared SM dataset for an "__all__" mapping and push it
    // into the overlay + this panel's SM store. `source` decides how it's loaded:
    //   - "s3": smRef is a PUBLIC S3 base URL → fromCustomS3()
    //   - "app": smRef is an app dataset id (sm_…) → fromS3() (presigned API),
    //            which is how app-uploaded overlays (no public bucket) load.
    const loadAllOverlay = async (smRef: string, source: "app" | "s3") => {
      if (!smRef) return;
      setSmLoading(true);
      try {
        const { SingleMoleculeDataset: SMDataset } = await import(
          "@/lib/SingleMoleculeDataset"
        );
        const smDs =
          source === "app"
            ? await SMDataset.fromS3(smRef)
            : await SMDataset.fromCustomS3(smRef);

        // This effect runs twice under StrictMode and on every dataset
        // change; only the latest run may install its overlay, or two layers
        // would be built and every gene downloaded twice.
        if (stale) return;
        smDatasetRef.current = smDs;
        setSmDataset(smDs);

        smStoreApi.getState().addDataset(smDs);
        // The saved-default-genes fallback keys off a public URL, so it's only
        // meaningful for the "s3" source — and only the main (left) viewer's
        // page reads it, from the global store.
        if (source === "s3" && smStoreApi === useSingleMoleculeStore) {
          useSingleMoleculeStore.getState().setOverlaySourceUrl(smRef);
        }
      } catch (error) {
        console.warn("Failed to load SM overlay dataset:", error);
      } finally {
        setSmLoading(false);
      }
    };

    clearSmOverlay();
    setOverlayImage(null);
    setOverlayShapes(null);

    if (!dataset) {
      linkConfigRef.current = null;

      return;
    }

    let stale = false;
    const loadOverlayImage = async (
      mapping?: Pick<DatasetLinkConfig, "image" | "imageSource"> | null,
    ) => {
      if (!mapping?.image) return;
      // The image is placed in raw coordinates; a normalized scene has none.
      if (dataset.normalized !== false) {
        console.warn("mapping.json image skipped: dataset is normalized");

        return;
      }
      try {
        const image = await openImageRef({
          image: mapping.image,
          source: mapping.imageSource,
        });

        if (stale) return;
        setOverlayImage(image);
      } catch (error) {
        console.warn("Failed to load overlay image:", error);
      }
    };

    // Molecules handed over with the cells take the place of a linked set
    if (moleculesRef.current) {
      smDatasetRef.current = moleculesRef.current;
      setSmDataset(moleculesRef.current);
      smStoreApi.getState().addDataset(moleculesRef.current);

      return () => {
        stale = true;
      };
    }

    // A folder dropped for a local preview brings its overlays with it,
    // already open: there is nothing to look up or fetch.
    const local = useLocalBundleStore.getState();

    if (local.anchorId === dataset.id) {
      if (local.image && dataset.normalized !== false) {
        toast.warn(
          "The image is not drawn: these cells were saved with normalized coordinates, and an image is placed in the original ones.",
          { autoClose: 12000, toastId: "local-image-normalized" },
        );
      } else if (local.image) {
        setOverlayImage(local.image);
      }
      if (local.shapes) setOverlayShapes(() => local.shapes);
      if (local.molecules) {
        smDatasetRef.current = local.molecules;
        setSmDataset(local.molecules);
        smStoreApi.getState().addDataset(local.molecules);
      }

      return () => {
        stale = true;
      };
    }

    // Check hardcoded registry first
    const registryConfig = getDatasetLinkConfig(dataset);

    linkConfigRef.current = registryConfig;

    // Try fetching mapping.json (only for custom S3 datasets)
    if (!registryConfig && dataset.metadata?.customS3BaseUrl) {
      fetchMappingConfig(dataset).then(async (mappingConfig) => {
        if (!mappingConfig) return;

        void loadOverlayImage(mappingConfig);
        // A mapping that only names an image links no molecules
        if (!mappingConfig.links) return;
        linkConfigRef.current = mappingConfig;

        // "__all__" means every cell links to the same SM dataset — load overlay
        if (mappingConfig.linkColumn === "__all__") {
          await loadAllOverlay(
            mappingConfig.links["__all__"],
            mappingConfig.source ?? "s3",
          );

          return;
        }

        // Check if link column is already loaded
        const alreadyLoaded = dataset.clusters?.some(
          (c) => c.column === mappingConfig.linkColumn,
        );

        if (alreadyLoaded) return;

        // Check if the column exists in the dataset
        if (
          dataset.allClusterColumnNames &&
          !dataset.allClusterColumnNames.includes(mappingConfig.linkColumn)
        ) {
          console.warn(
            `mapping.json linkColumn "${mappingConfig.linkColumn}" not found in dataset columns`,
          );

          return;
        }

        // Lazy-load the link column in the background
        try {
          const { getStandardizedDatasetWorker } = await import(
            "@/lib/workers/standardizedDatasetWorkerManager"
          );
          const worker = await getStandardizedDatasetWorker();

          const newClusters = await worker.loadClusterFromS3(
            dataset.id,
            [mappingConfig.linkColumn],
            dataset.metadata?.customS3BaseUrl,
          );

          if (newClusters && newClusters.length > 0) {
            dataset.addClusters(newClusters);
            incrementClusterVersion();
          }
        } catch (error) {
          console.warn(
            `Failed to preload link column "${mappingConfig.linkColumn}":`,
            error,
          );
        }
      });
    }

    // App-uploaded datasets open by id (no public custom S3 URL), so the fetch
    // above skips them. Try our own mapping.json proxy instead — currently only
    // the "__all__" overlay is produced this way (the combined SC+SM upload).
    if (!registryConfig && !dataset.metadata?.customS3BaseUrl && dataset.id) {
      // Revalidate: the route is cached for a minute, which would hide an
      // overlay that was just attached or removed until the cache expired.
      fetch(`/api/datasets/${dataset.id}/mapping`, { cache: "no-cache" })
        .then((r) => (r.ok ? r.json() : null))
        .then((mc: DatasetLinkConfig | null) => {
          void loadOverlayImage(mc);
          if (mc && mc.linkColumn === "__all__" && mc.links) {
            linkConfigRef.current = mc;
            void loadAllOverlay(mc.links["__all__"], mc.source ?? "s3");
          }
        })
        .catch(() => {});
    }

    return () => {
      stale = true;
    };
  }, [dataset, incrementClusterVersion]);

  // Store current mode and selection in refs to avoid closure issues
  const modeRef = useRef(mode);
  const selectedGeneRef = useRef(selectedGene);
  const selectedColumnRef = useRef(selectedColumn);
  const previousGeneRef = useRef<string | null>(null);
  const previousGene2Ref = useRef<string | null>(null);
  const previousColumnRef = useRef<string | null>(null);
  const previousColumnTypeRef = useRef<boolean>(false);

  // Update refs when store values change
  useEffect(() => {
    modeRef.current = mode;
    selectedGeneRef.current = selectedGene;
    selectedColumnRef.current = selectedColumn;
  }, [mode, selectedGene, selectedColumn]);

  // Helper function: Create tooltip element
  const createTooltip = (): HTMLDivElement => {
    const tooltip = document.createElement("div");

    tooltip.className =
      "absolute bg-black/80 text-white px-2.5 py-1.5 rounded text-sm font-sans pointer-events-none hidden z-[1000] shadow-lg min-w-[80px]";
    document.body.appendChild(tooltip);

    return tooltip;
  };

  // Helper function: Get color from point cloud at specific index
  const getPointColor = (index: number): string => {
    if (!pointCloudRef.current) return "#808080";
    const colorAttr = pointCloudRef.current.geometry.attributes
      .color as THREE.BufferAttribute;

    if (!colorAttr) return "#808080";

    const r = Math.round(colorAttr.getX(index) * 255);
    const g = Math.round(colorAttr.getY(index) * 255);
    const b = Math.round(colorAttr.getZ(index) * 255);

    return `#${r.toString(16).padStart(2, "0")}${g.toString(16).padStart(2, "0")}${b.toString(16).padStart(2, "0")}`;
  };

  // Helper function: Show tooltip
  const showTooltip = (position: THREE.Vector3, index: number) => {
    if (!tooltipRef.current || !cameraRef.current || !rendererRef.current)
      return;

    // Convert 3D position to screen coordinates
    const vector = position.clone();

    vector.project(cameraRef.current);

    // Get canvas position relative to screen coordinates
    const rect = rendererRef.current.domElement.getBoundingClientRect();

    const x = (vector.x * 0.5 + 0.5) * rect.width + rect.left;
    const y = (-vector.y * 0.5 + 0.5) * rect.height + rect.top;

    // Get point color (gene gradient color in gene mode)
    const pointColor = getPointColor(index);

    // Determine what to show based on current mode and cluster type
    const isNumerical = isNumericalClusterRef.current;
    const clusterValue = clusterRef.current
      ? getClusterValue(clusterRef.current, index)
      : undefined;
    const geneValue = geneExpressionRef.current?.[index];

    // Get current values from refs
    const currentMode = modeRef.current;
    const currentGene = selectedGeneRef.current;
    const currentColumn = selectedColumnRef.current;

    // Build pinned column rows
    const pinnedRows = pinnedClustersRef.current
      .map((p) => {
        const val = getClusterValue(p.cluster, index);

        if (p.isNumerical) {
          return `<div class="text-xs text-default-400">${p.column}: ${val}</div>`;
        }
        const color = p.palette[String(val)] || "#808080";

        return `<div class="flex items-center text-xs">
        <div style="width: 10px; height: 10px; border-radius: 50%; background-color: ${color}; margin-right: 5px;"></div>
        <span class="text-default-400">${p.column}: ${val}</span>
      </div>`;
      })
      .join("");

    let tooltipContent = "";

    if (currentMode.includes("gene") && currentGene) {
      // Gene mode
      if (isNumerical) {
        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div>${currentColumn}: ${clusterValue}</div>
            ${pinnedRows}
            <div>${currentGene}: ${geneValue?.toFixed(2) ?? "N/A"}</div>
          </div>
        `;
      } else {
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";

        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
              <span>${clusterValue}</span>
            </div>
            ${pinnedRows}
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${pointColor}; margin-right: 6px;"></div>
              <span>${geneValue?.toFixed(2) ?? "N/A"}</span>
            </div>
          </div>
        `;
      }
    } else {
      // Celltype mode or gene mode without gene selected
      if (isNumerical) {
        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div>${currentColumn}: ${clusterValue}</div>
            ${pinnedRows}
          </div>
        `;
      } else {
        const clusterColor =
          colorPaletteRef.current[String(clusterValue)] || "#808080";

        tooltipContent = `
          <div class="flex flex-col gap-1">
            <div class="flex items-center">
              <div style="width: 12px; height: 12px; border-radius: 50%; background-color: ${clusterColor}; margin-right: 6px;"></div>
              <span>${clusterValue}</span>
            </div>
            ${pinnedRows}
          </div>
        `;
      }
    }

    tooltipRef.current.innerHTML = tooltipContent;
    tooltipRef.current.style.left = `${x + 10}px`;
    tooltipRef.current.style.top = `${y + 10}px`;
    tooltipRef.current.classList.remove("hidden");
  };

  // Helper function: Hide tooltip
  const hideTooltip = () => {
    if (tooltipRef.current) {
      tooltipRef.current.classList.add("hidden");
    }
  };

  // Helper function: Check intersections with adaptive threshold
  const checkIntersections = () => {
    if (!pointCloudRef.current || !cameraRef.current || !rendererRef.current)
      return;

    // Get current camera position to check if we've moved
    const currentCameraPosition = cameraRef.current.position.clone();
    const cameraHasMoved = !currentCameraPosition.equals(
      lastCameraPositionRef.current,
    );

    lastCameraPositionRef.current.copy(currentCameraPosition);

    // Screen-space pixel threshold: convert a fixed pixel radius to world units
    // This adapts automatically to any coordinate range and zoom level
    const PIXEL_RADIUS = 5; // hover within 5 pixels of a point
    const camera = cameraRef.current as THREE.PerspectiveCamera;
    const target = controlsRef.current?.target;
    const cameraDistance = target
      ? camera.position.distanceTo(target)
      : camera.position.length();
    const canvasHeight = rendererRef.current.domElement.clientHeight;
    const fovRad = (camera.fov / 2) * (Math.PI / 180);
    const pixelSize = (2 * cameraDistance * Math.tan(fovRad)) / canvasHeight;
    const threshold = pixelSize * PIXEL_RADIUS;

    raycasterRef.current.params.Points!.threshold = threshold;

    // Update the raycaster with the current mouse position and camera
    raycasterRef.current.setFromCamera(mouseRef.current, cameraRef.current);

    // Check for intersections with the points mesh
    const intersects = raycasterRef.current.intersectObject(
      pointCloudRef.current,
    );

    // If we found an intersection
    if (intersects.length > 0) {
      // Sort intersections by distance if there are multiple
      if (intersects.length > 1) {
        intersects.sort((a, b) => a.distance - b.distance);
      }

      // Get the index of the closest point that was intersected
      const index = intersects[0].index!;

      // If this is a different point than the one we were previously hovering over
      if (hoveredPointRef.current !== index) {
        hoveredPointRef.current = index;
        const position = intersects[0].point;

        showTooltip(position, index);
      }
    } else {
      // If we're not hovering over any point, hide the tooltip
      if (hoveredPointRef.current !== null) {
        hoveredPointRef.current = null;
        hideTooltip();
      }
    }
  };

  // Effect 0: Tooltip cleanup on component unmount
  useEffect(() => {
    return () => {
      if (tooltipRef.current) {
        document.body.removeChild(tooltipRef.current);
        tooltipRef.current = null;
      }
    };
  }, []);

  // Effect 1: Scene initialization - runs when dataset changes
  useEffect(() => {
    if (!containerRef.current) return;

    // If dataset is provided, use its spatial coordinates
    if (dataset) {
      // Calculate bounding box and center of spatial data
      const bounds = {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity,
      };

      const coords = dataset.spatial.coordinates;
      const dims = dataset.spatial.dimensions;
      const numPts = dataset.getPointCount();

      for (let p = 0; p < numPts; p++) {
        const x =
          coords instanceof Float32Array
            ? coords[p * dims]
            : (coords as number[][])[p][0];
        const y =
          coords instanceof Float32Array
            ? coords[p * dims + 1]
            : (coords as number[][])[p][1];

        bounds.minX = Math.min(bounds.minX, x);
        bounds.maxX = Math.max(bounds.maxX, x);
        bounds.minY = Math.min(bounds.minY, y);
        bounds.maxY = Math.max(bounds.maxY, y);
        if (dims === 3) {
          const z =
            coords instanceof Float32Array
              ? coords[p * dims + 2]
              : (coords as number[][])[p][2];

          if (z !== undefined) {
            bounds.minZ = Math.min(bounds.minZ, z);
            bounds.maxZ = Math.max(bounds.maxZ, z);
          }
        }
      }

      // Scale factor: 500x for normalized [-1,1] data, 1x for raw coordinates
      const coordScale = dataset.normalized === false ? 1 : 500;

      // Calculate center point (scaled)
      const center = new THREE.Vector3(
        ((bounds.minX + bounds.maxX) / 2) * coordScale,
        ((bounds.minY + bounds.maxY) / 2) * coordScale,
        dataset.spatial.dimensions === 3
          ? ((bounds.minZ + bounds.maxZ) / 2) * coordScale
          : 0,
      );

      // Calculate size of data
      const size = Math.max(
        (bounds.maxX - bounds.minX) * coordScale,
        (bounds.maxY - bounds.minY) * coordScale,
        dataset.spatial.dimensions === 3
          ? (bounds.maxZ - bounds.minZ) * coordScale
          : 0,
      );

      // Position camera at appropriate distance
      const distance = size * 1.5;
      const cameraPos = new THREE.Vector3(
        center.x,
        center.y,
        center.z + distance,
      );

      // Set near/far planes based on data scale to avoid clipping
      const maxExtent = Math.max(size, distance);
      const near = maxExtent * 0.001;
      const far = maxExtent * 10;

      // Initialize Three.js scene with options. WebGL context creation can
      // fail (GPU blocklist, stale browser update, software GL) — surface
      // that instead of leaving a silently blank canvas.
      let init: ReturnType<typeof initializeScene> | null = null;

      try {
        init = initializeScene(containerRef.current, {
          is2D: viewMode === "2D",
          cameraPosition: cameraPos,
          lookAtPosition: center,
          near,
          far,
          onFrame: (now) => hostRef.current?.tick(now),
          idleFps: IDLE_FPS,
        });
        setWebglFailed(false);
      } catch (e) {
        console.error("[ThreeScene] WebGL init failed:", e);
        setWebglFailed(true);
      }
      if (!init) return;
      const {
        scene,
        camera,
        renderer,
        controls,
        animate,
        dispose,
        invalidate,
      } = init;

      invalidateRef.current = invalidate;

      // Store camera, renderer, and controls refs for raycasting + scale bar
      cameraRef.current = camera;
      rendererRef.current = renderer;
      controlsRef.current = controls;

      // Stash the initial lookAt so the camera reset can snap back to it.
      initialLookAtRef.current = center.clone();

      // Create tooltip
      if (!tooltipRef.current) {
        tooltipRef.current = createTooltip();
      }

      // Throttled intersection checking (50ms)
      let lastCheckTime = 0;
      const throttleDelay = 10;

      // Mouse move handler
      const handleMouseMove = (event: MouseEvent) => {
        const rect = renderer.domElement.getBoundingClientRect();

        mouseRef.current.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
        mouseRef.current.y =
          -((event.clientY - rect.top) / rect.height) * 2 + 1;

        // Throttle intersection checking
        const now = Date.now();

        if (now - lastCheckTime > throttleDelay) {
          lastCheckTime = now;
          checkIntersections();
        }
      };

      // Double-click handler
      const handleDoubleClick = () => {
        if (hoveredPointRef.current !== null && dataset.clusters) {
          const index = hoveredPointRef.current;
          const clusterValueStr = clusterRef.current
            ? getClusterValue(clusterRef.current, index)
            : "";

          // Only toggle celltype if it's not a numerical cluster
          if (!isNumericalClusterRef.current) {
            // Toggle the cluster in selectedCelltypes
            toggleCelltype(clusterValueStr);
            // Mode is now automatically updated by toggleCelltype
          }
        }
      };

      // Right-click handler: open split screen with linked SM dataset
      const handleContextMenu = (event: MouseEvent) => {
        // Only works on left panel (not inside a split right panel)
        if (panelId) return;
        // Need a hovered point
        if (hoveredPointRef.current === null) return;
        // Need a link config for this dataset
        const linkConfig = linkConfigRef.current;

        if (!linkConfig) return;

        event.preventDefault();

        const index = hoveredPointRef.current;

        // "__all__" linkColumn means every cell maps to the same SM dataset,
        // which is already rendered as an overlay on this scene — no need to
        // open a split panel on right-click.
        if (linkConfig.linkColumn === "__all__") return;

        // Find the link column in dataset clusters (always reads the configured column,
        // regardless of which column is currently selected for visualization)
        const linkCluster = dataset!.clusters?.find(
          (c) => c.column === linkConfig.linkColumn,
        );

        if (!linkCluster) {
          console.warn(
            `Link column "${linkConfig.linkColumn}" not found in dataset clusters`,
          );

          return;
        }

        const setValue = getClusterValue(linkCluster, index);
        const smUrl = linkConfig.links[setValue];

        if (!smUrl) {
          toast.warning(`No single molecule data available for "${setValue}"`);

          return;
        }

        enableSplit();
        setRightPanelS3(smUrl, "sm");
        toast.info(`Opening SM dataset for ${setValue}`);
      };

      // Cmd/Ctrl+click: pan the camera (and its orbit/pan pivot) so the
      // clicked SC cell lands at the screen center. Translation only —
      // camera↔target distance and angle are preserved. Works in both 2D
      // (OrbitControls, top-down) and 3D (TrackballControls) because both
      // expose the same .target/.position API. ~350ms gsap tween; a second
      // click during the tween overrides cleanly.
      const handleClickRecenter = (event: MouseEvent) => {
        if (!(event.metaKey || event.ctrlKey)) return;
        const idx = hoveredPointRef.current;
        const pc = pointCloudRef.current;
        const ctrls = controlsRef.current;
        const cam = cameraRef.current;

        if (idx == null || !pc || !ctrls || !cam) return;

        // Stop a running orbit before panning — otherwise the orbit ticker
        // keeps overwriting the camera position during the pan tween.
        if (orbitTickerRef.current) {
          gsap.ticker.remove(orbitTickerRef.current);
          orbitTickerRef.current = null;
        }

        const posAttr = pc.geometry.attributes
          .position as THREE.BufferAttribute;
        const local = new THREE.Vector3(
          posAttr.getX(idx),
          posAttr.getY(idx),
          posAttr.getZ(idx),
        );

        pc.updateMatrixWorld();
        const world = local.applyMatrix4(pc.matrixWorld);

        setLastClickCoords({ x: world.x, y: world.y, z: world.z });

        const targetStart = ctrls.target.clone();
        const positionStart = cam.position.clone();
        const delta = world.clone().sub(targetStart);
        const positionEnd = positionStart.clone().add(delta);
        const tweenObj = { t: 0 };

        gsap.to(tweenObj, {
          t: 1,
          duration: 0.35,
          ease: "power2.out",
          overwrite: true,
          onUpdate: () => {
            ctrls.target.lerpVectors(targetStart, world, tweenObj.t);
            cam.position.lerpVectors(positionStart, positionEnd, tweenObj.t);
            ctrls.update();
          },
        });
      };

      // Add event listeners
      renderer.domElement.addEventListener("mousemove", handleMouseMove);
      renderer.domElement.addEventListener("dblclick", handleDoubleClick);
      renderer.domElement.addEventListener("contextmenu", handleContextMenu);
      renderer.domElement.addEventListener("click", handleClickRecenter);

      // Build positions Float32Array directly from dataset coordinates
      const spatialCoords = dataset.spatial.coordinates;
      const spatialDims = dataset.spatial.dimensions;
      const ptCount = dataset.getPointCount();
      const positions = new Float32Array(ptCount * 3);

      // Scale factor: 500x for normalized [-1,1] data, 1x for raw coordinates
      const cs = dataset.normalized === false ? 1 : 500;

      if (spatialCoords instanceof Float32Array) {
        // Flat Float32Array path (optimized — no object creation)
        if (spatialDims === 3 && cs === 1) {
          // Best case: 3D raw coords — direct copy, no scaling needed
          positions.set(spatialCoords);
        } else {
          for (let p = 0; p < ptCount; p++) {
            positions[p * 3] = spatialCoords[p * spatialDims] * cs;
            positions[p * 3 + 1] = spatialCoords[p * spatialDims + 1] * cs;
            positions[p * 3 + 2] =
              spatialDims === 3 ? spatialCoords[p * spatialDims + 2] * cs : 0;
          }
        }
      } else {
        // number[][] path
        const coords = spatialCoords as number[][];

        for (let p = 0; p < ptCount; p++) {
          positions[p * 3] = coords[p][0] * cs;
          positions[p * 3 + 1] = coords[p][1] * cs;
          positions[p * 3 + 2] =
            spatialDims === 3 ? (coords[p][2] ?? 0) * cs : 0;
        }
      }

      // Stash base (untransformed) positions and the scene-space scale so the
      // per-sample transform effect can recompute live without re-deriving.
      basePositionsRef.current = new Float32Array(positions);
      sceneScaleRef.current = cs;

      // Create point cloud mesh with custom shaders
      // The shader computes: gl_PointSize = size * dotSize * proj[1][1] / -mvPosition.z
      // We want ~3px dots at the initial zoom level.
      // Back-calculate: dotSize = targetPx * distance / (baseSize * proj[1][1])
      // where distance ≈ size*1.5, proj[1][1] ≈ 1.3 (75° FOV), baseSize = POINT_BASE_SIZE
      const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2); // ~1.3

      cameraDistanceRef.current = distance;
      const baseDotSize =
        (targetPx * distance) / (VISUALIZATION_CONFIG.POINT_BASE_SIZE * proj11);

      baseDotSizeRef.current = baseDotSize;
      const pointCloud = createPointCloudFromBuffers(
        positions,
        ptCount,
        baseDotSize,
      );

      pointCloudRef.current = pointCloud; // Store reference
      pointCloud.visible = scLayerVisibleRef.current; // honor master SC toggle
      sceneRef.current = scene; // Store scene reference

      // Wrap point clouds in an outer/inner group pair so scene transforms
      // (rotation, flip) pivot on the data center instead of the world
      // origin. Outer sits at center and carries the transform; inner is
      // offset by -center and holds the clouds, so a point at absolute
      // coord P still renders at P (center + -center + P).
      const outerGroup = new THREE.Group();
      const innerGroup = new THREE.Group();

      outerGroup.position.copy(center);
      innerGroup.position.set(-center.x, -center.y, -center.z);
      outerGroup.add(innerGroup);
      innerGroup.add(pointCloud);
      scene.add(outerGroup);
      sceneGroupRef.current = outerGroup;
      innerGroupRef.current = innerGroup;
      setPointCloudVersion((v) => v + 1); // Trigger visualization update

      // Start animation
      animate();
      markSceneReady();
      const host = new LayerHost(
        innerGroup,
        camera,
        renderer,
        invalidate,
        () =>
          (pointCloudRef.current?.geometry.getAttribute("color") as
            | THREE.BufferAttribute
            | undefined) ?? null,
      );
      const sceneHandle = host.handle;

      hostRef.current = host;
      const removeLayers = onSceneReady(sceneHandle);

      sceneHandleRef.current = sceneHandle;

      // Cleanup on unmount
      return () => {
        removeLayers();
        host.dispose();
        hostRef.current = null;
        invalidateRef.current = () => {};
        sceneHandleRef.current = null;
        // Remove event listeners
        renderer.domElement.removeEventListener("mousemove", handleMouseMove);
        renderer.domElement.removeEventListener("dblclick", handleDoubleClick);
        renderer.domElement.removeEventListener(
          "contextmenu",
          handleContextMenu,
        );
        renderer.domElement.removeEventListener("click", handleClickRecenter);

        // Hide tooltip
        hideTooltip();

        if (sceneGroupRef.current) {
          scene.remove(sceneGroupRef.current);
          sceneGroupRef.current = null;
        }
        innerGroupRef.current = null;
        pointCloud.geometry.dispose();
        (pointCloud.material as any).dispose();
        pointCloudRef.current = null;
        sceneRef.current = null;
        basePositionsRef.current = null;
        cameraRef.current = null;
        rendererRef.current = null;
        controlsRef.current = null;
        dispose();
      };
    } else {
      // Initialize scene without dataset
      let init: ReturnType<typeof initializeScene> | null = null;

      try {
        init = initializeScene(containerRef.current, { idleFps: IDLE_FPS });
        setWebglFailed(false);
      } catch (e) {
        console.error("[ThreeScene] WebGL init failed:", e);
        setWebglFailed(true);
      }
      if (!init) return;
      const { animate, dispose } = init;

      // Start animation even without data
      animate();

      return () => {
        dispose();
      };
    }
  }, [dataset, viewMode]);

  // Camera reset: snap controls.target back to the initial lookAt captured at
  // scene init. Triggered by cameraResetSignal (Reset View button) or R key.
  useEffect(() => {
    if (cameraResetSignal === 0) return;
    const ctrls = controlsRef.current;
    const initial = initialLookAtRef.current;

    if (!ctrls || !initial) return;
    ctrls.target.copy(initial);
    ctrls.update();
  }, [cameraResetSignal]);

  // R key resets the camera (skipped while typing in an input).
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "r" && e.key !== "R") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      const ctrls = controlsRef.current;
      const initial = initialLookAtRef.current;

      if (!ctrls || !initial) return;
      ctrls.target.copy(initial);
      ctrls.update();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, []);

  // O key toggles a continuous orbit around controls.target in 3D mode.
  // The camera rotates around world Y at ~4.5s per revolution and also bobs
  // up/down (Y offset) on a slower phase, so the view sweeps the cell from
  // multiple angles rather than tracing a flat ring. Pressing O again stops.
  useEffect(() => {
    const stopOrbit = () => {
      if (orbitTickerRef.current) {
        gsap.ticker.remove(orbitTickerRef.current);
        orbitTickerRef.current = null;
      }
    };

    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "o" && e.key !== "O") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      if (viewMode !== "3D") return;

      // Toggle off if already orbiting.
      if (orbitTickerRef.current) {
        stopOrbit();

        return;
      }

      const ctrls = controlsRef.current;
      const cam = cameraRef.current;

      if (!ctrls || !cam) return;

      const target = ctrls.target.clone();
      const startRadius = cam.position.clone().sub(target);
      const radiusXZ = Math.sqrt(
        startRadius.x * startRadius.x + startRadius.z * startRadius.z,
      );
      const baseY = startRadius.y;
      // Base period at speed 1.0 is 4.5s/rev. Live-updates via the ref.
      const BASE_PERIOD_MS = 4500;
      let angle = 0;
      // The Y bob runs at half the azimuth frequency so XZ and Y waves drift
      // out of sync — gives the "shifting figure-8" look.
      let elev = 0;

      const ticker = (_time: number, deltaTime: number) => {
        const speed = Math.max(0.01, orbitSpeedRef.current);

        angle += (deltaTime / (BASE_PERIOD_MS / speed)) * Math.PI * 2;
        elev += (deltaTime / ((BASE_PERIOD_MS * 2) / speed)) * Math.PI * 2;
        const sin = Math.sin(angle);
        const cos = Math.cos(angle);
        const yOsc = Math.sin(elev) * radiusXZ * orbitYBobRef.current;

        cam.position.set(
          target.x + startRadius.x * cos + startRadius.z * sin,
          target.y + baseY + yOsc,
          target.z + -startRadius.x * sin + startRadius.z * cos,
        );
        ctrls.update();
      };

      orbitTickerRef.current = ticker;
      gsap.ticker.add(ticker);
    };

    window.addEventListener("keydown", onKey);

    return () => {
      window.removeEventListener("keydown", onKey);
      stopOrbit();
    };
  }, [viewMode]);

  // T key toggles tour playback. Works even while the UI is hidden (H), so a
  // tour can be started/stopped for a clean screen recording. Requires a
  // loaded tour config and 3D view (same gate as the Play tour button).
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.key !== "t" && e.key !== "T") return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const el = document.activeElement as HTMLElement | null;

      if (
        el &&
        (el.tagName === "INPUT" ||
          el.tagName === "TEXTAREA" ||
          el.isContentEditable)
      ) {
        return;
      }
      const tour = useTourStore.getState();

      if (tour.isPlaying) {
        tour.stop();

        return;
      }
      if (!tour.config || viewMode !== "3D") return;
      tour.start();
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, [viewMode]);

  // Effect 1b: Apply per-sample transforms to the live position buffer.
  // Reads base positions stashed during scene init and writes the result
  // into geometry.attributes.position. Marks needsUpdate so the GPU re-uploads.
  const centroidsCacheRef = useRef<{
    column: string | null;
    centroids: Map<string, [number, number]> | null;
  }>({ column: null, centroids: null });

  // Invalidate centroid cache when dataset changes.
  useEffect(() => {
    centroidsCacheRef.current = { column: null, centroids: null };
  }, [dataset]);

  useEffect(() => {
    const base = basePositionsRef.current;
    const pc = pointCloudRef.current;

    if (!base || !pc || !dataset) return;

    const positionAttr = pc.geometry.attributes
      .position as THREE.BufferAttribute;
    const out = positionAttr.array as Float32Array;

    if (!transformColumn || sampleTransforms.size === 0) {
      out.set(base);
      positionAttr.needsUpdate = true;
      invalidateRef.current();

      return;
    }

    const cluster =
      dataset.clusters?.find((c) => c.column === transformColumn) ?? null;

    if (!cluster) {
      // Column not loaded yet — fall back to base positions.
      out.set(base);
      positionAttr.needsUpdate = true;
      invalidateRef.current();

      return;
    }

    // Centroids depend only on (spatial, column); cache them across edits.
    let centroids = centroidsCacheRef.current.centroids;

    if (centroidsCacheRef.current.column !== transformColumn) {
      centroids = computeSampleCentroids(dataset.spatial, cluster);
      centroidsCacheRef.current = { column: transformColumn, centroids };
    }

    applyTransformsToPositions(
      base,
      out,
      cluster,
      centroids,
      sampleTransforms,
      sceneScaleRef.current,
    );
    positionAttr.needsUpdate = true;
    invalidateRef.current();
  }, [
    dataset,
    pointCloudVersion,
    transformColumn,
    transformVersion,
    clusterVersion,
  ]);

  // Effect 2: Update visualization based on mode array
  useEffect(() => {
    if (!pointCloudRef.current || !dataset) return;

    const updateVisualization = async () => {
      if (!pointCloudRef.current || !dataset) return;

      // Store cluster data for tooltip
      const selectedCluster = dataset.clusters?.find(
        (c) => c.column === selectedColumn,
      );

      if (selectedCluster) {
        clusterRef.current = selectedCluster;
        isNumericalClusterRef.current = selectedColumn
          ? getEffectiveColumnType(
              selectedColumn,
              dataset,
              columnTypeOverrides,
            ) === "numerical"
          : false;
        colorPaletteRef.current = selectedCluster.palette || colorPalette;
      }

      // Populate pinned columns for tooltip
      pinnedClustersRef.current = [];
      for (const col of pinnedTooltipColumns) {
        if (col === selectedColumn) continue; // active column already shown
        const cluster = dataset.clusters?.find((c) => c.column === col);

        if (!cluster) continue;
        const isNum =
          getEffectiveColumnType(col, dataset, columnTypeOverrides) ===
          "numerical";

        pinnedClustersRef.current.push({
          column: col,
          cluster,
          palette: cluster.palette || {},
          isNumerical: isNum,
        });
      }

      // Build advanced settings from store
      const adv: AdvancedVizSettings = {
        selectedSizeMultiplier,
        greyedOutSizeMultiplier,
        greyedOutAlpha,
        expressionAlphaMin,
        expressionAlphaMax,
        pointSizeMultiplierMin,
        pointSizeMultiplierMax,
      };

      // Determine which visualization to use based on mode array
      const hasGeneMode = mode.includes("gene");
      const hasCelltypeMode = mode.includes("celltype");

      // The 3D scene renders only the visible subset of the selection
      // (selectedCelltypes minus hiddenCelltypes, set via the per-badge eye
      // toggles). selectedCelltypes itself is untouched, so DEG / plots keep
      // the full selection.
      const effectiveCelltypes =
        hiddenCelltypes.size > 0
          ? new Set(
              [...selectedCelltypes].filter((ct) => !hiddenCelltypes.has(ct)),
            )
          : selectedCelltypes;

      // Check if the selected column is numerical (respects overrides)
      const isNumerical = selectedColumn
        ? getEffectiveColumnType(
            selectedColumn,
            dataset,
            columnTypeOverrides,
          ) === "numerical"
        : false;

      let result = null;

      // Check if gene has changed (for auto-scaling)
      const geneChanged = previousGeneRef.current !== selectedGene;

      if (geneChanged) {
        previousGeneRef.current = selectedGene;
      }

      // Sniff raw-counts vs log, apply the chosen display space, and report
      // whether the range should be re-auto-scaled (gene OR display-mode change).
      // Used by both the gene-only and combined paths so coloring, the scale
      // bar, and the hover tooltip share identical numbers.
      const toDisplayExpression = (expression: number[]) => {
        const kind = detectExpressionKind(expression);

        setDetectedGeneKind(kind);

        const effectiveMode = geneScaleModeUserSet
          ? geneScaleMode
          : nativeModeFor(kind);

        // Keep the toggle UI in sync with the auto-detected native default.
        if (!geneScaleModeUserSet && geneScaleMode !== effectiveMode) {
          setGeneScaleModeAuto(effectiveMode);
        }

        const modeChanged = previousModeRef.current !== effectiveMode;

        previousModeRef.current = effectiveMode;

        return {
          display: transformExpression(expression, kind, effectiveMode),
          autoScale: geneChanged || modeChanged,
        };
      };

      // Check if 2nd gene has changed (for coexpression auto-scaling).
      const gene2Changed = previousGene2Ref.current !== selectedGene2;

      if (gene2Changed) {
        previousGene2Ref.current = selectedGene2;
      }

      const coexpressActive =
        coexpressEnabled && !!selectedGene && !!selectedGene2;

      // Check if column or its effective type has changed (for auto-scaling numerical columns)
      const columnChanged = previousColumnRef.current !== selectedColumn;
      const typeChanged = previousColumnTypeRef.current !== isNumerical;

      if (columnChanged) {
        previousColumnRef.current = selectedColumn;
      }
      if (typeChanged) {
        previousColumnTypeRef.current = isNumerical;
      }

      const shouldAutoScale = columnChanged || typeChanged;

      if (coexpressActive) {
        // Two-gene coexpression: green + magenta additive blend.
        try {
          if (geneChanged || gene2Changed) {
            geneToastIdRef.current = toast.loading(
              `Loading "${selectedGene}" + "${selectedGene2}"…`,
            );
          }

          // Cache primary gene expression for tooltips (use gene 1).
          geneExpressionRef.current = await dataset.getGeneExpression(
            selectedGene as string,
          );

          const greyOutNonSelected =
            hasCelltypeMode && selectedCelltypes.size > 0 && !geneEverywhere;

          result = await updateCoexpressionVisualization(
            dataset,
            selectedGene as string,
            selectedGene2 as string,
            geneScaleMin,
            geneScaleMax,
            geneChanged ? setGeneScaleMin : undefined,
            geneChanged ? setGeneScaleMax : undefined,
            gene2ScaleMin,
            gene2ScaleMax,
            gene2Changed ? setGene2ScaleMin : undefined,
            gene2Changed ? setGene2ScaleMax : undefined,
            coexpressSwapped,
            alphaScale,
            sizeScale,
            selectedColumn,
            effectiveCelltypes,
            greyOutNonSelected,
            adv,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (
        hasGeneMode &&
        hasCelltypeMode &&
        selectedGene &&
        selectedCelltypes.size > 0 &&
        !geneEverywhere
      ) {
        // Combined mode: gene expression on selected celltypes
        try {
          if (geneChanged) {
            geneToastIdRef.current = toast.loading(
              `Loading expression for "${selectedGene}"...`,
            );
          }

          // Fetch gene expression, then apply the raw/log display transform.
          const raw = await dataset.getGeneExpression(selectedGene);
          const { display, autoScale } = raw
            ? toDisplayExpression(raw)
            : { display: raw, autoScale: geneChanged };

          geneExpressionRef.current = display;

          result = await updateCombinedVisualization(
            dataset,
            selectedGene,
            selectedColumn,
            effectiveCelltypes,
            alphaScale,
            sizeScale,
            geneScaleMin,
            geneScaleMax,
            autoScale ? setGeneScaleMin : undefined,
            autoScale ? setGeneScaleMax : undefined,
            adv,
            colormap,
            display,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (hasGeneMode && selectedGene) {
        // Gene mode only
        try {
          if (geneChanged) {
            geneToastIdRef.current = toast.loading(
              `Loading expression for "${selectedGene}"...`,
            );
          }

          // Fetch gene expression, then apply the raw/log display transform.
          const raw = await dataset.getGeneExpression(selectedGene);
          const { display, autoScale } = raw
            ? toDisplayExpression(raw)
            : { display: raw, autoScale: geneChanged };

          geneExpressionRef.current = display;

          result = await updateGeneVisualization(
            dataset,
            selectedGene,
            alphaScale,
            sizeScale,
            geneScaleMin,
            geneScaleMax,
            autoScale ? setGeneScaleMin : undefined,
            autoScale ? setGeneScaleMax : undefined,
            adv,
            colormap,
            display,
          );
        } finally {
          if (geneToastIdRef.current != null) {
            toast.dismiss(geneToastIdRef.current);
            geneToastIdRef.current = null;
          }
        }
      } else if (hasCelltypeMode) {
        // Celltype mode only

        // Use appropriate visualization function based on column type
        result = isNumerical
          ? updateNumericalCelltypeVisualization(
              dataset,
              selectedColumn,
              alphaScale,
              sizeScale,
              numericalScaleMin,
              numericalScaleMax,
              shouldAutoScale ? setNumericalScaleMin : undefined,
              shouldAutoScale ? setNumericalScaleMax : undefined,
              adv,
              colormap,
            )
          : updateCelltypeVisualization(
              dataset,
              selectedColumn,
              effectiveCelltypes,
              colorPalette,
              alphaScale,
              sizeScale,
              adv,
            );
      }

      // Update point cloud if we have a result
      if (result && pointCloudRef.current) {
        updatePointCloudAttributes(
          pointCloudRef.current,
          result.colors,
          result.sizes,
          result.alphas,
        );
        // Apply current sizeScale to dotSize uniform
        updateDotSize(
          pointCloudRef.current,
          baseDotSizeRef.current * sizeScale,
        );
        invalidateRef.current();
      }

      // E2E timing hook: mark this render complete (covers initial render and
      // every gene/celltype change). Stats are read straight off the dataset.
      markRenderComplete({
        pointCount: dataset.getPointCount(),
        geneCount: dataset.genes.length,
        // Genes with expression actually loaded — differs from geneCount when
        // an adapter listed gene names but skipped the matrix (e.g. the
        // Xenium h5 dense-size guard).
        exprGeneCount:
          dataset.matrix instanceof Map ? dataset.matrix.size : null,
        dimensions: dataset.spatial.dimensions,
        dataType: "single_cell",
      });
      if (mode.includes("gene") && selectedGene) {
        markGeneRendered(selectedGene, dataset.getPointCount());
      }
    };

    updateVisualization();
  }, [
    dataset,
    selectedGene,
    selectedColumn,
    selectedCelltypes,
    geneEverywhere,
    hiddenCelltypes,
    colorPalette,
    alphaScale,
    geneScaleMin,
    geneScaleMax,
    geneScaleMode,
    numericalScaleMin,
    numericalScaleMax,
    mode,
    clusterVersion,
    columnTypeOverrides,
    pointCloudVersion,
    selectedSizeMultiplier,
    greyedOutSizeMultiplier,
    greyedOutAlpha,
    expressionAlphaMin,
    expressionAlphaMax,
    pointSizeMultiplierMin,
    pointSizeMultiplierMax,
    pinnedTooltipColumns,
    colormap,
    coexpressEnabled,
    selectedGene2,
    coexpressSwapped,
    gene2ScaleMin,
    gene2ScaleMax,
  ]);

  // Effect 3: Update dotSize uniform when slider or targetPx changes (instant)
  useEffect(() => {
    if (pointCloudRef.current) {
      const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2);
      const newBaseDotSize =
        (targetPx * cameraDistanceRef.current) /
        (VISUALIZATION_CONFIG.POINT_BASE_SIZE * proj11);

      baseDotSizeRef.current = newBaseDotSize;
      updateDotSize(pointCloudRef.current, newBaseDotSize * sizeScale);
      invalidateRef.current();
    }
  }, [sizeScale, targetPx]);

  // Effect 4: Apply scene transforms (rotation, flip) — instant via group matrix
  useEffect(() => {
    const group = sceneGroupRef.current;

    if (!group) return;

    group.rotation.z = (sceneRotation * Math.PI) / 180;
    group.scale.x = flipX ? -1 : 1;
    group.scale.y = flipY ? -1 : 1;
    invalidateRef.current();
  }, [sceneRotation, flipX, flipY]);

  // Distance-filter uniforms: pushed every frame because controls.target
  // mutates outside React (cmd+click pan, mouse drag, R reset, orbit). Cost
  // is a handful of Vector3 copies per frame — negligible.
  useEffect(() => {
    const last = { x: NaN, y: NaN, z: NaN, radius: NaN, enabled: NaN };
    const pushFilterUniforms = () => {
      const ctrls = controlsRef.current;
      const pc = pointCloudRef.current;

      if (!ctrls) return;
      const enabled = targetFilterEnabledRef.current ? 1.0 : 0.0;
      const radius = targetFilterRadiusRef.current;
      const feather = Math.max(1.0, radius * 0.1);
      const center = ctrls.target;

      if (
        center.x === last.x &&
        center.y === last.y &&
        center.z === last.z &&
        radius === last.radius &&
        enabled === last.enabled
      ) {
        return;
      }
      last.x = center.x;
      last.y = center.y;
      last.z = center.z;
      last.radius = radius;
      last.enabled = enabled;

      if (pc) {
        const mat = pc.material as THREE.ShaderMaterial;

        if (mat.uniforms.uTargetCenter) {
          (mat.uniforms.uTargetCenter.value as THREE.Vector3).copy(center);
          mat.uniforms.uTargetRadius.value = radius;
          mat.uniforms.uTargetFeather.value = feather;
          mat.uniforms.uTargetFilterEnabled.value = enabled;
        }
      }

      smLayerRef.current?.setTargetFilter(
        center,
        radius,
        feather,
        targetFilterEnabledRef.current,
      );
      invalidateRef.current();
    };

    gsap.ticker.add(pushFilterUniforms);

    return () => gsap.ticker.remove(pushFilterUniforms);
  }, []);

  return (
    <>
      <div
        ref={containerRef}
        className="absolute inset-0 w-full h-full"
        data-testid="sc-scene-canvas"
        style={{ margin: 0, padding: 0 }}
      />

      {webglFailed && <WebGLErrorOverlay />}

      <LayerControls slot="panel" />
      <LayerBindings ctx={layerCtx} />

      {/* SM loading indicator */}
      {smLoading && (
        <div className="absolute top-28 right-6 z-50 flex items-center gap-2 px-3 py-2 rounded-lg bg-black/60 backdrop-blur-sm text-white text-sm">
          <div className="w-4 h-4 border-2 border-white/30 border-t-white rounded-full animate-spin" />
          Loading molecules...
        </div>
      )}

      {/* Last cmd/ctrl+click coordinates. */}
      {lastClickCoords && (
        <div className="absolute bottom-4 left-4 z-50 px-3 py-1.5 rounded-lg bg-black/60 backdrop-blur-sm text-white text-xs font-mono">
          ({lastClickCoords.x.toFixed(2)}, {lastClickCoords.y.toFixed(2)},{" "}
          {lastClickCoords.z.toFixed(2)})
        </div>
      )}

      {/* SC + SM legends share a single absolute container so the SM gene
          pills stack below the SC celltype/gene legend instead of overlapping
          it at top-right. */}
      <div className="absolute right-6 top-24 z-10 flex flex-col items-end gap-4 max-w-xs">
        <VisualizationLegends embedded />
        <SingleMoleculeLegends embedded />
      </div>

      {/* Spatial scale bar — only for datasets with reliable raw coordinates */}
      {dataset && !dataset.metadata?.wasNormalized && (
        <SpatialScaleBar
          cameraRef={
            cameraRef as React.RefObject<THREE.PerspectiveCamera | null>
          }
          controlsRef={controlsRef}
          rendererRef={rendererRef}
        />
      )}

      {/* Export box overlay — gated by viewMode + raw-coords (same as scale bar). */}
      <ExportBoxOverlay
        cameraRef={cameraRef as React.RefObject<THREE.PerspectiveCamera | null>}
        canShow={
          viewMode === "2D" && !!dataset && !dataset.metadata?.wasNormalized
        }
        controlsRef={controlsRef}
        rendererRef={rendererRef}
        setCenterPx={setExportBoxCenterPx}
        state={{
          enabled: exportBoxEnabled,
          widthMm: exportBoxWidthMm,
          heightMm: exportBoxHeightMm,
          centerPx: exportBoxCenterPx,
        }}
      />
    </>
  );
}
