"use client";

import type { StandardizedDataset } from "@/lib/StandardizedDataset";
import type { LmMenu } from "@/lib/stores/createLabelledMoleculeVisualizationStore";

import { useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";

import { resolveValueColor } from "@/lib/stores/createLabelledMoleculeVisualizationStore";
import {
  labelledMoleculeVisualizationStore,
  useLabelledMoleculeVisualizationStore,
} from "@/lib/stores/labelledMoleculeVisualizationStore";
import {
  buildPaletteLut,
  buildSelectionLut,
} from "@/lib/webgl/labelled-molecule-lut";
import {
  labelledMoleculeFragmentShader,
  labelledMoleculeVertexShader,
} from "@/lib/webgl/labelled-molecule-shaders";
import { glassPanel } from "@/components/primitives";
import { initializeScene } from "@/lib/webgl/scene-manager";

interface ColumnData {
  column: string;
  uniqueValues: string[];
  valueIndices: ArrayLike<number>;
  palette: Record<string, string> | null;
}

interface Props {
  dataset: StandardizedDataset;
  /** Bumped by the page when a lazily-loaded column arrives. */
  clusterVersion?: number;
}

/** A vertex attribute wide enough for the column's category count. */
function indexAttribute(
  indices: ArrayLike<number>,
  count: number,
  uniqueCount: number,
) {
  const Arr =
    uniqueCount <= 0xff
      ? Uint8Array
      : uniqueCount <= 0xffff
        ? Uint16Array
        : Float32Array;
  const arr = new Arr(count);

  for (let i = 0; i < count; i++) arr[i] = indices[i];

  // normalized = false, so an integer attribute reaches the shader as its
  // exact float value rather than being scaled into [0,1].
  return new THREE.BufferAttribute(arr, 1, false);
}

function makeLutTexture(data: Uint8Array, channels: 1 | 4): THREE.DataTexture {
  const width = data.length / channels;
  const tex = new THREE.DataTexture(
    data,
    width,
    1,
    channels === 1 ? THREE.RedFormat : THREE.RGBAFormat,
    THREE.UnsignedByteType,
  );

  // Nearest sampling: these are lookup tables, not images — any filtering
  // would blend neighbouring categories together.
  tex.magFilter = THREE.NearestFilter;
  tex.minFilter = THREE.NearestFilter;
  tex.generateMipmaps = false;
  // Single-channel rows are byte-packed at arbitrary widths (213 genes), so
  // the default 4-byte unpack alignment would misread them.
  tex.unpackAlignment = channels === 1 ? 1 : 4;
  tex.needsUpdate = true;

  return tex;
}

export default function LabelledMoleculeThreeScene({
  dataset,
  clusterVersion = 0,
}: Props) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const pointsRef = useRef<THREE.Points | null>(null);
  const materialRef = useRef<THREE.ShaderMaterial | null>(null);
  const texturesRef = useRef<THREE.DataTexture[]>([]);
  const cameraRef = useRef<THREE.PerspectiveCamera | null>(null);
  const rendererRef = useRef<THREE.WebGLRenderer | null>(null);
  const controlsTargetRef = useRef<THREE.Vector3 | null>(null);
  const raycasterRef = useRef(new THREE.Raycaster());
  const mouseRef = useRef(new THREE.Vector2());
  const [hover, setHover] = useState<{
    x: number;
    y: number;
    rows: { menu: LmMenu; value: string; color: string }[];
  } | null>(null);
  const [glError, setGlError] = useState<string | null>(null);
  // Bumped whenever the scene effect builds a NEW material. The LUT and scalar
  // effects below key off it: without it, rebuilding the scene (a 2D/3D toggle)
  // would leave the fresh material with null lookup textures, because none of
  // their own dependencies changed.
  const [materialVersion, setMaterialVersion] = useState(0);

  const {
    colorBy,
    domainVariant,
    selections,
    colorOverrides,
    sizeOverrides,
    geneColorSlots,
    unselectedMode,
    unselectedAlpha,
    globalScale,
    globalAlpha,
    viewMode,
  } = useLabelledMoleculeVisualizationStore();

  /** The column backing each menu, resolved through the domain variant. */
  const columnFor = useMemo(
    (): Record<LmMenu, string> => ({
      gene: "gene",
      domain: domainVariant,
      cell: "cell",
    }),
    [domainVariant],
  );

  /** Menu -> its loaded column, or null while it is still being fetched. */
  const columns = useMemo(() => {
    const find = (name: string): ColumnData | null => {
      const c = dataset.clusters?.find((cl) => cl.column === name);

      if (!c?.uniqueValues || !c.valueIndices) return null;

      return {
        column: c.column,
        uniqueValues: c.uniqueValues,
        valueIndices: c.valueIndices,
        palette: c.palette,
      };
    };

    return {
      gene: find(columnFor.gene),
      domain: find(columnFor.domain),
      cell: find(columnFor.cell),
    } as Record<LmMenu, ColumnData | null>;
    // clusterVersion is the signal that a lazily-loaded column has landed.
  }, [dataset, columnFor, clusterVersion]);

  const ready = !!(columns.gene && columns.domain && columns.cell);

  // ── Scene + geometry. Rebuilt only when the data or the camera mode
  //    changes; selection changes never touch this.
  useEffect(() => {
    const container = containerRef.current;

    if (!container || !ready) return;

    let setup;

    try {
      setup = initializeScene(container, { is2D: viewMode === "2D" });
    } catch (e) {
      setGlError(e instanceof Error ? e.message : String(e));

      return;
    }
    setGlError(null);

    const { scene, camera, controls, animate, dispose } = setup;
    const count = dataset.getPointCount();
    const coords = dataset.spatial.coordinates as Float32Array;
    const dims = dataset.spatial.dimensions;

    // Centre the cloud on the origin so the trackball orbits its middle.
    const positions = new Float32Array(count * 3);
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];

    for (let i = 0; i < count; i++) {
      for (let d = 0; d < 3; d++) {
        const v = d < dims ? coords[i * dims + d] : 0;

        positions[i * 3 + d] = v;
        if (v < min[d]) min[d] = v;
        if (v > max[d]) max[d] = v;
      }
    }
    const center = min.map((lo, d) => (lo + max[d]) / 2);
    const extent = Math.max(...max.map((hi, d) => hi - min[d]), 1);

    for (let i = 0; i < count; i++) {
      for (let d = 0; d < 3; d++) positions[i * 3 + d] -= center[d];
    }

    const geometry = new THREE.BufferGeometry();

    geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute(
      "aGene",
      indexAttribute(
        columns.gene!.valueIndices,
        count,
        columns.gene!.uniqueValues.length,
      ),
    );
    geometry.setAttribute(
      "aDomain",
      indexAttribute(
        columns.domain!.valueIndices,
        count,
        columns.domain!.uniqueValues.length,
      ),
    );
    geometry.setAttribute(
      "aCell",
      indexAttribute(
        columns.cell!.valueIndices,
        count,
        columns.cell!.uniqueValues.length,
      ),
    );

    // The shader computes gl_PointSize = aSize * sizeMul * uGlobalSize
    // * dotSize * proj[1][1] / -mvPosition.z, so dotSize has to be
    // back-calculated from the data extent to land on a few pixels at the
    // opening zoom — same derivation as three-scene.tsx. Getting this wrong by
    // orders of magnitude puts every point under the sub-pixel cull below.
    const cameraDistance = extent * 1.6;
    const proj11 = 1.0 / Math.tan((75 * Math.PI) / 180 / 2); // ~1.303
    const TARGET_PX = 3;
    const baseDotSize = (TARGET_PX * cameraDistance) / proj11;

    const material = new THREE.ShaderMaterial({
      vertexShader: labelledMoleculeVertexShader,
      fragmentShader: labelledMoleculeFragmentShader,
      transparent: true,
      depthWrite: false,
      uniforms: {
        uSelGene: { value: null },
        uSelDomain: { value: null },
        uSelCell: { value: null },
        uPalette: { value: null },
        uNGene: { value: columns.gene!.uniqueValues.length },
        uNDomain: { value: columns.domain!.uniqueValues.length },
        uNCell: { value: columns.cell!.uniqueValues.length },
        uNColorBy: { value: 1 },
        uColorBy: { value: 2 },
        dotSize: { value: baseDotSize },
        uGlobalSize: { value: 1 },
        uGlobalAlpha: { value: 1 },
        uUnselectedAlpha: { value: 0.2 },
        uUnselectedHidden: { value: 0 },
        uUnselectedSize: { value: 0.6 },
        uUnselectedColor: { value: new THREE.Color(0x555555) },
        uShape: { value: 0 },
      },
    });

    // aSize / aAlpha are reserved for genuinely per-molecule continuous values
    // (correlation, brightness, scoreA). Until one is in use they are constant,
    // so we supply them as default attribute values rather than paying for two
    // Float32Arrays of length N.
    material.defaultAttributeValues = {
      ...(material.defaultAttributeValues ?? {}),
      aSize: [1],
      aAlpha: [1],
    };

    const points = new THREE.Points(geometry, material);

    scene.add(points);
    pointsRef.current = points;
    materialRef.current = material;
    setMaterialVersion((v) => v + 1);

    camera.position.set(0, 0, cameraDistance);
    camera.lookAt(0, 0, 0);
    controls.target.set(0, 0, 0);
    controls.update();

    cameraRef.current = camera;
    rendererRef.current = setup.renderer;
    controlsTargetRef.current = controls.target;

    // ── Hover picking.
    // Only molecules that pass the filter are pickable: with grey-mode on there
    // are millions of dimmed points, and letting them win the depth test would
    // make the selection nearly impossible to inspect.
    const passesFilter = (i: number) => {
      const { selections: sel } = labelledMoleculeVisualizationStore.getState();

      for (const menu of ["gene", "domain", "cell"] as LmMenu[]) {
        const s = sel[menu];

        if (s.size === 0) continue;
        const col = columns[menu]!;

        if (!s.has(col.uniqueValues[col.valueIndices[i]])) return false;
      }

      return true;
    };

    const describe = (i: number) => {
      const st = labelledMoleculeVisualizationStore.getState();

      return (["gene", "domain", "cell"] as LmMenu[]).map((menu) => {
        const col = columns[menu]!;
        const value = col.uniqueValues[col.valueIndices[i]];

        return {
          menu,
          value,
          color: resolveValueColor(
            menu,
            value,
            st.colorOverrides[menu],
            st.geneColorSlots,
            col.palette,
          ),
        };
      });
    };

    let lastCheck = 0;

    const handleMouseMove = (event: MouseEvent) => {
      const el = setup.renderer.domElement;
      const rect = el.getBoundingClientRect();

      mouseRef.current.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
      mouseRef.current.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

      const now = Date.now();

      // Raycasting a 3M-point cloud is O(N); throttle so dragging stays smooth.
      if (now - lastCheck < 50) return;
      lastCheck = now;

      // Convert a fixed pixel radius to world units so the pick target stays
      // the same on screen at any zoom (same derivation as three-scene.tsx).
      const target = controlsTargetRef.current;
      const camDist = target
        ? camera.position.distanceTo(target)
        : camera.position.length();
      const fovRad = (camera.fov / 2) * (Math.PI / 180);
      const pixelSize = (2 * camDist * Math.tan(fovRad)) / rect.height;

      raycasterRef.current.params.Points!.threshold = pixelSize * 5;
      raycasterRef.current.setFromCamera(mouseRef.current, camera);

      const hits = raycasterRef.current.intersectObject(points);

      hits.sort((a, b) => a.distance - b.distance);

      const hit = hits.find((h) => h.index != null && passesFilter(h.index));

      setHover(
        hit
          ? {
              x: event.clientX - rect.left,
              y: event.clientY - rect.top,
              rows: describe(hit.index!),
            }
          : null,
      );
    };

    const handleMouseLeave = () => setHover(null);

    setup.renderer.domElement.addEventListener("mousemove", handleMouseMove);
    setup.renderer.domElement.addEventListener("mouseleave", handleMouseLeave);

    animate();

    return () => {
      setup.renderer.domElement.removeEventListener(
        "mousemove",
        handleMouseMove,
      );
      setup.renderer.domElement.removeEventListener(
        "mouseleave",
        handleMouseLeave,
      );
      dispose();
      geometry.dispose();
      material.dispose();
      for (const t of texturesRef.current) t.dispose();
      texturesRef.current = [];
      pointsRef.current = null;
      materialRef.current = null;
      cameraRef.current = null;
      rendererRef.current = null;
      setHover(null);
    };
  }, [dataset, ready, viewMode, columnFor, clusterVersion]);

  // ── Selection LUTs. This is the hot path: a checkbox toggle re-uploads
  //    three textures of one byte per category and nothing else.
  useEffect(() => {
    const material = materialRef.current;

    if (!material || !ready) return;

    const menus: LmMenu[] = ["gene", "domain", "cell"];
    const uniforms = ["uSelGene", "uSelDomain", "uSelCell"] as const;

    menus.forEach((menu, i) => {
      const col = columns[menu]!;
      const lut = buildSelectionLut(col.uniqueValues, selections[menu]);
      const tex = makeLutTexture(lut, 1);

      material.uniforms[uniforms[i]].value?.dispose?.();
      material.uniforms[uniforms[i]].value = tex;
    });
  }, [columns, selections, ready, materialVersion]);

  // ── Palette LUT for the colouring column.
  useEffect(() => {
    const material = materialRef.current;

    if (!material || !ready) return;

    const col = columns[colorBy]!;
    const lut = buildPaletteLut(
      col.uniqueValues,
      (v) =>
        resolveValueColor(
          colorBy,
          v,
          colorOverrides[colorBy],
          geneColorSlots,
          col.palette,
        ),
      (v) => sizeOverrides[colorBy][v] ?? 1,
    );

    material.uniforms.uPalette.value?.dispose?.();
    material.uniforms.uPalette.value = makeLutTexture(lut, 4);
    material.uniforms.uNColorBy.value = col.uniqueValues.length;
    material.uniforms.uColorBy.value =
      colorBy === "gene" ? 0 : colorBy === "domain" ? 1 : 2;
  }, [
    columns,
    colorBy,
    colorOverrides,
    sizeOverrides,
    geneColorSlots,
    ready,
    materialVersion,
  ]);

  // ── Scalars.
  useEffect(() => {
    const material = materialRef.current;

    if (!material) return;
    material.uniforms.uGlobalSize.value = globalScale;
    material.uniforms.uGlobalAlpha.value = globalAlpha;
    material.uniforms.uUnselectedAlpha.value = unselectedAlpha;
    material.uniforms.uUnselectedHidden.value =
      unselectedMode === "hidden" ? 1 : 0;
  }, [
    globalScale,
    globalAlpha,
    unselectedAlpha,
    unselectedMode,
    materialVersion,
  ]);

  if (glError) {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black p-8 text-center">
        <div className="max-w-md">
          <p className="mb-2 text-lg text-white">3D view unavailable</p>
          <p className="text-sm text-default-400">
            Your browser couldn&apos;t create a WebGL context, so the scene
            can&apos;t render. The dataset itself loaded fine.
          </p>
          <p className="mt-3 font-mono text-xs text-default-500">{glError}</p>
        </div>
      </div>
    );
  }

  const MENU_LABEL: Record<LmMenu, string> = {
    gene: "Gene",
    domain: "Domain",
    cell: "Cell",
  };

  return (
    <div className="absolute inset-0 bg-black">
      <div ref={containerRef} className="h-full w-full" />

      {hover && (
        <div
          className={`pointer-events-none absolute z-[var(--z-legends)] rounded-xl px-3 py-2 text-xs ${glassPanel()}`}
          data-testid="lm-hover-tooltip"
          style={{ left: hover.x + 12, top: hover.y + 12 }}
        >
          {hover.rows.map((r) => (
            <div key={r.menu} className="flex items-center gap-2 py-0.5">
              <span
                className="h-2.5 w-2.5 shrink-0 rounded-full"
                style={{ backgroundColor: r.color }}
              />
              <span className="text-default-500">{MENU_LABEL[r.menu]}</span>
              <span className="max-w-[220px] truncate font-medium">
                {r.menu === "gene" ? r.value.split(" (")[0] : r.value}
              </span>
            </div>
          ))}
        </div>
      )}
      {!ready && (
        <div className="absolute inset-0 flex items-center justify-center">
          <p className="text-sm text-default-400">Loading columns…</p>
        </div>
      )}
    </div>
  );
}
