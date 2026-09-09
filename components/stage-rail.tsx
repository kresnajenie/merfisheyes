"use client";

import type { DatasetPreview } from "@/lib/webgl/preview";

import { Button } from "@heroui/button";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";

import { glassPanel } from "@/components/primitives";
import { loadPreview } from "@/lib/webgl/preview";

export interface ProjectDatasetSummary {
  id: string;
  title: string;
  datasetType: string | null;
  numCells: number;
  s3BaseUrl: string;
  stage: string | null;
  sortOrder: number;
}

interface Stage {
  stage: string;
  members: ProjectDatasetSummary[];
  /** The one whose preview the tile renders: first in project order. */
  rep: ProjectDatasetSummary;
}

interface Props {
  projectId: string;
  /** s3BaseUrl of the dataset currently open, to mark the active stage. */
  currentUrl?: string | null;
  onOpen: (d: ProjectDatasetSummary) => void;
  onSplit: (d: ProjectDatasetSummary) => void;
}

const TILE = 96;
/** The stage caption under each tile: 12px line plus its 2px top margin. */
const LABEL_H = 14;
/** Vertical padding of the bar, top + bottom. */
const PAD_Y = 12;

/**
 * Idle height of the rail; the viewer lifts its bottom-left chrome by this.
 *
 * Deliberately the *idle* height, not the magnified one: the bar grows upward
 * on hover and the chrome above it must not jump every time the pointer
 * crosses the strip.
 */
export const STAGE_RAIL_HEIGHT = TILE + LABEL_H + PAD_Y;

/** Dock magnification: the tile under the pointer, and how far it spreads. */
const MAG_MAX = 1.4;
const MAG_NEIGHBOURS = 2;
const MAG_RAIL_HEIGHT = Math.round(TILE * MAG_MAX) + LABEL_H + PAD_Y;

/**
 * Mac-dock falloff. Squared rather than linear so the tile under the pointer
 * stands clearly proud of its neighbours instead of the whole run swelling.
 */
function tileScale(index: number, focus: number) {
  const d = Math.abs(index - focus);

  if (d > MAG_NEIGHBOURS) return 1;
  const t = 1 - d / (MAG_NEIGHBOURS + 1);

  return 1 + (MAG_MAX - 1) * t * t;
}

/**
 * Allen-atlas style rail: one tile per developmental stage, each a live
 * preview rather than a screenshot.
 *
 * All tiles share ONE WebGL context, drawn with per-tile scissor rectangles.
 * A context each would exhaust the browser's limit (~16) and cost far more
 * memory; one canvas over the strip, with viewports recomputed from the tiles'
 * DOM rects each frame, also means scrolling needs no special handling.
 */
export default function StageRail({
  projectId,
  currentUrl,
  onOpen,
  onSplit,
}: Props) {
  const [stages, setStages] = useState<Stage[] | null>(null);
  const [hovered, setHovered] = useState<string | null>(null);
  // A pinned stage keeps its card open after the pointer leaves, so the embryo
  // list can be read and clicked without chasing a hover.
  const [pinned, setPinned] = useState<string | null>(null);
  // Which embryo within the open stage the card previews. Falls back to the
  // stage's representative whenever it names one from a different stage, so
  // changing stage needs no reset.
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const shown = pinned ?? hovered;
  const open = shown ? (stages?.find((s) => s.stage === shown) ?? null) : null;
  const selected =
    open?.members.find((m) => m.id === selectedId) ?? open?.rep ?? null;
  const selectedUrl = selected?.s3BaseUrl ?? null;
  // Read inside the render loop; as state it would re-run the effect below on
  // every hover, tearing down and re-fetching all previews.
  const hoveredRef = useRef<string | null>(null);
  /** Hovered or pinned — whichever tile is being looked at, so it stops spinning. */
  const holdRef = useRef<string | null>(null);

  // Written during render: idempotent, and it keeps the render loop's notion of
  // "held" in step with both the pointer and the pin.
  holdRef.current = hovered ?? pinned;
  const [previewsReady, setPreviewsReady] = useState(0);
  const [loaded, setLoaded] = useState<Set<string>>(new Set());
  /** Keyed by s3BaseUrl — a stage's tile shows its rep, the card any member. */
  const previewsRef = useRef(new Map<string, DatasetPreview>());
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const stripRef = useRef<HTMLDivElement | null>(null);
  const tileRefs = useRef(new Map<string, HTMLDivElement>());
  const scenesRef = useRef(
    new Map<
      string,
      {
        scene: THREE.Scene;
        camera: THREE.PerspectiveCamera;
        points: THREE.Points;
      }
    >(),
  );
  const bigCanvasRef = useRef<HTMLCanvasElement | null>(null);

  // ── Project members, grouped into stages in project order.
  useEffect(() => {
    let alive = true;

    fetch(`/api/projects/${projectId}/public`)
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => {
        if (!alive || !j) return;
        const byStage = new Map<string, ProjectDatasetSummary[]>();

        for (const d of j.datasets as ProjectDatasetSummary[]) {
          if (!d.stage) continue;
          if (!byStage.has(d.stage)) byStage.set(d.stage, []);
          byStage.get(d.stage)!.push(d);
        }
        setStages(
          [...byStage.entries()].map(([stage, members]) => ({
            stage,
            members,
            rep: members[0],
          })),
        );
      })
      .catch(() => {});

    return () => {
      alive = false;
    };
  }, [projectId]);

  // ── Preview data. Independent of the renderer, so a browser without WebGL
  //    still fetches nothing it cannot use and one without previews degrades to
  //    plain labelled buttons.
  useEffect(() => {
    if (!stages?.length) return;
    const controllers: AbortController[] = [];

    for (const s of stages) {
      const ac = new AbortController();

      controllers.push(ac);
      loadPreview(s.rep.s3BaseUrl, ac.signal).then((p) => {
        if (!p) return;
        previewsRef.current.set(s.rep.s3BaseUrl, p);
        setLoaded((prev) => new Set(prev).add(s.rep.s3BaseUrl));
        setPreviewsReady((n) => n + 1);
      });
    }

    return () => {
      for (const c of controllers) c.abort();
      previewsRef.current.clear();
      setLoaded(new Set());
    };
  }, [stages]);

  // ── One renderer for every tile.
  useEffect(() => {
    const canvas = canvasRef.current;

    if (!canvas || !stages?.length) return;

    let renderer: THREE.WebGLRenderer;

    try {
      renderer = new THREE.WebGLRenderer({
        canvas,
        antialias: true,
        alpha: true,
      });
    } catch {
      return; // no WebGL: the tiles still work as labelled buttons
    }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setScissorTest(true);

    const scenes = scenesRef.current;

    for (const [url, p] of previewsRef.current) {
      if (scenes.has(url)) continue;
      const geom = new THREE.BufferGeometry();

      geom.setAttribute("position", new THREE.BufferAttribute(p.positions, 3));
      geom.setAttribute("color", new THREE.BufferAttribute(p.colors, 3));

      const points = new THREE.Points(
        geom,
        new THREE.PointsMaterial({ size: 0.012, vertexColors: true }),
      );
      const scene = new THREE.Scene();

      scene.add(points);

      const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);

      camera.position.set(0, 0, 1.6);
      camera.lookAt(0, 0, 0);
      scenes.set(url, { scene, camera, points });
    }

    let raf = 0;
    let last = performance.now();
    const sized = { w: 0, h: 0 };

    const frame = () => {
      raf = requestAnimationFrame(frame);
      const strip = stripRef.current;

      if (!strip) return;
      const now = performance.now();
      const dt = (now - last) / 1000;

      last = now;

      const stripRect = strip.getBoundingClientRect();

      // Compare against the last CSS size, not canvas.width: that is scaled by
      // the pixel ratio, so on a HiDPI display it never matches and the drawing
      // buffer was being reallocated every frame.
      if (sized.w !== stripRect.width || sized.h !== stripRect.height) {
        sized.w = stripRect.width;
        sized.h = stripRect.height;
        renderer.setSize(stripRect.width, stripRect.height, false);
      }

      renderer.setClearColor(0x000000, 0);
      renderer.clear();

      for (const s of stages) {
        const entry = scenes.get(s.rep.s3BaseUrl);
        const el = tileRefs.current.get(s.stage);

        if (!entry || !el) continue;
        const r = el.getBoundingClientRect();

        // Cull tiles scrolled out of the strip.
        if (r.right < stripRect.left || r.left > stripRect.right) continue;

        // Slow spin so the preview reads as 3D; paused while hovered so a
        // drag-free look is possible.
        if (holdRef.current !== s.stage) entry.points.rotation.y += dt * 0.35;

        const x = r.left - stripRect.left;
        const y = stripRect.bottom - r.bottom; // WebGL origin is bottom-left

        renderer.setViewport(x, y, r.width, r.height);
        renderer.setScissor(x, y, r.width, r.height);
        entry.camera.aspect = r.width / r.height;
        entry.camera.updateProjectionMatrix();
        renderer.render(entry.scene, entry.camera);
      }
    };

    frame();

    return () => {
      cancelAnimationFrame(raf);
      for (const { scene, points } of scenes.values()) {
        points.geometry.dispose();
        (points.material as THREE.Material).dispose();
        scene.clear();
      }
      scenes.clear();
      renderer.dispose();
    };
  }, [stages, previewsReady]);

  // ── A non-representative member has no preview yet; fetch on selection.
  useEffect(() => {
    if (!selectedUrl || previewsRef.current.has(selectedUrl)) return;
    const ac = new AbortController();

    loadPreview(selectedUrl, ac.signal).then((p) => {
      if (!p) return;
      previewsRef.current.set(selectedUrl, p);
      setLoaded((prev) => new Set(prev).add(selectedUrl));
      setPreviewsReady((n) => n + 1);
    });

    return () => ac.abort();
  }, [selectedUrl]);

  // ── Enlarged preview inside the hover card.
  //
  // Its own context: the strip's canvas lives inside the strip and cannot reach
  // the card above it. Two contexts is still far under the browser's limit and
  // much simpler than reshaping the shared canvas to span both.
  useEffect(() => {
    const canvas = bigCanvasRef.current;
    const entry = selectedUrl ? scenesRef.current.get(selectedUrl) : null;

    if (!canvas || !entry) return;

    let renderer: THREE.WebGLRenderer;

    try {
      renderer = new THREE.WebGLRenderer({
        canvas,
        antialias: true,
        alpha: true,
      });
    } catch {
      return;
    }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    renderer.setSize(canvas.clientWidth, canvas.clientHeight, false);

    // Its own camera: the tile's is framed for a 96px square.
    const camera = new THREE.PerspectiveCamera(
      45,
      canvas.clientWidth / canvas.clientHeight,
      0.01,
      100,
    );

    camera.position.set(0, 0, 1.5);
    camera.lookAt(0, 0, 0);

    let raf = 0;
    const frame = () => {
      raf = requestAnimationFrame(frame);
      renderer.setClearColor(0x000000, 0);
      renderer.render(entry.scene, camera);
    };

    frame();

    return () => {
      cancelAnimationFrame(raf);
      renderer.dispose();
    };
  }, [selectedUrl, previewsReady]);

  const hover = useCallback((stage: string | null) => {
    hoveredRef.current = stage;
    setHovered(stage);
  }, []);

  const setTileRef = useCallback((stage: string, el: HTMLDivElement | null) => {
    if (el) tileRefs.current.set(stage, el);
    else tileRefs.current.delete(stage);
  }, []);

  const activeStage = useMemo(() => {
    if (!currentUrl || !stages) return null;
    const norm = currentUrl.replace(/\/+$/, "");

    return (
      stages.find((s) =>
        s.members.some((m) => m.s3BaseUrl.replace(/\/+$/, "") === norm),
      )?.stage ?? null
    );
  }, [currentUrl, stages]);

  // Rendered before the project resolves so the rail doesn't pop in once the
  // scene is already up. Sized identically to the real strip.
  if (!stages) {
    return (
      <div
        data-ui-overlay
        className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
      >
        <div
          className={`flex items-center justify-center gap-2 overflow-hidden px-2 py-1.5 ${glassPanel()} rounded-none`}
          style={{ height: STAGE_RAIL_HEIGHT }}
        >
          {Array.from({ length: 8 }).map((_, i) => (
            <div key={i} className="shrink-0">
              <div
                className="animate-pulse rounded-lg bg-default-200/40"
                style={{ width: TILE, height: TILE }}
              />
              <div className="mx-auto mt-1 h-2 w-12 animate-pulse rounded bg-default-200/30" />
            </div>
          ))}
        </div>
      </div>
    );
  }

  if (!stages.length) return null;

  const focus = hovered ? stages.findIndex((s) => s.stage === hovered) : -1;

  return (
    <div
      data-ui-overlay
      className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
      onMouseLeave={() => hover(null)}
    >
      {/* Stage card: the embryos at this stage on the left, the selected
          one previewed on the right. Pinned by a click on the tile. */}
      {open && selected && (
        <div
          className={`absolute bottom-full left-1/2 mb-2 -translate-x-1/2 p-3 ${glassPanel()}`}
        >
          <div className="mb-2 flex items-baseline gap-2">
            <span className="text-sm font-medium">{open.stage}</span>
            <span className="text-xs text-default-500">
              {open.members.length} embryo{open.members.length > 1 ? "s" : ""}
            </span>
            {pinned && (
              <button
                aria-label="Close"
                className="ml-auto text-xs text-default-500 hover:text-default-700"
                type="button"
                onClick={() => {
                  setPinned(null);
                  hover(null);
                }}
              >
                ✕
              </button>
            )}
          </div>

          <div className="flex gap-3">
            {/* Left: the embryos at this stage. */}
            <div className="flex w-48 max-h-52 flex-col gap-0.5 overflow-y-auto">
              {open.members.map((m) => (
                <button
                  key={m.id}
                  className={`flex items-center justify-between rounded px-2 py-1 text-left text-xs transition-colors ${
                    m.id === selected.id
                      ? "bg-primary/20 text-primary"
                      : "hover:bg-default-100"
                  }`}
                  type="button"
                  onClick={() => setSelectedId(m.id)}
                  onDoubleClick={() => onOpen(m)}
                >
                  <span className="truncate">{m.title}</span>
                  <span className="shrink-0 tabular-nums text-[10px] text-default-500">
                    {(m.numCells / 1e6).toFixed(1)}M
                  </span>
                </button>
              ))}
            </div>

            {/* Right: the selected embryo. */}
            <canvas
              ref={bigCanvasRef}
              className="h-52 w-52 shrink-0 rounded-lg bg-black/40"
            />
          </div>

          <div className="mt-2 flex gap-1">
            <Button
              className="flex-1"
              color="primary"
              size="sm"
              variant="flat"
              onPress={() => onOpen(selected)}
            >
              Open {selected.title}
            </Button>
            <Button size="sm" variant="ghost" onPress={() => onSplit(selected)}>
              Split screen
            </Button>
          </div>
        </div>
      )}

      <div
        ref={stripRef}
        className={`relative flex items-end gap-2 overflow-x-auto px-2 py-1.5 transition-[height] duration-150 [justify-content:safe_center] ${glassPanel()} rounded-none`}
        style={{
          height: hovered ? MAG_RAIL_HEIGHT : STAGE_RAIL_HEIGHT,
        }}
      >
        {/* One canvas for every tile; positioned under them, drawn per-viewport. */}
        <canvas
          ref={canvasRef}
          className="pointer-events-none absolute inset-0 h-full w-full"
        />

        {stages.map((s, i) => {
          const size = Math.round(
            TILE * (focus === -1 ? 1 : tileScale(i, focus)),
          );

          return (
            <button
              key={s.stage}
              className="relative z-10 shrink-0 cursor-pointer"
              style={{ width: size }}
              type="button"
              onClick={() => setPinned(pinned === s.stage ? null : s.stage)}
              // Always tracked, even while pinned: the card stays put because
              // `shown` prefers the pin, but magnification must follow the pointer.
              onMouseEnter={() => hover(s.stage)}
            >
              <div
                ref={(el) => setTileRef(s.stage, el)}
                className={`relative rounded-lg border-2 transition-[width,height,border-color] duration-150 ${
                  activeStage === s.stage
                    ? "border-primary"
                    : shown === s.stage
                      ? "border-default-400"
                      : "border-transparent"
                }`}
                style={{ width: size, height: size }}
              >
                {!loaded.has(s.rep.s3BaseUrl) && (
                  <div className="absolute inset-2 animate-pulse rounded bg-default-200/30" />
                )}
              </div>
              <div
                className={`mt-0.5 truncate text-center text-[10px] leading-3 ${
                  activeStage === s.stage ? "text-primary" : "text-default-500"
                }`}
              >
                {s.stage}
              </div>
            </button>
          );
        })}
      </div>
    </div>
  );
}
