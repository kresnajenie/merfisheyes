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
  // Read inside the render loop; as state it would re-run the effect below on
  // every hover, tearing down and re-fetching all previews.
  const hoveredRef = useRef<string | null>(null);
  const [previewsReady, setPreviewsReady] = useState(0);
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
        previewsRef.current.set(s.stage, p);
        setPreviewsReady((n) => n + 1);
      });
    }

    return () => {
      for (const c of controllers) c.abort();
      previewsRef.current.clear();
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

    for (const s of stages) {
      const p = previewsRef.current.get(s.stage);

      if (!p || scenes.has(s.stage)) continue;
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
      scenes.set(s.stage, { scene, camera, points });
    }

    let raf = 0;
    let last = performance.now();

    const frame = () => {
      raf = requestAnimationFrame(frame);
      const strip = stripRef.current;

      if (!strip) return;
      const now = performance.now();
      const dt = (now - last) / 1000;

      last = now;

      const stripRect = strip.getBoundingClientRect();

      if (
        canvas.width !== Math.floor(stripRect.width) ||
        canvas.height !== Math.floor(stripRect.height)
      ) {
        renderer.setSize(stripRect.width, stripRect.height, false);
      }

      renderer.setClearColor(0x000000, 0);
      renderer.clear();

      for (const s of stages) {
        const entry = scenes.get(s.stage);
        const el = tileRefs.current.get(s.stage);

        if (!entry || !el) continue;
        const r = el.getBoundingClientRect();

        // Cull tiles scrolled out of the strip.
        if (r.right < stripRect.left || r.left > stripRect.right) continue;

        // Slow spin so the preview reads as 3D; paused while hovered so a
        // drag-free look is possible.
        if (hoveredRef.current !== s.stage)
          entry.points.rotation.y += dt * 0.35;

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

  if (!stages?.length) return null;

  const open = hovered ? stages.find((s) => s.stage === hovered) : null;

  return (
    <div
      data-ui-overlay
      className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
      onMouseLeave={() => hover(null)}
    >
      {/* Hover card: what this stage is, and where it can go. */}
      {open && (
        <div
          className={`absolute bottom-full mb-2 w-72 p-3 ${glassPanel()}`}
          style={{
            left: Math.max(
              8,
              (tileRefs.current.get(open.stage)?.offsetLeft ?? 0) -
                (stripRef.current?.scrollLeft ?? 0) -
                40,
            ),
          }}
        >
          <div className="mb-1 text-sm font-medium">{open.stage}</div>
          <div className="mb-2 text-xs text-default-500">
            {open.members.length} embryo{open.members.length > 1 ? "s" : ""} ·
            showing {open.rep.title}
          </div>

          <div className="mb-2 flex gap-1">
            <Button
              className="flex-1"
              color="primary"
              size="sm"
              variant="flat"
              onPress={() => onOpen(open.rep)}
            >
              Open
            </Button>
            <Button size="sm" variant="ghost" onPress={() => onSplit(open.rep)}>
              Split screen
            </Button>
          </div>

          {/* A stage holds several embryos; pick one directly. */}
          {open.members.length > 1 && (
            <div className="max-h-40 overflow-y-auto">
              {open.members.map((m) => (
                <button
                  key={m.id}
                  className="flex w-full items-center justify-between rounded px-2 py-1 text-left text-xs hover:bg-default-100"
                  type="button"
                  onClick={() => onOpen(m)}
                >
                  <span className="truncate">{m.title}</span>
                  <span className="shrink-0 tabular-nums text-[10px] text-default-500">
                    {(m.numCells / 1e6).toFixed(1)}M
                  </span>
                </button>
              ))}
            </div>
          )}
        </div>
      )}

      <div
        ref={stripRef}
        className={`relative flex gap-2 overflow-x-auto px-3 py-2 ${glassPanel()} rounded-none`}
        style={{ height: TILE + 34 }}
      >
        {/* One canvas for every tile; positioned under them, drawn per-viewport. */}
        <canvas
          ref={canvasRef}
          className="pointer-events-none absolute inset-0 h-full w-full"
        />

        {stages.map((s) => (
          <button
            key={s.stage}
            className="relative z-10 shrink-0 cursor-pointer"
            type="button"
            onClick={() => onOpen(s.rep)}
            onMouseEnter={() => hover(s.stage)}
          >
            <div
              ref={(el) => setTileRef(s.stage, el)}
              className={`rounded-lg border-2 transition-colors ${
                activeStage === s.stage
                  ? "border-primary"
                  : hovered === s.stage
                    ? "border-default-400"
                    : "border-transparent"
              }`}
              style={{ width: TILE, height: TILE }}
            />
            <div
              className={`mt-0.5 text-center text-[10px] ${
                activeStage === s.stage ? "text-primary" : "text-default-500"
              }`}
            >
              {s.stage}
            </div>
          </button>
        ))}
      </div>
    </div>
  );
}
