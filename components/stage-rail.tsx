"use client";

import type { DatasetPreview } from "@/lib/webgl/preview";

import { Button } from "@heroui/button";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";

import { glassButton, glassPanel } from "@/components/primitives";
import { STAGE_RAIL_HEIGHT, TILE } from "@/lib/ui/stage-rail-geometry";
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

interface Props {
  projectId: string;
  /** s3BaseUrl of the dataset currently open, to mark its tile. */
  currentUrl?: string | null;
  onOpen: (d: ProjectDatasetSummary) => void;
  onSplit: (d: ProjectDatasetSummary) => void;
}

// Re-exported so existing importers are unaffected; the numbers live in their
// own module because the footer needs them without pulling three.js in.
export { STAGE_RAIL_HEIGHT };

/** Width of the hover card, needed to clamp it inside the strip. */
const CARD_W = 260;

const norm = (u: string) => u.replace(/\/+$/, "");

/**
 * Allen-atlas style rail: one tile per embryo, in developmental order, each a
 * live preview rather than a screenshot. Clicking a tile opens it.
 *
 * Every embryo is its own tile. Grouping by stage put two clicks between the
 * user and a dataset the strip was already showing — one to open the stage,
 * another to pick an embryo inside it.
 *
 * All tiles share ONE WebGL context, drawn with per-tile scissor rectangles.
 * A context each would exhaust the browser's limit (~16) and cost far more
 * memory; one canvas over the strip, with viewports recomputed from the tiles'
 * DOM rects, also means scrolling needs no special handling.
 */
export default function StageRail({
  projectId,
  currentUrl,
  onOpen,
  onSplit,
}: Props) {
  const [items, setItems] = useState<ProjectDatasetSummary[] | null>(null);
  const [hovered, setHovered] = useState<string | null>(null);
  /** Centre of the hovered tile, relative to the strip, to anchor the card. */
  const [cardX, setCardX] = useState(0);

  const [previewsReady, setPreviewsReady] = useState(0);
  const [loaded, setLoaded] = useState<Set<string>>(new Set());
  /** Keyed by s3BaseUrl. */
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

  // ── Project members, in project order (which is developmental order).
  useEffect(() => {
    let alive = true;

    fetch(`/api/projects/${projectId}/public`)
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => {
        if (!alive || !j) return;
        const ds = (j.datasets as ProjectDatasetSummary[])
          .filter((d) => d.s3BaseUrl)
          .sort((a, b) => a.sortOrder - b.sortOrder);

        setItems(ds);
      })
      .catch(() => {});

    return () => {
      alive = false;
    };
  }, [projectId]);

  // ── Preview data. Independent of the renderer, so a browser without WebGL
  //    fetches nothing it cannot use and degrades to plain labelled buttons.
  useEffect(() => {
    if (!items?.length) return;
    const controllers: AbortController[] = [];

    for (const d of items) {
      const ac = new AbortController();

      controllers.push(ac);
      loadPreview(d.s3BaseUrl, ac.signal).then((p) => {
        if (!p) return;
        previewsRef.current.set(d.s3BaseUrl, p);
        setLoaded((prev) => new Set(prev).add(d.s3BaseUrl));
        setPreviewsReady((n) => n + 1);
      });
    }

    return () => {
      for (const c of controllers) c.abort();
      previewsRef.current.clear();
      setLoaded(new Set());
    };
  }, [items]);

  /**
   * Scene for one preview, built on first use and kept until the project
   * changes.
   *
   * Built lazily rather than in a pass keyed on `previewsReady`: that rebuilt
   * every scene whenever one more preview arrived, resetting all the tiles'
   * rotations.
   */
  const sceneFor = useCallback((url: string) => {
    const existing = scenesRef.current.get(url);

    if (existing) return existing;

    const p = previewsRef.current.get(url);

    if (!p) return null;

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

    const entry = { scene, camera, points };

    scenesRef.current.set(url, entry);

    return entry;
  }, []);

  // ── One renderer for every tile.
  useEffect(() => {
    const canvas = canvasRef.current;

    if (!canvas || !items?.length) return;

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

    const strip = stripRef.current;

    if (!strip) return;

    // Tile rectangles are measured once and reused. Reading them every frame
    // meant a forced layout flush per tile per frame, and with 45 embryos that
    // is 45 of them.
    /** Last CSS size the drawing buffer was sized for. */
    const sized = { w: 0, h: 0 };
    let layoutDirty = true;
    let layout: { url: string; x: number; y: number; w: number; h: number }[] =
      [];

    const measure = () => {
      const sr = strip.getBoundingClientRect();

      if (sized.w !== sr.width || sized.h !== sr.height) {
        sized.w = sr.width;
        sized.h = sr.height;
        renderer.setSize(sr.width, sr.height, false);
      }

      layout = [];
      for (const d of items) {
        const el = tileRefs.current.get(d.id);

        if (!el) continue;
        const r = el.getBoundingClientRect();

        // Cull tiles scrolled out of the strip. With 45 embryos most are off
        // screen at any moment, so this is what keeps the cost flat.
        if (r.right < sr.left || r.left > sr.right) continue;
        layout.push({
          url: d.s3BaseUrl,
          x: r.left - sr.left,
          y: sr.bottom - r.bottom, // WebGL origin is bottom-left
          w: r.width,
          h: r.height,
        });
      }
      layoutDirty = false;
    };

    const markDirty = () => {
      layoutDirty = true;
    };

    strip.addEventListener("scroll", markDirty, { passive: true });
    window.addEventListener("resize", markDirty);

    const ro = new ResizeObserver(markDirty);

    ro.observe(strip);

    // 30fps: the spin is slow enough that the extra 30 frames buy nothing, and
    // this halves what the rail takes from the main scene.
    const FRAME_MS = 1000 / 30;
    let raf = 0;
    let last = performance.now();

    const frame = () => {
      raf = requestAnimationFrame(frame);

      const now = performance.now();

      if (now - last < FRAME_MS) return;
      // Clamped: a backgrounded tab resumes with a huge gap and would snap
      // every tile to a new angle.
      const dt = Math.min((now - last) / 1000, 0.1);

      last = now;

      if (layoutDirty) measure();

      renderer.setClearColor(0x000000, 0);
      renderer.clear();

      for (const t of layout) {
        const entry = sceneFor(t.url);

        if (!entry) continue;

        // Slow spin so the preview reads as 3D.
        entry.points.rotation.y += dt * 0.35;

        renderer.setViewport(t.x, t.y, t.w, t.h);
        renderer.setScissor(t.x, t.y, t.w, t.h);
        entry.camera.aspect = t.w / t.h;
        entry.camera.updateProjectionMatrix();
        renderer.render(entry.scene, entry.camera);
      }
    };

    frame();

    return () => {
      cancelAnimationFrame(raf);
      strip.removeEventListener("scroll", markDirty);
      window.removeEventListener("resize", markDirty);
      ro.disconnect();
      for (const { scene, points } of scenesRef.current.values()) {
        points.geometry.dispose();
        (points.material as THREE.Material).dispose();
        scene.clear();
      }
      scenesRef.current.clear();
      renderer.dispose();
    };
  }, [items, sceneFor]);

  const open = useMemo(
    () => items?.find((d) => d.id === hovered) ?? null,
    [items, hovered],
  );

  // ── Enlarged preview inside the hover card.
  //
  // Its own context: the strip's canvas lives inside the strip and cannot
  // reach the card above it. Two contexts is still far under the browser's
  // limit and much simpler than reshaping the shared canvas to span both.
  useEffect(() => {
    const canvas = bigCanvasRef.current;
    const entry = open ? sceneFor(open.s3BaseUrl) : null;

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
  }, [open, previewsReady, sceneFor]);

  const setTileRef = useCallback((id: string, el: HTMLDivElement | null) => {
    if (el) tileRefs.current.set(id, el);
    else tileRefs.current.delete(id);
  }, []);

  /** The embryo currently open, for the title and the active ring. */
  const activeId = useMemo(() => {
    if (!currentUrl || !items) return null;

    return (
      items.find((d) => norm(d.s3BaseUrl) === norm(currentUrl))?.id ?? null
    );
  }, [currentUrl, items]);

  const active = items?.find((d) => d.id === activeId) ?? null;

  /** Anchor the card over the tile the pointer is on, clamped to the strip. */
  const hover = useCallback((d: ProjectDatasetSummary | null) => {
    setHovered(d?.id ?? null);
    if (!d) return;

    const el = tileRefs.current.get(d.id);
    const strip = stripRef.current;

    if (!el || !strip) return;
    const r = el.getBoundingClientRect();
    const sr = strip.getBoundingClientRect();
    const centre = r.left + r.width / 2 - sr.left;

    setCardX(
      Math.max(CARD_W / 2 + 8, Math.min(sr.width - CARD_W / 2 - 8, centre)),
    );
  }, []);

  if (!items) {
    return (
      <div
        data-ui-overlay
        className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
      >
        <div
          className={`mx-auto flex w-fit max-w-full items-end gap-2 overflow-hidden px-2 py-1.5 ${glassPanel()} rounded-b-none`}
          style={{ height: STAGE_RAIL_HEIGHT }}
        >
          {Array.from({ length: 10 }).map((_, i) => (
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

  if (!items.length) return null;

  return (
    <>
      {/* Which embryo is on screen. The rail shows where it sits in the
          series; this says what it actually is. */}
      {active && (
        <div
          data-ui-overlay
          className={`absolute left-1/2 top-6 z-[var(--z-legends)] -translate-x-1/2 whitespace-nowrap rounded-full px-4 py-1.5 ${glassButton()}`}
        >
          <span className="text-sm font-medium">{active.stage}</span>
          <span className="text-sm text-default-500">
            {" · "}
            {active.title}
          </span>
        </div>
      )}

      <div
        data-ui-overlay
        className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
        onMouseLeave={() => hover(null)}
      >
        {/* Hover card, anchored over its tile. Clicking the tile opens the
            embryo, so this only has to name it and offer the split. */}
        {open && (
          <div
            className={`absolute bottom-full mb-2 -translate-x-1/2 p-3 ${glassPanel()}`}
            style={{ left: cardX, width: CARD_W }}
          >
            <canvas
              ref={bigCanvasRef}
              className="mb-2 h-40 w-full rounded-lg bg-black/40"
            />
            <div className="text-sm font-medium">{open.title}</div>
            <div className="mb-2 text-xs text-default-500">
              {open.stage} · {(open.numCells / 1e6).toFixed(1)}M molecules
            </div>
            <Button
              className="w-full"
              size="sm"
              variant="ghost"
              onPress={() => {
                onSplit(open);
                hover(null);
              }}
            >
              Split screen
            </Button>
          </div>
        )}

        <div
          ref={stripRef}
          // Hugs its tiles, and `max-w-full` keeps the scroll for narrow
          // windows — with 45 embryos it is always scrolling.
          className={`relative mx-auto flex w-fit max-w-full items-end gap-2 overflow-x-auto px-2 py-1.5 ${glassPanel()} rounded-b-none`}
          style={{ height: STAGE_RAIL_HEIGHT }}
        >
          {/* One canvas for every tile; under them, drawn per-viewport. */}
          <canvas
            ref={canvasRef}
            className="pointer-events-none absolute inset-0 h-full w-full"
          />

          {items.map((d) => (
            <button
              key={d.id}
              className="relative z-10 shrink-0 cursor-pointer"
              style={{ width: TILE }}
              type="button"
              onClick={() => onOpen(d)}
              onMouseEnter={() => hover(d)}
            >
              <div
                ref={(el) => setTileRef(d.id, el)}
                className={`relative rounded-lg border-2 transition-colors ${
                  activeId === d.id
                    ? "border-primary"
                    : hovered === d.id
                      ? "border-default-400"
                      : "border-transparent"
                }`}
                style={{ width: TILE, height: TILE }}
              >
                {!loaded.has(d.s3BaseUrl) && (
                  <div className="absolute inset-2 animate-pulse rounded bg-default-200/30" />
                )}
              </div>
              <div
                className={`mt-0.5 truncate text-center text-[10px] leading-3 ${
                  activeId === d.id ? "text-primary" : "text-default-500"
                }`}
              >
                {d.title}
              </div>
            </button>
          ))}
        </div>
      </div>
    </>
  );
}
