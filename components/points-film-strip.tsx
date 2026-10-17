"use client";

import type { PointsElement } from "@/lib/spatialdata/points";

import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as THREE from "three";

import { glassButton, glassPanel } from "@/components/primitives";
import { usePanelSingleMoleculeVisualizationStore } from "@/lib/hooks/usePanelStores";
import {
  findPointsElements,
  openPointsElement,
} from "@/lib/spatialdata/points";
import { openSpatialDataStoreFromUrl } from "@/lib/spatialdata/store";
import { useFilmStripStore } from "@/lib/stores/filmStripStore";
import { STAGE_RAIL_HEIGHT, TILE } from "@/lib/ui/stage-rail-geometry";

export interface FilmStripDataset {
  id: string;
  title: string;
  numCells: number;
  s3BaseUrl: string;
  sortOrder: number;
}

// Molecules per gene per tile: a sparse level is plenty at thumbnail size.
const TILE_GENE_POINTS = 30_000;
// Points of the grey outline drawn per tile (all genes, coarsest level).
const OUTLINE_POINTS = 25_000;

const norm = (u: string) => u.replace(/\/+$/, "");

interface Tile {
  element: PointsElement;
  scene: THREE.Scene;
  camera: THREE.OrthographicCamera;
  genes: Map<string, { points: THREE.Points; level: number }>;
  loading: Set<string>;
  molecules: number;
}

/**
 * The sections of a project, one tile each, in project order — like the
 * lm-viewer's stage rail. Each tile draws its section's molecules: a grey
 * outline, and every gene selected in the main view, in its colour, from a
 * sparse level. Clicking a tile opens that section.
 *
 * All tiles share one WebGL context, drawn with per-tile scissor rectangles,
 * and are only redrawn when something changes.
 */
export function PointsFilmStrip({
  projectId,
  currentUrl,
  onOpen,
}: {
  projectId: string;
  currentUrl: string | null;
  onOpen: (d: FilmStripDataset) => void;
}) {
  const [project, setProject] = useState<{
    title: string;
    datasets: FilmStripDataset[];
  } | null>(null);
  const [ready, setReady] = useState<Set<string>>(new Set());
  const [busy, setBusy] = useState<Set<string>>(new Set());
  const [hovered, setHovered] = useState<string | null>(null);
  const tilesRef = useRef(new Map<string, Tile>());
  const tileEls = useRef(new Map<string, HTMLDivElement>());
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const stripRef = useRef<HTMLDivElement | null>(null);
  const dirtyRef = useRef(true);
  const selected = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenes,
  );
  const legend = usePanelSingleMoleculeVisualizationStore(
    (s) => s.selectedGenesLegend,
  );
  const geneCache = usePanelSingleMoleculeVisualizationStore(
    (s) => s.geneDataCache,
  );
  const setShown = useFilmStripStore((s) => s.setShown);

  useEffect(() => {
    setShown(true);

    return () => setShown(false);
  }, [setShown]);

  // ── The project's sections, in order
  useEffect(() => {
    let alive = true;

    fetch(`/api/projects/${projectId}/public`)
      .then((r) => (r.ok ? r.json() : null))
      .then((j) => {
        if (!alive || !j) return;
        setProject({
          title: j.title,
          datasets: (j.datasets as FilmStripDataset[])
            .filter((d) => d.s3BaseUrl)
            .sort((a, b) => a.sortOrder - b.sortOrder),
        });
      })
      .catch(() => {});

    return () => {
      alive = false;
    };
  }, [projectId]);

  // ── Each section's molecules: its element, and a grey outline
  useEffect(() => {
    if (!project) return;
    let alive = true;

    for (const d of project.datasets) {
      (async () => {
        const store = await openSpatialDataStoreFromUrl(d.s3BaseUrl);
        const [name] = await findPointsElements(store);

        if (!name || !alive) return;
        const element = await openPointsElement(store, name);
        const coarsest = element.levels.length - 1;
        const xyz = await element.readLevel(coarsest);

        if (!alive) return;
        const step = Math.max(1, Math.ceil(xyz.length / 3 / OUTLINE_POINTS));
        const thin = new Float32Array(Math.ceil(xyz.length / 3 / step) * 3);

        for (let i = 0, j = 0; i < xyz.length; i += 3 * step, j += 3) {
          thin[j] = xyz[i];
          thin[j + 1] = xyz[i + 1];
        }
        const scene = new THREE.Scene();
        const outline = new THREE.Points(
          new THREE.BufferGeometry().setAttribute(
            "position",
            new THREE.BufferAttribute(thin, 3),
          ),
          new THREE.PointsMaterial({
            color: 0x8a8f98,
            size: 1,
            sizeAttenuation: false,
            transparent: true,
            opacity: 0.18,
          }),
        );

        scene.add(outline);
        const camera = frame(element, xyz);
        const counts = await element.featureCounts();

        tilesRef.current.set(d.s3BaseUrl, {
          element,
          scene,
          camera,
          genes: new Map(),
          loading: new Set(),
          molecules: counts.reduce((a, b) => a + b, 0),
        });
        dirtyRef.current = true;
        setReady((prev) => new Set(prev).add(d.s3BaseUrl));
      })().catch((e) => console.warn(`Film strip: ${d.title}:`, e));
    }

    return () => {
      alive = false;
      for (const t of tilesRef.current.values()) {
        t.scene.traverse((o) => {
          if (o instanceof THREE.Points) {
            o.geometry.dispose();
            (o.material as THREE.Material).dispose();
          }
        });
      }
      tilesRef.current.clear();
      setReady(new Set());
    };
  }, [project]);

  // ── The selected genes, in every tile, from a sparse level
  useEffect(() => {
    // Hidden genes keep their points, just not drawn, so showing them again
    // fetches nothing; a gene is let go only when it is deselected
    // The legend's genes: a hidden one (eye, or another gene's solo) is out
    // of the selection but still in the legend
    const viz = (g: string) => selected.get(g) ?? geneCache.get(g);
    const shown = (g: string) =>
      selected.has(g) && !!selected.get(g)?.showAssigned;
    const wanted = new Map(
      Array.from(legend)
        .filter((g) => viz(g))
        .map((g) => [g, viz(g)!.color]),
    );

    for (const [url, tile] of tilesRef.current) {
      for (const [gene, entry] of tile.genes) {
        const color = wanted.get(gene);

        if (color === undefined) {
          tile.scene.remove(entry.points);
          entry.points.geometry.dispose();
          (entry.points.material as THREE.Material).dispose();
          tile.genes.delete(gene);
        } else {
          (entry.points.material as THREE.PointsMaterial).color.set(color);
          entry.points.visible = shown(gene);
        }
      }
      for (const gene of wanted.keys()) {
        if (!shown(gene)) continue;
        if (tile.genes.has(gene) || tile.loading.has(gene)) continue;
        const f = tile.element.features.indexOf(gene);

        if (f < 0) continue;
        tile.loading.add(gene);
        setBusy((prev) => new Set(prev).add(url));
        loadGene(tile, f)
          .then((got) => {
            const color = viz(gene)?.color;

            if (
              !got ||
              color === undefined ||
              tilesRef.current.get(url) !== tile
            )
              return;
            const points = new THREE.Points(
              new THREE.BufferGeometry().setAttribute(
                "position",
                new THREE.BufferAttribute(got.xyz, 3),
              ),
              new THREE.PointsMaterial({
                color,
                size: 2,
                sizeAttenuation: false,
              }),
            );

            tile.scene.add(points);
            tile.genes.set(gene, { points, level: got.level });
          })
          .catch((e) => console.warn(`Film strip: ${gene}:`, e))
          .finally(() => {
            tile.loading.delete(gene);
            dirtyRef.current = true;
            if (tile.loading.size === 0) {
              setBusy((prev) => {
                const next = new Set(prev);

                next.delete(url);

                return next;
              });
            }
          });
      }
    }
    dirtyRef.current = true;
  }, [selected, legend, geneCache, ready]);

  // ── One renderer for every tile, redrawn when something changed
  useEffect(() => {
    const canvas = canvasRef.current;
    const strip = stripRef.current;

    if (!canvas || !strip || !project?.datasets.length) return;
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
    const markDirty = () => {
      dirtyRef.current = true;
    };

    strip.addEventListener("scroll", markDirty, { passive: true });
    window.addEventListener("resize", markDirty);
    const ro = new ResizeObserver(markDirty);

    ro.observe(strip);
    let raf = 0;
    const draw = () => {
      raf = requestAnimationFrame(draw);
      if (!dirtyRef.current) return;
      dirtyRef.current = false;
      const sr = strip.getBoundingClientRect();

      renderer.setSize(sr.width, sr.height, false);
      renderer.setScissor(0, 0, sr.width, sr.height);
      renderer.setClearColor(0x000000, 0);
      renderer.clear();
      for (const d of project.datasets) {
        const tile = tilesRef.current.get(d.s3BaseUrl);
        const el = tileEls.current.get(d.id);

        if (!tile || !el) continue;
        const r = el.getBoundingClientRect();

        if (r.right < sr.left || r.left > sr.right) continue;
        const x = r.left - sr.left;
        const y = sr.bottom - r.bottom; // WebGL origin is bottom-left

        renderer.setViewport(x, y, r.width, r.height);
        renderer.setScissor(x, y, r.width, r.height);
        renderer.render(tile.scene, tile.camera);
      }
    };

    draw();

    return () => {
      cancelAnimationFrame(raf);
      strip.removeEventListener("scroll", markDirty);
      window.removeEventListener("resize", markDirty);
      ro.disconnect();
      renderer.dispose();
    };
  }, [project]);

  const setTileEl = useCallback((id: string, el: HTMLDivElement | null) => {
    if (el) tileEls.current.set(id, el);
    else tileEls.current.delete(id);
  }, []);

  const activeId = useMemo(
    () =>
      currentUrl && project
        ? (project.datasets.find((d) => norm(d.s3BaseUrl) === norm(currentUrl))
            ?.id ?? null)
        : null,
    [currentUrl, project],
  );
  const hover = project?.datasets.find((d) => d.id === hovered) ?? null;
  const hoverTile = hover ? tilesRef.current.get(hover.s3BaseUrl) : null;

  if (!project?.datasets.length) return null;

  return (
    <div
      data-ui-overlay
      className="absolute bottom-0 left-0 right-0 z-[var(--z-rail)]"
      data-testid="points-film-strip"
      onMouseLeave={() => setHovered(null)}
    >
      {hover && (
        <div
          className={`absolute bottom-full left-1/2 mb-2 -translate-x-1/2 whitespace-nowrap px-3 py-1.5 text-xs ${glassButton()}`}
        >
          <span className="font-medium">{hover.title}</span>
          <span className="text-default-500">
            {" · "}
            {hover.numCells.toLocaleString()} cells
            {hoverTile &&
              ` · ${(hoverTile.molecules / 1e6).toFixed(0)} M molecules`}
            {hoverTile &&
              Array.from(hoverTile.genes.entries())
                .map(
                  ([g, e]) =>
                    ` · ${g} ${e.level ? `1 in ${4 ** e.level}` : "all"}`,
                )
                .join("")}
          </span>
        </div>
      )}
      <div
        className={`relative mx-auto w-fit max-w-full ${glassPanel()} rounded-b-none`}
        style={{ height: STAGE_RAIL_HEIGHT }}
      >
        <canvas
          ref={canvasRef}
          className="pointer-events-none absolute inset-0 h-full w-full"
        />
        <div
          ref={stripRef}
          className="relative z-10 flex h-full flex-col gap-1 overflow-x-auto px-2 py-1.5"
        >
          <div className="truncate px-1 text-[10px] leading-3 text-default-500">
            {project.title}
          </div>
          <div className="flex gap-2">
            {project.datasets.map((d) => (
              <button
                key={d.id}
                className="relative shrink-0 cursor-pointer"
                style={{ width: TILE }}
                type="button"
                onClick={() => onOpen(d)}
                onMouseEnter={() => setHovered(d.id)}
              >
                <div
                  ref={(el) => setTileEl(d.id, el)}
                  className={`relative rounded-lg border-2 transition-colors ${
                    activeId === d.id
                      ? "border-primary"
                      : hovered === d.id
                        ? "border-default-400"
                        : "border-transparent"
                  }`}
                  style={{ width: TILE, height: TILE }}
                >
                  {!ready.has(d.s3BaseUrl) && (
                    <div className="absolute inset-2 animate-pulse rounded bg-default-200/30" />
                  )}
                  {busy.has(d.s3BaseUrl) && (
                    <span className="absolute right-1 top-1 h-2.5 w-2.5 animate-spin rounded-full border-2 border-white/25 border-t-white" />
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
      </div>
    </div>
  );
}

/** A camera framing the whole section in a square tile. */
function frame(element: PointsElement, xyz: Float32Array) {
  let [[x0, y0], [x1, y1]] = element.extent ?? [
    [Infinity, Infinity],
    [-Infinity, -Infinity],
  ];

  if (!element.extent) {
    for (let i = 0; i < xyz.length; i += 3) {
      x0 = Math.min(x0, xyz[i]);
      x1 = Math.max(x1, xyz[i]);
      y0 = Math.min(y0, xyz[i + 1]);
      y1 = Math.max(y1, xyz[i + 1]);
    }
  }
  const cx = (x0 + x1) / 2;
  const cy = (y0 + y1) / 2;
  const half = (Math.max(x1 - x0, y1 - y0) / 2) * 1.05;
  const camera = new THREE.OrthographicCamera(
    cx - half,
    cx + half,
    cy + half,
    cy - half,
    -1000,
    1000,
  );

  camera.position.set(0, 0, 10);

  return camera;
}

/** A gene's molecules in one section, from the finest level that is sparse enough. */
async function loadGene(tile: Tile, f: number) {
  const { element } = tile;

  for (let level = 0; level < element.levels.length; level++) {
    const counts = await element.featureCounts(level);

    if (counts[f] > TILE_GENE_POINTS && level < element.levels.length - 1) {
      continue;
    }
    if (counts[f] === 0) return null;
    const ft = await element.featureTiles(level, f);
    const lo = ft.starts[0];

    return { level, xyz: await element.readXyz(level, lo, lo + counts[f]) };
  }

  return null;
}
