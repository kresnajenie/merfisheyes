"use client";

import type { CatalogDatasetItem } from "./types";

import { useCallback, useEffect, useRef, useState } from "react";

import { ExploreDatasetCard } from "./explore-dataset-card";

interface FeaturedDatasetsProps {
  datasets: CatalogDatasetItem[];
  onViewAll?: () => void;
  onCardClick?: (dataset: CatalogDatasetItem) => void;
}

/**
 * Featured row. This used to render a 4-up grid and hide the rest behind a
 * "view all" link, so anything past the fourth card was invisible. It is a
 * horizontal carousel instead: every featured dataset is reachable by
 * scrolling, and the arrows page through a viewport at a time.
 *
 * Cards use `usePopover` so their entry dropdown renders in a portal — an
 * inline dropdown would be clipped by the scroll container's overflow.
 */
export function FeaturedDatasets({
  datasets,
  onViewAll,
  onCardClick,
}: FeaturedDatasetsProps) {
  const scrollerRef = useRef<HTMLDivElement>(null);
  const [atStart, setAtStart] = useState(true);
  const [atEnd, setAtEnd] = useState(true);

  const sync = useCallback(() => {
    const el = scrollerRef.current;

    if (!el) return;
    // 1px slack: fractional widths mean scrollLeft rarely lands exactly on the
    // maximum, which would otherwise leave the "next" arrow permanently live.
    const max = el.scrollWidth - el.clientWidth;

    setAtStart(el.scrollLeft <= 1);
    setAtEnd(el.scrollLeft >= max - 1);
  }, []);

  useEffect(() => {
    const el = scrollerRef.current;

    if (!el) return;
    sync();
    const observer = new ResizeObserver(sync);

    observer.observe(el);

    return () => observer.disconnect();
  }, [sync, datasets.length]);

  const page = (direction: 1 | -1) => {
    const el = scrollerRef.current;

    if (!el) return;
    const reduced = window.matchMedia(
      "(prefers-reduced-motion: reduce)",
    ).matches;

    el.scrollBy({
      left: direction * el.clientWidth * 0.9,
      behavior: reduced ? "auto" : "smooth",
    });
  };

  if (datasets.length === 0) return null;

  const arrow =
    "grid place-items-center w-8 h-8 rounded-full border border-default-200 " +
    "bg-content1/80 backdrop-blur transition-opacity hover:bg-default-100 " +
    "disabled:opacity-30 disabled:cursor-default";

  return (
    <section>
      <div className="flex items-center justify-between mb-4 gap-3">
        <h2 className="text-lg font-semibold">
          Featured Datasets
          <span className="ml-2 text-sm font-normal text-default-500 tabular-nums">
            {datasets.length}
          </span>
        </h2>
        <div className="flex items-center gap-2">
          {onViewAll && (
            <button
              className="text-sm text-primary hover:underline"
              onClick={onViewAll}
            >
              View all
            </button>
          )}
          {/* Arrows are redundant when everything already fits. */}
          {!(atStart && atEnd) && (
            <div className="flex gap-1">
              <button
                aria-label="Scroll featured datasets left"
                className={arrow}
                disabled={atStart}
                type="button"
                onClick={() => page(-1)}
              >
                <span aria-hidden="true">←</span>
              </button>
              <button
                aria-label="Scroll featured datasets right"
                className={arrow}
                disabled={atEnd}
                type="button"
                onClick={() => page(1)}
              >
                <span aria-hidden="true">→</span>
              </button>
            </div>
          )}
        </div>
      </div>

      <div
        ref={scrollerRef}
        className="flex gap-4 overflow-x-auto snap-x snap-mandatory scroll-smooth
                   pb-2 -mx-1 px-1 [scrollbar-width:thin]"
        onScroll={sync}
      >
        {datasets.map((dataset) => (
          <div
            key={dataset.id}
            // Inline width rather than an arbitrary Tailwind class: `min()`
            // carries a comma, which is exactly the shape of arbitrary value
            // that silently fails to compile and collapses the card to 0px.
            className="snap-start shrink-0"
            style={{ width: "min(85vw, 272px)" }}
          >
            <ExploreDatasetCard
              usePopover
              dataset={dataset}
              onCardClick={onCardClick}
            />
          </div>
        ))}
      </div>
    </section>
  );
}
