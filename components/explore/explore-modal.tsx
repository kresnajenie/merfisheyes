"use client";

import type {
  CatalogDatasetEntry,
  CatalogDatasetItem,
  ExploreFilters,
} from "./types";

import { useEffect, useRef, useState } from "react";
import { createPortal } from "react-dom";
import { Button, Spinner } from "@heroui/react";

import { ExplorePageClient } from "./explore-page-client";
import { DatasetDetailPage } from "./dataset-detail-page";

interface ExploreModalProps {
  isOpen: boolean;
  onClose: () => void;
  /** Called when the user finally picks a viewable entry. Modal closes after. */
  onSelectEntry: (
    dataset: CatalogDatasetItem,
    entry: CatalogDatasetEntry,
  ) => void;
}

interface InitialData {
  items: CatalogDatasetItem[];
  total: number;
  featured: CatalogDatasetItem[];
  bil: CatalogDatasetItem[];
  filters: ExploreFilters;
}

/**
 * Full-viewport modal that renders the /explore page inside the split-screen
 * viewer. Card clicks are intercepted: single-entry datasets load straight
 * into the right panel; multi-entry / BIL cards switch into an in-modal
 * detail view (same component used by /explore/[id]) where the user picks
 * the specific entry.
 */
export function ExploreModal({
  isOpen,
  onClose,
  onSelectEntry,
}: ExploreModalProps) {
  const [data, setData] = useState<InitialData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // View state: list (the explore page) or a dataset's detail.
  const [detailDataset, setDetailDataset] = useState<CatalogDatasetItem | null>(
    null,
  );

  // Single-shot guard. Using a ref instead of the `loading`/`data` state
  // avoids the bug where `setLoading(true)` re-runs the effect, runs the
  // cleanup (`cancelled = true`), and the still-pending fetch silently
  // bails in its finally block — leaving `loading` stuck true forever.
  const hasFetchedRef = useRef(false);

  // Fetch initial data once when the modal first opens. Internal (admin-only)
  // datasets are lazy-loaded by ExplorePageClient when the Internal tab opens,
  // so there's nothing admin-specific to fetch here.
  useEffect(() => {
    if (!isOpen || hasFetchedRef.current) {
      return;
    }
    hasFetchedRef.current = true;
    let cancelled = false;

    (async () => {
      setLoading(true);
      setError(null);
      try {
        const resAll = await fetch("/api/explore?limit=50");

        if (!resAll.ok) throw new Error("Failed to load catalog");
        const all = await resAll.json();

        if (cancelled) return;
        setData({
          items: all.items ?? [],
          total: all.total ?? 0,
          featured: all.featured ?? [],
          bil: all.bil ?? [],
          filters: all.filters ?? { species: [], tissues: [], platforms: [] },
        });
      } catch (e) {
        if (!cancelled) setError((e as Error).message);
      } finally {
        if (!cancelled) setLoading(false);
      }
    })();

    return () => {
      cancelled = true;
    };
  }, [isOpen]);

  // Reset detail view + fetch guard when the modal closes so the next open
  // refetches.
  useEffect(() => {
    if (!isOpen) {
      setDetailDataset(null);
      hasFetchedRef.current = false;
    }
  }, [isOpen]);

  // ESC to close.
  useEffect(() => {
    if (!isOpen) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") {
        if (detailDataset) setDetailDataset(null);
        else onClose();
      }
    };

    window.addEventListener("keydown", onKey);

    return () => window.removeEventListener("keydown", onKey);
  }, [isOpen, detailDataset, onClose]);

  if (!isOpen) return null;
  if (typeof document === "undefined") return null;

  // List view → switch to detail for multi-entry / BIL; load directly when
  // exactly one viewable entry exists.
  const handleCardClick = (dataset: CatalogDatasetItem) => {
    const entries = (dataset.entries ?? []).filter(
      (e) => e.s3BaseUrl || e.datasetId,
    );
    const isMulti = entries.length > 1 || !!dataset.bilCode;

    if (isMulti) {
      setDetailDataset(dataset);

      return;
    }
    if (entries.length === 1) {
      onSelectEntry(dataset, entries[0]);
      onClose();

      return;
    }
    // No viewable entries — fall back to opening detail (where the user can
    // at least see metadata + any external links).
    setDetailDataset(dataset);
  };

  const handleDetailEntrySelect = (entry: CatalogDatasetEntry) => {
    if (!detailDataset) return;
    onSelectEntry(detailDataset, entry);
    onClose();
  };

  return createPortal(
    // Outer layer — fills the viewport so clicks anywhere outside the card
    // close the modal. Backdrop-blur softens the viewer behind without
    // darkening it; the viewer's colors still show through the blur.
    <div
      className="fixed inset-0 z-[100] flex items-center justify-center p-6 sm:p-12 md:p-16 lg:p-20 backdrop-blur-md"
      role="presentation"
      onClick={onClose}
    >
      {/* Inner card — capped width, centered, bordered, rounded, own scroll.
          Stops click propagation so interacting inside doesn't dismiss. */}
      <div
        className="relative w-full h-full max-w-7xl rounded-2xl bg-background shadow-2xl border border-default-200 overflow-y-auto"
        role="dialog"
        onClick={(e) => e.stopPropagation()}
      >
        {/* Close button — sits inside the card, top-right corner. Solid
            danger color for prominence. */}
        <button
          aria-label="Close"
          className="absolute top-4 right-4 z-[101] w-10 h-10 rounded-full bg-danger hover:bg-danger-600 text-white flex items-center justify-center transition-colors shadow-lg"
          type="button"
          onClick={onClose}
        >
          <svg
            className="w-5 h-5"
            fill="none"
            stroke="currentColor"
            strokeWidth={2.5}
            viewBox="0 0 24 24"
          >
            <path
              d="M6 18L18 6M6 6l12 12"
              strokeLinecap="round"
              strokeLinejoin="round"
            />
          </svg>
        </button>

        <div className="container mx-auto max-w-7xl px-6 py-8">
        {loading && !data ? (
          <div className="flex items-center justify-center py-24">
            <Spinner size="lg" />
          </div>
        ) : error ? (
          <div className="flex flex-col items-center gap-3 py-24">
            <p className="text-danger">{error}</p>
            <Button
              size="sm"
              variant="bordered"
              onPress={() => {
                setData(null);
                setError(null);
              }}
            >
              Retry
            </Button>
          </div>
        ) : data ? (
          detailDataset ? (
            <DatasetDetailPage
              dataset={detailDataset}
              onBack={() => setDetailDataset(null)}
              onSelectEntry={handleDetailEntrySelect}
            />
          ) : (
            <ExplorePageClient
              disableUrlSync
              initialBil={data.bil}
              initialFeatured={data.featured}
              initialFilters={data.filters}
              initialItems={data.items}
              initialTotal={data.total}
              onCardClick={handleCardClick}
            />
          )
        ) : null}
        </div>
      </div>
    </div>,
    document.body,
  );
}
