"use client";

import type {
  CatalogDatasetEntry,
  CatalogDatasetItem,
  ExploreFilters,
} from "./types";

import { useEffect, useState } from "react";
import { createPortal } from "react-dom";
import { Button, Spinner } from "@heroui/react";
import { useSession } from "next-auth/react";

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
  internal: CatalogDatasetItem[];
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
  const { data: session } = useSession();
  const isAdmin =
    session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";

  const [data, setData] = useState<InitialData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // View state: list (the explore page) or a dataset's detail.
  const [detailDataset, setDetailDataset] = useState<CatalogDatasetItem | null>(
    null,
  );

  // Fetch initial data once when the modal first opens.
  useEffect(() => {
    if (!isOpen || data || loading) return;
    let cancelled = false;

    (async () => {
      setLoading(true);
      setError(null);
      try {
        const [resAll, resInternal] = await Promise.all([
          fetch("/api/explore?limit=50"),
          isAdmin
            ? fetch("/api/explore?tab=internal&limit=50")
            : Promise.resolve(null),
        ]);

        if (!resAll.ok) throw new Error("Failed to load catalog");
        const all = await resAll.json();
        const internal =
          resInternal && resInternal.ok
            ? ((await resInternal.json()).items ?? [])
            : [];

        if (cancelled) return;
        setData({
          items: all.items ?? [],
          total: all.total ?? 0,
          featured: all.featured ?? [],
          bil: all.bil ?? [],
          internal,
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
  }, [isOpen, data, loading, isAdmin]);

  // Reset detail view when the modal closes.
  useEffect(() => {
    if (!isOpen) setDetailDataset(null);
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
    <div className="fixed inset-0 z-[100] bg-background overflow-y-auto">
      {/* Close button — top right, fixed so it stays in view when scrolling */}
      <button
        aria-label="Close"
        className="fixed top-4 right-4 z-[101] w-9 h-9 rounded-full bg-default-100 hover:bg-default-200 text-default-700 flex items-center justify-center transition-colors shadow"
        type="button"
        onClick={onClose}
      >
        <svg
          className="w-5 h-5"
          fill="none"
          stroke="currentColor"
          strokeWidth={2}
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
              initialInternal={data.internal}
              initialItems={data.items}
              initialTotal={data.total}
              isAdmin={isAdmin}
              onCardClick={handleCardClick}
            />
          )
        ) : null}
      </div>
    </div>,
    document.body,
  );
}
