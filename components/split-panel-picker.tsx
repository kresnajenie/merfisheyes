"use client";

import type { PanelType } from "@/lib/stores/splitScreenStore";
import type {
  CatalogDatasetEntry,
  CatalogDatasetItem,
} from "@/components/explore/types";

import { useCallback, useEffect, useRef, useState } from "react";
import { Button, Input, Skeleton } from "@heroui/react";
import { usePathname, useSearchParams } from "next/navigation";

import { useSplitScreenStore } from "@/lib/stores/splitScreenStore";
import { ExploreDatasetCard } from "@/components/explore/explore-dataset-card";

function SkeletonCard() {
  return (
    <div className="rounded-xl bg-default-100/50 overflow-hidden">
      <Skeleton className="h-32 w-full rounded-none" />
      <div className="p-4 flex flex-col gap-2">
        <div className="flex gap-1">
          <Skeleton className="h-5 w-8 rounded-full" />
        </div>
        <Skeleton className="h-4 w-3/4 rounded-md" />
        <Skeleton className="h-3 w-full rounded-md" />
        <Skeleton className="h-3 w-1/2 rounded-md" />
      </div>
    </div>
  );
}

function SkeletonGrid({ count = 4 }: { count?: number }) {
  return (
    <div className="grid grid-cols-2 gap-3">
      {Array.from({ length: count }).map((_, i) => (
        <SkeletonCard key={i} />
      ))}
    </div>
  );
}

type PickerView = "main" | "link" | "catalog";

export function SplitPanelPicker() {
  const { setRightPanel, setRightPanelS3 } = useSplitScreenStore();
  const pathname = usePathname();
  const searchParams = useSearchParams();
  const [linkInput, setLinkInput] = useState("");
  const [view, setView] = useState<PickerView>("main");

  // Catalog state
  const [catalogItems, setCatalogItems] = useState<CatalogDatasetItem[]>([]);
  const [featuredItems, setFeaturedItems] = useState<CatalogDatasetItem[]>([]);
  const [catalogSearch, setCatalogSearch] = useState("");
  const [catalogLoading, setCatalogLoading] = useState(false);
  const [featuredLoaded, setFeaturedLoaded] = useState(false);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const isFromS3Page = pathname.includes("/from-s3");

  // Extract current dataset info from URL
  const getCurrentDatasetInfo = (): {
    id: string | null;
    s3Url: string | null;
    type: PanelType;
  } => {
    const isSm = pathname.startsWith("/sm-viewer");
    const type: PanelType = isSm ? "sm" : "cell";

    // For from-s3 pages, return the S3 URL
    if (isFromS3Page) {
      const s3Url = searchParams.get("url");

      return {
        id: null,
        s3Url: s3Url ? decodeURIComponent(s3Url) : null,
        type,
      };
    }

    // Extract ID from /viewer/[id] or /sm-viewer/[id]
    const parts = pathname.split("/");
    const id = parts.length >= 3 ? parts[parts.length - 1] : null;

    // Don't return ID for non-ID pages like /viewer or /sm-viewer
    if (id === "viewer" || id === "sm-viewer") {
      return { id: null, s3Url: null, type };
    }

    return { id, s3Url: null, type };
  };

  const handleSameDataset = () => {
    const { id, s3Url, type } = getCurrentDatasetInfo();

    if (s3Url) {
      setRightPanelS3(s3Url, type);
    } else if (id) {
      setRightPanel(id, type);
    }
  };

  // Detect if a URL is a raw S3 URL (e.g., https://bucket.s3.region.amazonaws.com/...)
  const isRawS3Url = (url: string): boolean => {
    try {
      const parsed = new URL(url);

      return (
        parsed.hostname.includes("s3") &&
        parsed.hostname.includes("amazonaws.com")
      );
    } catch {
      return false;
    }
  };

  const handlePasteLink = () => {
    if (!linkInput.trim()) return;

    const input = linkInput.trim();

    // Check if it's a raw S3 URL
    if (isRawS3Url(input)) {
      const isSm = pathname.startsWith("/sm-viewer");

      setRightPanelS3(input, isSm ? "sm" : "cell");

      return;
    }

    try {
      const url = new URL(input);
      const pathParts = url.pathname.split("/").filter(Boolean);

      let type: PanelType = "cell";
      let id: string | null = null;

      if (pathParts[0] === "sm-viewer") {
        type = "sm";
        if (pathParts[1] === "from-s3") {
          const s3UrlParam = url.searchParams.get("url");

          if (s3UrlParam) {
            setRightPanelS3(decodeURIComponent(s3UrlParam), type);

            return;
          }
        } else if (pathParts.length >= 2) {
          id = pathParts[1];
        }
      } else if (pathParts[0] === "viewer") {
        type = "cell";
        if (pathParts[1] === "from-s3") {
          const s3UrlParam = url.searchParams.get("url");

          if (s3UrlParam) {
            setRightPanelS3(decodeURIComponent(s3UrlParam), type);

            return;
          }
        } else if (pathParts.length >= 2) {
          id = pathParts[1];
        }
      }

      if (id) {
        setRightPanel(id, type);
      }
    } catch {
      // Invalid URL — try treating as dataset ID
      if (input.length > 0) {
        const isSm = pathname.startsWith("/sm-viewer");

        setRightPanel(input, isSm ? "sm" : "cell");
      }
    }
  };

  // Fetch featured datasets on mount
  useEffect(() => {
    if (featuredLoaded) return;
    (async () => {
      try {
        const res = await fetch("/api/explore?limit=1");

        if (!res.ok) return;
        const data = await res.json();

        setFeaturedItems(data.featured ?? []);
      } catch {
        // ignore
      } finally {
        setFeaturedLoaded(true);
      }
    })();
  }, [featuredLoaded]);

  // Fetch catalog datasets (for catalog view search)
  const fetchCatalog = useCallback(async (search: string) => {
    setCatalogLoading(true);
    try {
      const params = new URLSearchParams({ limit: "50" });

      if (search) params.set("search", search);
      const res = await fetch(`/api/explore?${params}`);

      if (!res.ok) throw new Error("Failed to fetch catalog");
      const data = await res.json();

      setCatalogItems(data.items ?? []);
    } catch {
      setCatalogItems([]);
    } finally {
      setCatalogLoading(false);
    }
  }, []);

  // Fetch on catalog open + debounced search
  useEffect(() => {
    if (view !== "catalog") return;

    if (debounceRef.current) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(
      () => fetchCatalog(catalogSearch),
      catalogSearch ? 300 : 0,
    );

    return () => {
      if (debounceRef.current) clearTimeout(debounceRef.current);
    };
  }, [view, catalogSearch, fetchCatalog]);

  // Handle selecting a catalog entry
  const handleSelectEntry = (
    _dataset: CatalogDatasetItem,
    entry: CatalogDatasetEntry,
  ) => {
    const type: PanelType =
      entry.datasetType === "single_molecule" ? "sm" : "cell";

    if (entry.s3BaseUrl) {
      setRightPanelS3(entry.s3BaseUrl, type);
    } else if (entry.datasetId) {
      setRightPanel(entry.datasetId, type);
    }
  };

  const { id: currentId, s3Url: currentS3Url } = getCurrentDatasetInfo();
  const hasCurrentDataset = !!(currentId || currentS3Url);

  // Are we currently searching?
  const isSearching = catalogSearch.trim().length > 0;

  // Catalog browser view
  if (view === "catalog") {
    return (
      <div className="absolute inset-0 flex flex-col bg-black">
        {/* Header — mt-16 clears navbar */}
        <div className="flex items-center gap-2 px-4 pt-4 pb-2 mt-16">
          <Button
            isIconOnly
            className="text-white/70"
            size="sm"
            variant="light"
            onPress={() => {
              setView("main");
              setCatalogSearch("");
            }}
          >
            <svg
              className="w-5 h-5"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M15.75 19.5L8.25 12l7.5-7.5"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          </Button>
          <h3 className="text-sm font-semibold text-white">Browse Catalog</h3>
        </div>

        {/* Search */}
        <div className="px-4 pb-3">
          <Input
            classNames={{
              input: "text-white",
              inputWrapper: "bg-white/5 border-white/20 hover:bg-white/10",
            }}
            placeholder="Search datasets..."
            size="sm"
            value={catalogSearch}
            onChange={(e) => setCatalogSearch(e.target.value)}
          />
        </div>

        {/* Dataset list */}
        <div className="flex-1 overflow-y-auto px-4 pb-4">
          {catalogLoading && catalogItems.length === 0 ? (
            <div className="flex flex-col gap-4">
              <SkeletonGrid count={4} />
            </div>
          ) : (
            <div className="flex flex-col gap-4">
              {/* Featured section — shown when not searching */}
              {!isSearching && featuredItems.length > 0 && (
                <div>
                  <p className="text-xs font-medium text-white/40 uppercase tracking-wider mb-2">
                    Featured
                  </p>
                  <div className="grid grid-cols-2 gap-2">
                    {featuredItems.map((item) => (
                      <ExploreDatasetCard
                        key={item.id}
                        dataset={item}
                        onSelect={handleSelectEntry}
                      />
                    ))}
                  </div>
                </div>
              )}

              {/* All / search results */}
              <div>
                {!isSearching && catalogItems.length > 0 && (
                  <p className="text-xs font-medium text-white/40 uppercase tracking-wider mb-2">
                    All Datasets
                  </p>
                )}
                {catalogLoading ? (
                  <SkeletonGrid count={4} />
                ) : catalogItems.length === 0 ? (
                  <p className="text-sm text-white/40 text-center py-8">
                    No datasets found
                  </p>
                ) : (
                  <div className="grid grid-cols-2 gap-2">
                    {catalogItems.map((item) => (
                      <ExploreDatasetCard
                        key={item.id}
                        dataset={item}
                        onSelect={handleSelectEntry}
                      />
                    ))}
                  </div>
                )}
              </div>
            </div>
          )}
        </div>
      </div>
    );
  }

  // Link input view
  if (view === "link") {
    return (
      <div className="absolute inset-0 flex items-center justify-center bg-black">
        <div className="flex flex-col items-center gap-6 max-w-sm w-full px-8">
          <div className="text-center mb-2">
            <h3 className="text-lg font-semibold text-white">Paste Link</h3>
            <p className="text-sm text-white/50 mt-1">
              Enter a dataset URL or ID
            </p>
          </div>
          <div className="flex flex-col gap-2 w-full">
            <Input
              classNames={{
                input: "text-white",
                inputWrapper: "bg-white/5 border-white/20 hover:bg-white/10",
              }}
              placeholder="Paste dataset URL or ID..."
              size="sm"
              value={linkInput}
              onChange={(e) => setLinkInput(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") handlePasteLink();
              }}
            />
            <div className="flex gap-2">
              <Button
                className="flex-1"
                color="primary"
                isDisabled={!linkInput.trim()}
                size="sm"
                onPress={handlePasteLink}
              >
                Load
              </Button>
              <Button
                className="flex-1"
                size="sm"
                variant="bordered"
                onPress={() => {
                  setView("main");
                  setLinkInput("");
                }}
              >
                Cancel
              </Button>
            </div>
          </div>
        </div>
      </div>
    );
  }

  // Main picker view
  return (
    <div className="absolute inset-0 flex flex-col bg-black overflow-y-auto">
      <div className="flex flex-col items-center gap-6 w-full px-6 pt-24 pb-8">
        <div className="text-center mb-2">
          <h3 className="text-lg font-semibold text-white">Choose a Dataset</h3>
          <p className="text-sm text-white/50 mt-1">
            Select what to display in this panel
          </p>
        </div>

        <div className="flex flex-col gap-3 w-full max-w-sm">
          {/* Same Dataset */}
          {hasCurrentDataset && (
            <Button
              className="w-full justify-start"
              color="primary"
              variant="bordered"
              onPress={handleSameDataset}
            >
              <svg
                className="w-5 h-5 mr-2 flex-shrink-0"
                fill="none"
                stroke="currentColor"
                strokeWidth={1.5}
                viewBox="0 0 24 24"
              >
                <path
                  d="M15.75 17.25v3.375c0 .621-.504 1.125-1.125 1.125h-9.75a1.125 1.125 0 01-1.125-1.125V7.875c0-.621.504-1.125 1.125-1.125H6.75a9.06 9.06 0 011.5.124m7.5 10.376h3.375c.621 0 1.125-.504 1.125-1.125V11.25c0-4.46-3.243-8.161-7.5-8.876a9.06 9.06 0 00-1.5-.124H9.375c-.621 0-1.125.504-1.125 1.125v3.5m7.5 10.375H9.375a1.125 1.125 0 01-1.125-1.125v-9.25m12 6.625v-1.875a3.375 3.375 0 00-3.375-3.375h-1.5a1.125 1.125 0 01-1.125-1.125v-1.5a3.375 3.375 0 00-3.375-3.375H9.75"
                  strokeLinecap="round"
                  strokeLinejoin="round"
                />
              </svg>
              Same Dataset
            </Button>
          )}

          {/* Paste Link */}
          <Button
            className="w-full justify-start"
            variant="bordered"
            onPress={() => setView("link")}
          >
            <svg
              className="w-5 h-5 mr-2 flex-shrink-0"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M13.19 8.688a4.5 4.5 0 011.242 7.244l-4.5 4.5a4.5 4.5 0 01-6.364-6.364l1.757-1.757m9.86-3.06a4.5 4.5 0 00-1.242-7.244l-4.5-4.5a4.5 4.5 0 00-6.364 6.364l1.757 1.757"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
            Paste Link
          </Button>

          {/* Browse Catalog */}
          <Button
            className="w-full justify-start"
            variant="bordered"
            onPress={() => setView("catalog")}
          >
            <svg
              className="w-5 h-5 mr-2 flex-shrink-0"
              fill="none"
              stroke="currentColor"
              strokeWidth={1.5}
              viewBox="0 0 24 24"
            >
              <path
                d="M12 21v-8.25M15.75 21v-8.25M8.25 21v-8.25M3 9l9-6 9 6m-1.5 12V10.332A48.36 48.36 0 0012 9.75c-2.551 0-5.056.2-7.5.582V21M3 21h18M12 6.75h.008v.008H12V6.75z"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
            Browse Catalog
          </Button>
        </div>

        {/* Featured datasets shown on main view */}
        <div className="w-full mt-4">
          <p className="text-xs font-medium text-white/40 uppercase tracking-wider mb-3">
            Featured Datasets
          </p>
          {!featuredLoaded ? (
            <SkeletonGrid count={4} />
          ) : featuredItems.length > 0 ? (
            <div className="grid grid-cols-2 gap-3">
              {featuredItems.map((item) => (
                <ExploreDatasetCard
                  key={item.id}
                  dataset={item}
                  onSelect={handleSelectEntry}
                />
              ))}
            </div>
          ) : null}
        </div>
      </div>
    </div>
  );
}
