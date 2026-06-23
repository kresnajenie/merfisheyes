"use client";

import type {
  CatalogDatasetItem,
  ExploreFilters,
  ExploreApiResponse,
} from "./types";

import { useState, useEffect, useCallback, useRef } from "react";
import { Tabs, Tab } from "@heroui/tabs";
import { useRouter, useSearchParams } from "next/navigation";

import { FeaturedDatasets } from "./featured-datasets";
import { ExploreSearchBar } from "./explore-search-bar";
import { ExploreDatasetGrid } from "./explore-dataset-grid";

const PAGE_LIMIT = 50;

type ExploreTab = "all" | "featured" | "bil" | "internal";

interface ExplorePageClientProps {
  initialItems: CatalogDatasetItem[];
  initialTotal: number;
  initialFeatured: CatalogDatasetItem[];
  initialBil: CatalogDatasetItem[];
  initialInternal: CatalogDatasetItem[];
  initialFilters: ExploreFilters;
  isAdmin: boolean;
  /**
   * When provided, intercepts dataset card clicks instead of letting the
   * card link/navigate. Used by the in-viewer explore modal to switch into
   * the in-modal detail view (multi-entry / BIL) or directly load
   * single-entry datasets into the right panel.
   */
  onCardClick?: (dataset: CatalogDatasetItem) => void;
  /**
   * When true, skip both reading initial filter state from `?q=...` etc.
   * AND syncing state back into the URL. Used by the modal so opening it
   * on `/viewer/X` doesn't try to push `?q=...` onto the viewer URL or
   * read leftover params from it.
   */
  disableUrlSync?: boolean;
}

export function ExplorePageClient({
  initialItems,
  initialTotal,
  initialFeatured,
  initialBil,
  initialInternal,
  initialFilters,
  isAdmin,
  onCardClick,
  disableUrlSync = false,
}: ExplorePageClientProps) {
  const router = useRouter();
  const searchParams = useSearchParams();

  // Initialize state from URL params (skip when URL sync is disabled —
  // otherwise the viewer's own params would bleed into the modal).
  const sp = disableUrlSync ? new URLSearchParams() : searchParams;
  const [activeTab, setActiveTab] = useState<ExploreTab>(
    (sp.get("tab") as ExploreTab) || "all",
  );
  const [search, setSearch] = useState(sp.get("q") || "");
  // Debounced mirror of `search`: the input stays bound to `search` for
  // instant typing feedback, but URL sync + fetch use the debounced value so
  // we don't spam the API / history on every keystroke.
  const [debouncedSearch, setDebouncedSearch] = useState(sp.get("q") || "");

  useEffect(() => {
    const t = setTimeout(() => setDebouncedSearch(search), 250);

    return () => clearTimeout(t);
  }, [search]);
  const [species, setSpecies] = useState(sp.get("species") || "");
  const [tissue, setTissue] = useState(sp.get("tissue") || "");
  const [platform, setPlatform] = useState(sp.get("platform") || "");
  const [geneSearch, setGeneSearch] = useState(sp.get("gene") || "");
  const [geneChips, setGeneChips] = useState<string[]>(
    sp.get("geneExact")?.split(",").filter(Boolean) || [],
  );
  const [page, setPage] = useState(Number(sp.get("page")) || 1);

  // Per-tab data
  const [allItems, setAllItems] = useState(initialItems);
  const [allTotal, setAllTotal] = useState(initialTotal);
  const [featuredItems, setFeaturedItems] = useState(initialFeatured);
  const [featuredTotal, setFeaturedTotal] = useState(initialFeatured.length);
  const [bilItems, setBilItems] = useState(initialBil);
  const [bilTotal, setBilTotal] = useState(initialBil.length);
  const [internalItems, setInternalItems] = useState(initialInternal);
  const [internalTotal, setInternalTotal] = useState(initialInternal.length);

  const [filters, setFilters] = useState(initialFilters);
  const [loading, setLoading] = useState(false);

  const [hasFetched, setHasFetched] = useState(false);
  const isInitialMount = useRef(true);

  const hasActiveFilters =
    debouncedSearch ||
    species ||
    tissue ||
    platform ||
    geneSearch ||
    geneChips.length > 0;

  // Sync state to URL (replace, not push, to avoid polluting history)
  useEffect(() => {
    if (disableUrlSync) return;
    // Skip on initial mount to avoid overwriting the URL we just read from
    if (isInitialMount.current) {
      isInitialMount.current = false;

      return;
    }

    const params = new URLSearchParams();

    if (activeTab !== "all") params.set("tab", activeTab);
    if (debouncedSearch) params.set("q", debouncedSearch);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (geneSearch) params.set("gene", geneSearch);
    if (geneChips.length > 0) params.set("geneExact", geneChips.join(","));
    if (page > 1) params.set("page", String(page));

    const qs = params.toString();

    router.replace(`/explore${qs ? `?${qs}` : ""}`, { scroll: false });
  }, [
    activeTab,
    debouncedSearch,
    species,
    tissue,
    platform,
    geneSearch,
    geneChips,
    page,
    router,
    disableUrlSync,
  ]);

  const fetchData = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams();

    if (debouncedSearch) params.set("search", debouncedSearch);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (geneSearch) params.set("genes", geneSearch);
    if (geneChips.length > 0) params.set("genesExact", geneChips.join(","));
    params.set("page", String(page));
    params.set("limit", String(PAGE_LIMIT));

    if (activeTab !== "all") {
      params.set("tab", activeTab);
    }

    const res = await fetch(`/api/explore?${params}`);
    const data: ExploreApiResponse = await res.json();

    switch (activeTab) {
      case "all":
        setAllItems(data.items);
        setAllTotal(data.total);
        break;
      case "featured":
        setFeaturedItems(data.items);
        setFeaturedTotal(data.total);
        break;
      case "bil":
        setBilItems(data.items);
        setBilTotal(data.total);
        break;
      case "internal":
        setInternalItems(data.items);
        setInternalTotal(data.total);
        break;
    }

    setFilters(data.filters);
    setLoading(false);
    setHasFetched(true);
  }, [
    debouncedSearch,
    species,
    tissue,
    platform,
    geneSearch,
    geneChips,
    page,
    activeTab,
  ]);

  // Refetch when filters/page/tab change
  useEffect(() => {
    if (hasActiveFilters || page > 1 || hasFetched) {
      fetchData();
    }
  }, [fetchData, hasActiveFilters, page]);

  // Fetch when switching to a non-all tab for the first time
  useEffect(() => {
    if (activeTab !== "all") {
      fetchData();
    } else if (!hasActiveFilters && page === 1) {
      setAllItems(initialItems);
      setAllTotal(initialTotal);
      setHasFetched(false);
    }
  }, [activeTab]);

  // Reset page when filters or tab change
  useEffect(() => {
    setPage(1);
  }, [
    debouncedSearch,
    species,
    tissue,
    platform,
    geneSearch,
    geneChips,
    activeTab,
  ]);

  const renderSearchAndGrid = (
    datasetItems: CatalogDatasetItem[],
    datasetTotal: number,
  ) => (
    <>
      <ExploreSearchBar
        filters={filters}
        geneChips={geneChips}
        geneSearch={geneSearch}
        platform={platform}
        search={search}
        species={species}
        tissue={tissue}
        onGeneChipsChange={setGeneChips}
        onGeneSearchChange={setGeneSearch}
        onPlatformChange={setPlatform}
        onSearchChange={setSearch}
        onSpeciesChange={setSpecies}
        onTissueChange={setTissue}
      />

      <ExploreDatasetGrid
        datasets={datasetItems}
        geneHighlight={geneSearch}
        limit={PAGE_LIMIT}
        loading={loading}
        page={page}
        total={datasetTotal}
        onCardClick={onCardClick}
        onPageChange={setPage}
      />
    </>
  );

  return (
    <>
      <Tabs
        aria-label="Explore tabs"
        className="mb-6"
        selectedKey={activeTab}
        variant="underlined"
        onSelectionChange={(key) => setActiveTab(key as ExploreTab)}
      >
        <Tab key="all" title="All">
          <div className="mb-8">
            <FeaturedDatasets
              datasets={featuredItems}
              onCardClick={onCardClick}
              onViewAll={() => setActiveTab("featured")}
            />
          </div>
          {renderSearchAndGrid(allItems, allTotal)}
        </Tab>
        <Tab key="featured" title="Featured">
          {renderSearchAndGrid(featuredItems, featuredTotal)}
        </Tab>
        <Tab key="bil" title="BIL">
          {renderSearchAndGrid(bilItems, bilTotal)}
        </Tab>
        {isAdmin && (
          <Tab key="internal" title="Internal">
            {renderSearchAndGrid(internalItems, internalTotal)}
          </Tab>
        )}
      </Tabs>
    </>
  );
}
