"use client";

import { useState, useEffect, useCallback } from "react";
import { Tabs, Tab } from "@heroui/tabs";

import { FeaturedDatasets } from "./featured-datasets";
import { BilDatasets } from "./bil-datasets";
import { ExploreSearchBar } from "./explore-search-bar";
import { ExploreFilterChips } from "./explore-filter-chips";
import { ExploreDatasetGrid } from "./explore-dataset-grid";
import type { CatalogDatasetItem, ExploreFilters, ExploreApiResponse } from "./types";

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
}

export function ExplorePageClient({
  initialItems,
  initialTotal,
  initialFeatured,
  initialBil,
  initialInternal,
  initialFilters,
  isAdmin,
}: ExplorePageClientProps) {
  const [activeTab, setActiveTab] = useState<ExploreTab>("all");

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

  const [search, setSearch] = useState("");
  const [species, setSpecies] = useState("");
  const [tissue, setTissue] = useState("");
  const [platform, setPlatform] = useState("");
  const [page, setPage] = useState(1);

  // Track whether a tab has been fetched with current filters (to avoid refetching SSR data)
  const [hasFetched, setHasFetched] = useState(false);

  const hasActiveFilters = search || species || tissue || platform;

  const fetchData = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams();
    if (search) params.set("search", search);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
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
  }, [search, species, tissue, platform, page, activeTab]);

  // Refetch when filters/page/tab change
  useEffect(() => {
    if (hasActiveFilters || page > 1 || hasFetched) {
      fetchData();
    }
  }, [fetchData, hasActiveFilters, page]);

  // Fetch when switching to a non-all tab for the first time (SSR only pre-fetches "all")
  useEffect(() => {
    if (activeTab !== "all") {
      fetchData();
    } else if (!hasActiveFilters && page === 1) {
      // Reset to SSR data when switching back to "all" with no filters
      setAllItems(initialItems);
      setAllTotal(initialTotal);
      setHasFetched(false);
    }
  }, [activeTab]); // eslint-disable-line react-hooks/exhaustive-deps

  // Reset page when filters or tab change
  useEffect(() => {
    setPage(1);
  }, [search, species, tissue, platform, activeTab]);

  const activeFilters = [
    species && { key: "species", label: "Species", value: species, onClear: () => setSpecies("") },
    tissue && { key: "tissue", label: "Tissue", value: tissue, onClear: () => setTissue("") },
    platform && { key: "platform", label: "Platform", value: platform, onClear: () => setPlatform("") },
  ].filter(Boolean) as { key: string; label: string; value: string; onClear: () => void }[];

  const renderSearchAndGrid = (
    datasetItems: CatalogDatasetItem[],
    datasetTotal: number,
  ) => (
    <>
      <ExploreSearchBar
        filters={filters}
        platform={platform}
        search={search}
        species={species}
        tissue={tissue}
        onPlatformChange={setPlatform}
        onSearchChange={setSearch}
        onSpeciesChange={setSpecies}
        onTissueChange={setTissue}
      />

      <ExploreFilterChips filters={activeFilters} />

      <ExploreDatasetGrid
        datasets={datasetItems}
        limit={PAGE_LIMIT}
        loading={loading}
        page={page}
        total={datasetTotal}
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
          <FeaturedDatasets datasets={featuredItems} />
          <BilDatasets datasets={bilItems} />
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
