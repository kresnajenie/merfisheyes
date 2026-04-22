"use client";

import { useState, useEffect, useCallback, useRef } from "react";
import { Tabs, Tab } from "@heroui/tabs";
import { useRouter, useSearchParams } from "next/navigation";

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
  const router = useRouter();
  const searchParams = useSearchParams();

  // Initialize state from URL params
  const [activeTab, setActiveTab] = useState<ExploreTab>(
    (searchParams.get("tab") as ExploreTab) || "all",
  );
  const [search, setSearch] = useState(searchParams.get("q") || "");
  const [species, setSpecies] = useState(searchParams.get("species") || "");
  const [tissue, setTissue] = useState(searchParams.get("tissue") || "");
  const [platform, setPlatform] = useState(searchParams.get("platform") || "");
  const [geneSearch, setGeneSearch] = useState(searchParams.get("gene") || "");
  const [geneChips, setGeneChips] = useState<string[]>(
    searchParams.get("geneExact")?.split(",").filter(Boolean) || [],
  );
  const [page, setPage] = useState(Number(searchParams.get("page")) || 1);

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

  const hasActiveFilters = search || species || tissue || platform || geneSearch || geneChips.length > 0;

  // Sync state to URL (replace, not push, to avoid polluting history)
  useEffect(() => {
    // Skip on initial mount to avoid overwriting the URL we just read from
    if (isInitialMount.current) {
      isInitialMount.current = false;
      return;
    }

    const params = new URLSearchParams();
    if (activeTab !== "all") params.set("tab", activeTab);
    if (search) params.set("q", search);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (geneSearch) params.set("gene", geneSearch);
    if (geneChips.length > 0) params.set("geneExact", geneChips.join(","));
    if (page > 1) params.set("page", String(page));

    const qs = params.toString();
    router.replace(`/explore${qs ? `?${qs}` : ""}`, { scroll: false });
  }, [activeTab, search, species, tissue, platform, geneSearch, geneChips, page, router]);

  const fetchData = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams();
    if (search) params.set("search", search);
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
  }, [search, species, tissue, platform, geneSearch, geneChips, page, activeTab]);

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
  }, [activeTab]); // eslint-disable-line react-hooks/exhaustive-deps

  // Reset page when filters or tab change
  useEffect(() => {
    setPage(1);
  }, [search, species, tissue, platform, geneSearch, geneChips, activeTab]);

  const activeFilters = [
    species && { key: "species", label: "Species", value: species, onClear: () => setSpecies("") },
    tissue && { key: "tissue", label: "Tissue", value: tissue, onClear: () => setTissue("") },
    platform && { key: "platform", label: "Platform", value: platform, onClear: () => setPlatform("") },
    geneSearch && { key: "gene", label: "Gene (includes)", value: geneSearch, onClear: () => setGeneSearch("") },
    geneChips.length > 0 && { key: "geneExact", label: "Gene (exact)", value: geneChips.join(", "), onClear: () => setGeneChips([]) },
  ].filter(Boolean) as { key: string; label: string; value: string; onClear: () => void }[];

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

      <ExploreFilterChips filters={activeFilters} />

      <ExploreDatasetGrid
        datasets={datasetItems}
        geneHighlight={geneSearch}
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
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-8">
            <FeaturedDatasets datasets={featuredItems} onViewAll={() => setActiveTab("featured")} />
            <BilDatasets datasets={bilItems} onViewAll={() => setActiveTab("bil")} />
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
