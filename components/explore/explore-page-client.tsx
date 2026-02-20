"use client";

import { useState, useEffect, useCallback } from "react";

import { FeaturedDatasets } from "./featured-datasets";
import { BilDatasets } from "./bil-datasets";
import { ExploreSearchBar } from "./explore-search-bar";
import { ExploreFilterChips } from "./explore-filter-chips";
import { ExploreDatasetGrid } from "./explore-dataset-grid";
import type { CatalogDatasetItem, ExploreFilters, ExploreApiResponse } from "./types";

interface ExplorePageClientProps {
  initialItems: CatalogDatasetItem[];
  initialTotal: number;
  initialFeatured: CatalogDatasetItem[];
  initialBil: CatalogDatasetItem[];
  initialFilters: ExploreFilters;
}

export function ExplorePageClient({
  initialItems,
  initialTotal,
  initialFeatured,
  initialBil,
  initialFilters,
}: ExplorePageClientProps) {
  const [items, setItems] = useState(initialItems);
  const [total, setTotal] = useState(initialTotal);
  const [featured] = useState(initialFeatured);
  const [bil] = useState(initialBil);
  const [filters, setFilters] = useState(initialFilters);
  const [loading, setLoading] = useState(false);

  const [search, setSearch] = useState("");
  const [species, setSpecies] = useState("");
  const [tissue, setTissue] = useState("");
  const [platform, setPlatform] = useState("");
  const [page, setPage] = useState(1);

  const hasActiveFilters = search || species || tissue || platform;

  const fetchData = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams();
    if (search) params.set("search", search);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    params.set("page", String(page));

    const res = await fetch(`/api/explore?${params}`);
    const data: ExploreApiResponse = await res.json();
    setItems(data.items);
    setTotal(data.total);
    setFilters(data.filters);
    setLoading(false);
  }, [search, species, tissue, platform, page]);

  // Refetch when filters change (but not on initial render with SSR data)
  useEffect(() => {
    if (hasActiveFilters || page > 1) {
      fetchData();
    }
  }, [fetchData, hasActiveFilters, page]);

  // Reset page when filters change
  useEffect(() => {
    setPage(1);
  }, [search, species, tissue, platform]);

  const activeFilters = [
    species && { key: "species", label: "Species", value: species, onClear: () => setSpecies("") },
    tissue && { key: "tissue", label: "Tissue", value: tissue, onClear: () => setTissue("") },
    platform && { key: "platform", label: "Platform", value: platform, onClear: () => setPlatform("") },
  ].filter(Boolean) as { key: string; label: string; value: string; onClear: () => void }[];

  return (
    <>
      <FeaturedDatasets datasets={featured} />

      <BilDatasets datasets={bil} />

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
        datasets={items}
        limit={20}
        loading={loading}
        page={page}
        total={total}
        onPageChange={setPage}
      />
    </>
  );
}
