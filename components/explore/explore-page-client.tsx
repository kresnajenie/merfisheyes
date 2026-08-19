"use client";

import type {
  CatalogDatasetItem,
  ExploreFilters,
  ExploreApiResponse,
} from "./types";

import { useState, useEffect, useCallback, useRef } from "react";
import { useRouter, useSearchParams } from "next/navigation";
import { useSession } from "next-auth/react";

import { FeaturedDatasets } from "./featured-datasets";
import { ExploreSearchBar } from "./explore-search-bar";
import { ExploreDatasetGrid } from "./explore-dataset-grid";

// One page of the explore grid. Exported so other entry points that seed
// ExplorePageClient (e.g. the in-viewer explore modal) fetch exactly one
// page — a mismatched limit renders more cards than the counter claims.
export const PAGE_LIMIT = 20;

// Category filter values (empty string = "All"). Replaces the old tab bar.
const CATEGORY_OPTIONS = ["featured", "bil", "community"];

interface ExplorePageClientProps {
  initialItems: CatalogDatasetItem[];
  initialTotal: number;
  initialFeatured: CatalogDatasetItem[];
  initialFilters: ExploreFilters;
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
  /**
   * Set by the server page when its initial DB load failed (e.g. a transient
   * DB blip during ISR regeneration), so the initial* props are empty. When
   * true, the client hydrates the full public payload from /api/explore
   * instead of trusting the empty SSR snapshot.
   */
  ssrFailed?: boolean;
}

export function ExplorePageClient({
  initialItems,
  initialTotal,
  initialFeatured,
  initialFilters,
  onCardClick,
  disableUrlSync = false,
  ssrFailed = false,
}: ExplorePageClientProps) {
  const router = useRouter();
  const searchParams = useSearchParams();

  // Admin status is read client-side (not passed as a prop) so the page itself
  // can be statically cached — calling auth() during render would force
  // per-request rendering. Gates the Internal category; the API enforces admin
  // access independently.
  const { data: session } = useSession();
  const isAdmin =
    session?.user?.role === "ADMIN" || session?.user?.role === "SUPER_ADMIN";
  const categoryOptions = isAdmin
    ? [...CATEGORY_OPTIONS, "internal"]
    : CATEGORY_OPTIONS;

  // Initialize state from URL params (skip when URL sync is disabled —
  // otherwise the viewer's own params would bleed into the modal).
  const sp = disableUrlSync ? new URLSearchParams() : searchParams;
  // "" = All; otherwise featured/bil/community/internal.
  const [category, setCategory] = useState((sp.get("tab") || "").toLowerCase());
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
  // Debounced mirror of `geneSearch` — same contract as `debouncedSearch`:
  // the input stays live, but URL sync + fetches wait out the typing burst.
  const [debouncedGeneSearch, setDebouncedGeneSearch] = useState(
    sp.get("gene") || "",
  );

  useEffect(() => {
    const t = setTimeout(() => setDebouncedGeneSearch(geneSearch), 250);

    return () => clearTimeout(t);
  }, [geneSearch]);
  const [geneChips, setGeneChips] = useState<string[]>(
    sp.get("geneExact")?.split(",").filter(Boolean) || [],
  );
  const initialPage = Number(sp.get("page")) || 1;
  const [page, setPage] = useState(initialPage);

  // Whether the very first render is the clean default view. The SSR seed is
  // always page 1 of "All", so when the URL asks for a different page/category/
  // filter we must NOT paint that seed (it would flash page 1 before the real
  // data loads) — start empty + loading and let the mount fetch fill in.
  const initialIsDefault =
    !sp.get("tab") &&
    initialPage === 1 &&
    !(
      sp.get("q") ||
      sp.get("species") ||
      sp.get("tissue") ||
      sp.get("platform") ||
      sp.get("gene") ||
      sp.get("geneExact")
    );

  const [items, setItems] = useState(initialIsDefault ? initialItems : []);
  const [total, setTotal] = useState(initialIsDefault ? initialTotal : 0);
  // Featured carousel (shown only on the default "All" view).
  const [featuredItems, setFeaturedItems] = useState(initialFeatured);

  const [filters, setFilters] = useState(initialFilters);
  const [loading, setLoading] = useState(!initialIsDefault);

  const isInitialMount = useRef(true);
  // Snapshot of the filter/category values, so the page-reset effect fires only
  // on a real change — never on mount or a React Strict Mode remount (which
  // would clobber a deep-linked ?page=N back to 1).
  const prevFilterKey = useRef<string | null>(null);

  const hasActiveFilters = Boolean(
    debouncedSearch ||
      species ||
      tissue ||
      platform ||
      debouncedGeneSearch ||
      geneChips.length > 0,
  );

  // The clean default view (All, no filters, first page) is served straight
  // from the SSR seed — no fetch needed.
  const isDefaultView = category === "" && !hasActiveFilters && page === 1;

  // Sync state to URL (replace, not push, to avoid polluting history)
  useEffect(() => {
    if (disableUrlSync) return;
    // Skip on initial mount to avoid overwriting the URL we just read from
    if (isInitialMount.current) {
      isInitialMount.current = false;

      return;
    }

    const params = new URLSearchParams();

    if (category) params.set("tab", category);
    if (debouncedSearch) params.set("q", debouncedSearch);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (debouncedGeneSearch) params.set("gene", debouncedGeneSearch);
    if (geneChips.length > 0) params.set("geneExact", geneChips.join(","));
    if (page > 1) params.set("page", String(page));

    const qs = params.toString();

    router.replace(`/explore${qs ? `?${qs}` : ""}`, { scroll: false });
  }, [
    category,
    debouncedSearch,
    species,
    tissue,
    platform,
    debouncedGeneSearch,
    geneChips,
    page,
    router,
    disableUrlSync,
  ]);

  // In-flight list request — aborted when a newer one starts, so a slow stale
  // response can never overwrite fresher results (type "nt" → "ntrk" fast and
  // the "nt" query may finish last).
  const fetchAbortRef = useRef<AbortController | null>(null);

  const fetchData = useCallback(async () => {
    fetchAbortRef.current?.abort();
    const controller = new AbortController();

    fetchAbortRef.current = controller;

    setLoading(true);
    const params = new URLSearchParams();

    if (debouncedSearch) params.set("search", debouncedSearch);
    if (species) params.set("species", species);
    if (tissue) params.set("tissue", tissue);
    if (platform) params.set("platform", platform);
    if (debouncedGeneSearch) params.set("genes", debouncedGeneSearch);
    if (geneChips.length > 0) params.set("genesExact", geneChips.join(","));
    params.set("page", String(page));
    params.set("limit", String(PAGE_LIMIT));
    if (category) params.set("tab", category);

    try {
      const res = await fetch(`/api/explore?${params}`, {
        signal: controller.signal,
      });
      const data: ExploreApiResponse = await res.json();

      setItems(data.items);
      setTotal(data.total);
      setLoading(false);
    } catch (err) {
      // Superseded by a newer request — leave state to that request.
      if ((err as Error).name !== "AbortError") setLoading(false);
    }
  }, [
    debouncedSearch,
    species,
    tissue,
    platform,
    debouncedGeneSearch,
    geneChips,
    page,
    category,
  ]);

  // Recovery path when SSR failed: fetch the full public payload (All grid +
  // featured + filters) in one request so an empty shell fills in.
  const hydrateAll = useCallback(async () => {
    setLoading(true);
    try {
      const res = await fetch(`/api/explore?limit=${PAGE_LIMIT}&include=meta`);
      const data: ExploreApiResponse = await res.json();

      setItems(data.items);
      setTotal(data.total);
      setFeaturedItems(data.featured ?? []);
      if (data.filters) setFilters(data.filters);
    } finally {
      setLoading(false);
    }
  }, []);

  // Load data when filters / category / page change. The clean default view is
  // restored instantly from the SSR seed instead of refetching.
  useEffect(() => {
    if (isDefaultView) {
      if (ssrFailed) {
        hydrateAll();
      } else {
        setItems(initialItems);
        setTotal(initialTotal);
        setLoading(false);
      }

      return;
    }
    fetchData();
  }, [fetchData, isDefaultView, ssrFailed]);

  // Reset page to 1 when filters or category change — but only on a *real*
  // change. Comparing against a snapshot (not a run-count) means mount and
  // Strict Mode remounts don't reset a deep-linked ?page=N.
  useEffect(() => {
    const key = JSON.stringify([
      debouncedSearch,
      species,
      tissue,
      platform,
      debouncedGeneSearch,
      geneChips,
      category,
    ]);

    if (prevFilterKey.current === null) {
      prevFilterKey.current = key;

      return;
    }
    if (prevFilterKey.current !== key) {
      prevFilterKey.current = key;
      setPage(1);
    }
  }, [
    debouncedSearch,
    species,
    tissue,
    platform,
    debouncedGeneSearch,
    geneChips,
    category,
  ]);

  return (
    <>
      {/* Featured carousel — persistently on top of the "All" view (any
          filters/page), hidden once a specific category is selected. */}
      {category === "" && featuredItems.length > 0 && (
        <div className="mb-8">
          <FeaturedDatasets
            datasets={featuredItems}
            onCardClick={onCardClick}
            onViewAll={() => setCategory("featured")}
          />
        </div>
      )}

      <ExploreSearchBar
        category={category}
        categoryOptions={categoryOptions}
        filters={filters}
        geneChips={geneChips}
        geneSearch={geneSearch}
        platform={platform}
        search={search}
        species={species}
        tissue={tissue}
        onCategoryChange={setCategory}
        onGeneChipsChange={setGeneChips}
        onGeneSearchChange={setGeneSearch}
        onPlatformChange={setPlatform}
        onSearchChange={setSearch}
        onSpeciesChange={setSpecies}
        onTissueChange={setTissue}
      />

      <ExploreDatasetGrid
        datasets={items}
        geneHighlight={geneSearch}
        limit={PAGE_LIMIT}
        loading={loading}
        page={page}
        total={total}
        onCardClick={onCardClick}
        onPageChange={setPage}
      />
    </>
  );
}
