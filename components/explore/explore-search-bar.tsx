"use client";

import type { ExploreFilters } from "./types";

import { Input } from "@heroui/input";
import { Chip } from "@heroui/chip";

import { FilterPill } from "./filter-pill";

import { SearchIcon } from "@/components/icons";

interface ExploreSearchBarProps {
  search: string;
  onSearchChange: (value: string) => void;
  species: string;
  onSpeciesChange: (value: string) => void;
  tissue: string;
  onTissueChange: (value: string) => void;
  platform: string;
  onPlatformChange: (value: string) => void;
  geneSearch: string;
  onGeneSearchChange: (gene: string) => void;
  geneChips: string[];
  onGeneChipsChange: (genes: string[]) => void;
  filters: ExploreFilters;
  /** Category filter (was tabs): "" = All, else featured/bil/community/internal. */
  category: string;
  onCategoryChange: (value: string) => void;
  /** Available category values (excludes "All"); Internal only for admins. */
  categoryOptions: string[];
}

const CATEGORY_LABELS: Record<string, string> = {
  featured: "Featured",
  bil: "BIL",
  community: "Community",
  internal: "Internal",
};

export function ExploreSearchBar({
  search,
  onSearchChange,
  species,
  onSpeciesChange,
  tissue,
  onTissueChange,
  platform,
  onPlatformChange,
  geneSearch,
  onGeneSearchChange,
  geneChips,
  onGeneChipsChange,
  filters,
  category,
  onCategoryChange,
  categoryOptions,
}: ExploreSearchBarProps) {
  return (
    <div className="flex flex-col gap-3 mb-6">
      <div className="flex flex-col sm:flex-row gap-3">
        <Input
          className="flex-1"
          classNames={{ inputWrapper: "bg-default-100" }}
          placeholder="Search datasets..."
          startContent={
            <SearchIcon className="text-default-400 pointer-events-none flex-shrink-0" />
          }
          value={search}
          onValueChange={onSearchChange}
        />

        <Input
          className="flex-1 sm:max-w-xs"
          classNames={{ inputWrapper: "bg-default-100" }}
          placeholder="Search by gene (Enter for exact)"
          value={geneSearch}
          onKeyDown={(e) => {
            if (e.key === "Enter") {
              e.preventDefault();
              const gene = geneSearch.trim();

              if (
                gene &&
                !geneChips.some((g) => g.toLowerCase() === gene.toLowerCase())
              ) {
                onGeneChipsChange([...geneChips, gene]);
              }
              onGeneSearchChange("");
            }
          }}
          onValueChange={onGeneSearchChange}
        />
      </div>

      {/* LinkedIn-style filter pill bar — wraps when it runs out of room. */}
      <div className="flex flex-wrap gap-2">
        <FilterPill
          hideSearch
          label="Category"
          options={categoryOptions}
          renderOption={(v) => CATEGORY_LABELS[v] ?? v}
          value={category}
          onChange={onCategoryChange}
        />
        {filters.species.length > 0 && (
          <FilterPill
            label="Species"
            options={filters.species}
            value={species}
            onChange={onSpeciesChange}
          />
        )}
        {filters.tissues.length > 0 && (
          <FilterPill
            label="Tissue"
            options={filters.tissues}
            value={tissue}
            onChange={onTissueChange}
          />
        )}
        {filters.platforms.length > 0 && (
          <FilterPill
            label="Platform"
            options={filters.platforms}
            value={platform}
            onChange={onPlatformChange}
          />
        )}
      </div>

      {/* Gene exact-match chips */}
      {geneChips.length > 0 && (
        <div className="flex flex-wrap gap-1.5">
          {geneChips.map((gene) => (
            <Chip
              key={gene}
              color="secondary"
              size="sm"
              variant="flat"
              onClose={() =>
                onGeneChipsChange(geneChips.filter((g) => g !== gene))
              }
            >
              {gene}
            </Chip>
          ))}
        </div>
      )}
    </div>
  );
}
