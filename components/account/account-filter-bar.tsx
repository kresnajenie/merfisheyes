"use client";

import { Input } from "@heroui/input";

import { FilterPill } from "@/components/explore/filter-pill";

export interface AccountFacets {
  species: string[];
  tissues: string[];
  platforms: string[];
}

const TYPE_LABELS: Record<string, string> = {
  single_cell: "Single cell",
  single_molecule: "Single molecule",
};

const STATUS_OPTIONS = ["COMPLETE", "PROCESSING", "QUEUED", "FAILED"];
const STATUS_LABELS: Record<string, string> = {
  COMPLETE: "Complete",
  PROCESSING: "Processing",
  QUEUED: "Queued",
  FAILED: "Failed",
};

interface AccountFilterBarProps {
  search: string;
  onSearchChange: (v: string) => void;
  type: string;
  onTypeChange: (v: string) => void;
  status: string;
  onStatusChange: (v: string) => void;
  species: string;
  onSpeciesChange: (v: string) => void;
  tissue: string;
  onTissueChange: (v: string) => void;
  platform: string;
  onPlatformChange: (v: string) => void;
  facets: AccountFacets;
}

/**
 * Explore-style filter bar for the account datasets grid: the same search box
 * and LinkedIn-style FilterPills as /explore, plus owner-only Type and Status
 * pills.
 */
export function AccountFilterBar({
  search,
  onSearchChange,
  type,
  onTypeChange,
  status,
  onStatusChange,
  species,
  onSpeciesChange,
  tissue,
  onTissueChange,
  platform,
  onPlatformChange,
  facets,
}: AccountFilterBarProps) {
  return (
    <div className="mb-6 flex flex-col gap-3">
      <Input
        classNames={{ inputWrapper: "bg-default-100" }}
        placeholder="Search your datasets…"
        value={search}
        onValueChange={onSearchChange}
      />

      <div className="flex flex-wrap gap-2">
        <FilterPill
          hideSearch
          label="Type"
          options={["single_cell", "single_molecule"]}
          renderOption={(v) => TYPE_LABELS[v] ?? v}
          value={type}
          onChange={onTypeChange}
        />
        <FilterPill
          hideSearch
          label="Status"
          options={STATUS_OPTIONS}
          renderOption={(v) => STATUS_LABELS[v] ?? v}
          value={status}
          onChange={onStatusChange}
        />
        {facets.species.length > 0 && (
          <FilterPill
            label="Species"
            options={facets.species}
            value={species}
            onChange={onSpeciesChange}
          />
        )}
        {facets.tissues.length > 0 && (
          <FilterPill
            label="Tissue"
            options={facets.tissues}
            value={tissue}
            onChange={onTissueChange}
          />
        )}
        {facets.platforms.length > 0 && (
          <FilterPill
            label="Platform"
            options={facets.platforms}
            value={platform}
            onChange={onPlatformChange}
          />
        )}
      </div>
    </div>
  );
}
