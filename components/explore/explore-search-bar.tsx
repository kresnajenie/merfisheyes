"use client";

import { Input } from "@heroui/input";
import { Select, SelectItem } from "@heroui/select";
import { Chip } from "@heroui/chip";
import { SearchIcon } from "@/components/icons";

import type { ExploreFilters } from "./types";

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
}

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
          onValueChange={onGeneSearchChange}
          onKeyDown={(e) => {
            if (e.key === "Enter") {
              e.preventDefault();
              const gene = geneSearch.trim();
              if (gene && !geneChips.some((g) => g.toLowerCase() === gene.toLowerCase())) {
                onGeneChipsChange([...geneChips, gene]);
              }
              onGeneSearchChange("");
            }
          }}
        />

        {filters.species.length > 0 && (
          <Select
            className="w-full sm:w-40"
            label="Species"
            selectedKeys={species ? [species] : []}
            size="sm"
            onSelectionChange={(keys) => {
              const val = Array.from(keys)[0] as string | undefined;
              onSpeciesChange(val ?? "");
            }}
          >
            {filters.species.map((s) => (
              <SelectItem key={s}>{s}</SelectItem>
            ))}
          </Select>
        )}

        {filters.tissues.length > 0 && (
          <Select
            className="w-full sm:w-40"
            label="Tissue"
            selectedKeys={tissue ? [tissue] : []}
            size="sm"
            onSelectionChange={(keys) => {
              const val = Array.from(keys)[0] as string | undefined;
              onTissueChange(val ?? "");
            }}
          >
            {filters.tissues.map((t) => (
              <SelectItem key={t}>{t}</SelectItem>
            ))}
          </Select>
        )}

        {filters.platforms.length > 0 && (
          <Select
            className="w-full sm:w-40"
            label="Platform"
            selectedKeys={platform ? [platform] : []}
            size="sm"
            onSelectionChange={(keys) => {
              const val = Array.from(keys)[0] as string | undefined;
              onPlatformChange(val ?? "");
            }}
          >
            {filters.platforms.map((p) => (
              <SelectItem key={p}>{p}</SelectItem>
            ))}
          </Select>
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
