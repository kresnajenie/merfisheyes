"use client";

import { Input } from "@heroui/input";
import { Select, SelectItem } from "@heroui/select";
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
  filters,
}: ExploreSearchBarProps) {
  return (
    <div className="flex flex-col sm:flex-row gap-3 mb-6">
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
  );
}
