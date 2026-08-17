"use client";

import type { DatasetRow } from "./types";

import { useMemo, useState } from "react";
import { Input } from "@heroui/input";
import { Select, SelectItem } from "@heroui/select";
import { Button } from "@heroui/button";
import { Chip } from "@heroui/chip";

export type SortKey = "date" | "name" | "cells" | "genes" | "views";

export interface DatasetFilters {
  search: string;
  type: "all" | "single_cell" | "single_molecule";
  species: string;
  tissue: string;
  platform: string;
  status: string;
  tag: string;
  minGenes: string;
  from: string;
  to: string;
  sort: SortKey;
  dir: "asc" | "desc";
}

export const DEFAULT_FILTERS: DatasetFilters = {
  search: "",
  type: "all",
  species: "",
  tissue: "",
  platform: "",
  status: "all",
  tag: "",
  minGenes: "",
  from: "",
  to: "",
  sort: "date",
  dir: "desc",
};

const SORT_OPTIONS: { key: SortKey; label: string }[] = [
  { key: "date", label: "Date added" },
  { key: "name", label: "Name" },
  { key: "cells", label: "Cells / molecules" },
  { key: "genes", label: "Genes" },
  { key: "views", label: "Views" },
];

/** Distinct, sorted, non-empty values of a field across the dataset list. */
function distinct(rows: DatasetRow[], pick: (d: DatasetRow) => string | null) {
  return Array.from(
    new Set(rows.map(pick).filter((v): v is string => !!v)),
  ).sort((a, b) => a.localeCompare(b));
}

/** Count of filters that are narrowing the list (for the "Filters" badge). */
export function activeFilterCount(f: DatasetFilters): number {
  let n = 0;

  if (f.type !== "all") n++;
  if (f.species) n++;
  if (f.tissue) n++;
  if (f.platform) n++;
  if (f.status !== "all") n++;
  if (f.tag) n++;
  if (f.minGenes) n++;
  if (f.from) n++;
  if (f.to) n++;

  return n;
}

interface AccountFilterBarProps {
  datasets: DatasetRow[];
  filters: DatasetFilters;
  onChange: (next: DatasetFilters) => void;
  /** Count after filtering, for the "Showing X of Y" line. */
  shown: number;
}

export function AccountFilterBar({
  datasets,
  filters,
  onChange,
  shown,
}: AccountFilterBarProps) {
  const [expanded, setExpanded] = useState(false);

  const speciesOpts = useMemo(
    () => distinct(datasets, (d) => d.species),
    [datasets],
  );
  const tissueOpts = useMemo(
    () => distinct(datasets, (d) => d.tissue),
    [datasets],
  );
  const platformOpts = useMemo(
    () => distinct(datasets, (d) => d.platform),
    [datasets],
  );
  const statusOpts = useMemo(
    () => distinct(datasets, (d) => d.status),
    [datasets],
  );
  const tagOpts = useMemo(
    () =>
      Array.from(new Set(datasets.flatMap((d) => d.tags ?? []))).sort((a, b) =>
        a.localeCompare(b),
      ),
    [datasets],
  );

  const set = <K extends keyof DatasetFilters>(
    key: K,
    value: DatasetFilters[K],
  ) => onChange({ ...filters, [key]: value });

  const activeCount = activeFilterCount(filters);

  return (
    <div className="flex flex-col gap-3 mb-4">
      {/* Primary row: search + type + sort */}
      <div className="flex flex-wrap items-center gap-2">
        <Input
          isClearable
          className="min-w-[180px] flex-1"
          placeholder="Search name, species, tissue, tags…"
          size="sm"
          startContent={<span className="text-default-400 text-sm">🔍</span>}
          value={filters.search}
          onClear={() => set("search", "")}
          onValueChange={(v) => set("search", v)}
        />

        <Select
          aria-label="Type"
          className="w-36"
          selectedKeys={[filters.type]}
          size="sm"
          onChange={(e) =>
            set("type", (e.target.value || "all") as DatasetFilters["type"])
          }
        >
          <SelectItem key="all">All types</SelectItem>
          <SelectItem key="single_cell">Single cell</SelectItem>
          <SelectItem key="single_molecule">Single molecule</SelectItem>
        </Select>

        <Select
          aria-label="Sort by"
          className="w-44"
          selectedKeys={[filters.sort]}
          size="sm"
          startContent={<span className="text-default-400 text-xs">Sort</span>}
          onChange={(e) => set("sort", (e.target.value || "date") as SortKey)}
        >
          {SORT_OPTIONS.map((o) => (
            <SelectItem key={o.key}>{o.label}</SelectItem>
          ))}
        </Select>

        <Button
          isIconOnly
          aria-label="Toggle sort direction"
          size="sm"
          variant="flat"
          onPress={() => set("dir", filters.dir === "asc" ? "desc" : "asc")}
        >
          {filters.dir === "asc" ? "↑" : "↓"}
        </Button>

        <Button
          size="sm"
          variant={expanded || activeCount > 0 ? "solid" : "flat"}
          onPress={() => setExpanded((v) => !v)}
        >
          Filters
          {activeCount > 0 && (
            <Chip className="ml-1" color="primary" size="sm" variant="solid">
              {activeCount}
            </Chip>
          )}
        </Button>
      </div>

      {/* Secondary row: advanced filters */}
      {expanded && (
        <div className="flex flex-wrap items-center gap-2 rounded-xl bg-default-100/60 p-3">
          <Select
            aria-label="Species"
            className="w-40"
            isDisabled={speciesOpts.length === 0}
            placeholder="Species"
            selectedKeys={filters.species ? [filters.species] : []}
            size="sm"
            onChange={(e) => set("species", e.target.value)}
          >
            {speciesOpts.map((s) => (
              <SelectItem key={s}>{s}</SelectItem>
            ))}
          </Select>

          <Select
            aria-label="Tissue"
            className="w-40"
            isDisabled={tissueOpts.length === 0}
            placeholder="Tissue"
            selectedKeys={filters.tissue ? [filters.tissue] : []}
            size="sm"
            onChange={(e) => set("tissue", e.target.value)}
          >
            {tissueOpts.map((s) => (
              <SelectItem key={s}>{s}</SelectItem>
            ))}
          </Select>

          <Select
            aria-label="Platform"
            className="w-40"
            isDisabled={platformOpts.length === 0}
            placeholder="Platform"
            selectedKeys={filters.platform ? [filters.platform] : []}
            size="sm"
            onChange={(e) => set("platform", e.target.value)}
          >
            {platformOpts.map((s) => (
              <SelectItem key={s}>{s}</SelectItem>
            ))}
          </Select>

          <Select
            aria-label="Tag"
            className="w-36"
            isDisabled={tagOpts.length === 0}
            placeholder="Tag"
            selectedKeys={filters.tag ? [filters.tag] : []}
            size="sm"
            onChange={(e) => set("tag", e.target.value)}
          >
            {tagOpts.map((s) => (
              <SelectItem key={s}>{s}</SelectItem>
            ))}
          </Select>

          <Select
            aria-label="Status"
            className="w-36"
            selectedKeys={[filters.status]}
            size="sm"
            onChange={(e) => set("status", e.target.value || "all")}
          >
            {[
              <SelectItem key="all">All statuses</SelectItem>,
              ...statusOpts.map((s) => (
                <SelectItem key={s}>{s.toLowerCase()}</SelectItem>
              )),
            ]}
          </Select>

          <Input
            aria-label="Minimum genes"
            className="w-32"
            min={0}
            placeholder="Min genes"
            size="sm"
            type="number"
            value={filters.minGenes}
            onValueChange={(v) => set("minGenes", v)}
          />

          <Input
            aria-label="Added from"
            className="w-40"
            label="From"
            labelPlacement="outside-left"
            size="sm"
            type="date"
            value={filters.from}
            onValueChange={(v) => set("from", v)}
          />
          <Input
            aria-label="Added to"
            className="w-40"
            label="To"
            labelPlacement="outside-left"
            size="sm"
            type="date"
            value={filters.to}
            onValueChange={(v) => set("to", v)}
          />

          <Button
            className="ml-auto"
            size="sm"
            variant="light"
            onPress={() =>
              onChange({
                ...DEFAULT_FILTERS,
                sort: filters.sort,
                dir: filters.dir,
                search: filters.search,
              })
            }
          >
            Clear filters
          </Button>
        </div>
      )}

      <p className="text-xs text-default-400">
        Showing {shown} of {datasets.length} dataset
        {datasets.length === 1 ? "" : "s"}
      </p>
    </div>
  );
}

/** Apply the active filters + sort to the dataset list. Pure, client-side. */
export function applyDatasetFilters(
  datasets: DatasetRow[],
  f: DatasetFilters,
): DatasetRow[] {
  const term = f.search.trim().toLowerCase();
  const minGenes = f.minGenes ? Number(f.minGenes) : null;
  const fromTs = f.from ? new Date(f.from).getTime() : null;
  // Include the whole "to" day.
  const toTs = f.to ? new Date(f.to).getTime() + 86_400_000 : null;

  const filtered = datasets.filter((d) => {
    if (f.type !== "all" && d.datasetType !== f.type) return false;
    if (f.species && d.species !== f.species) return false;
    if (f.tissue && d.tissue !== f.tissue) return false;
    if (f.platform && d.platform !== f.platform) return false;
    if (f.status !== "all" && d.status !== f.status) return false;
    if (f.tag && !(d.tags ?? []).includes(f.tag)) return false;
    if (minGenes != null && !Number.isNaN(minGenes) && d.numGenes < minGenes)
      return false;

    const created = new Date(d.createdAt).getTime();

    if (fromTs != null && created < fromTs) return false;
    if (toTs != null && created >= toTs) return false;

    if (term) {
      const hay = [
        d.title,
        d.species,
        d.tissue,
        d.platform,
        d.disease,
        d.institute,
        d.description,
        ...(d.tags ?? []),
      ]
        .filter(Boolean)
        .join(" ")
        .toLowerCase();

      if (!hay.includes(term)) return false;
    }

    return true;
  });

  const dir = f.dir === "asc" ? 1 : -1;

  return filtered.sort((a, b) => {
    switch (f.sort) {
      case "name":
        return dir * (a.title ?? "").localeCompare(b.title ?? "");
      case "cells":
        return dir * (a.numCells - b.numCells);
      case "genes":
        return dir * (a.numGenes - b.numGenes);
      case "views":
        return dir * (a.viewCount - b.viewCount);
      case "date":
      default:
        return (
          dir *
          (new Date(a.createdAt).getTime() - new Date(b.createdAt).getTime())
        );
    }
  });
}
