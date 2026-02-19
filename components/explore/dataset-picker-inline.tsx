"use client";

import { useState, useEffect, useCallback } from "react";
import { Input } from "@heroui/input";
import { Spinner } from "@heroui/spinner";
import { ScrollShadow } from "@heroui/scroll-shadow";
import { SearchIcon } from "@/components/icons";

import { CompactDatasetCard } from "./compact-dataset-card";
import type { CatalogDatasetItem } from "./types";

interface DatasetPickerInlineProps {
  /** Called when user selects a dataset */
  onSelect: (dataset: CatalogDatasetItem) => void;
}

export function DatasetPickerInline({ onSelect }: DatasetPickerInlineProps) {
  const [search, setSearch] = useState("");
  const [items, setItems] = useState<CatalogDatasetItem[]>([]);
  const [loading, setLoading] = useState(true);

  const fetchData = useCallback(async () => {
    setLoading(true);
    const params = new URLSearchParams({ limit: "20" });
    if (search) params.set("search", search);

    const res = await fetch(`/api/explore?${params}`);
    const data = await res.json();
    setItems(data.items ?? []);
    setLoading(false);
  }, [search]);

  useEffect(() => {
    const timer = setTimeout(fetchData, 300);
    return () => clearTimeout(timer);
  }, [fetchData]);

  return (
    <div className="flex flex-col gap-3">
      <Input
        classNames={{ inputWrapper: "bg-default-100" }}
        placeholder="Search datasets..."
        size="sm"
        startContent={
          <SearchIcon className="text-default-400 pointer-events-none flex-shrink-0 w-4 h-4" />
        }
        value={search}
        onValueChange={setSearch}
      />

      <ScrollShadow className="max-h-80">
        {loading ? (
          <div className="flex justify-center py-6">
            <Spinner size="sm" />
          </div>
        ) : items.length === 0 ? (
          <p className="text-sm text-default-500 text-center py-6">
            No published datasets found
          </p>
        ) : (
          <div className="flex flex-col gap-2">
            {items.map((item) => (
              <CompactDatasetCard
                key={item.id}
                dataset={item}
                onSelect={onSelect}
              />
            ))}
          </div>
        )}
      </ScrollShadow>
    </div>
  );
}
