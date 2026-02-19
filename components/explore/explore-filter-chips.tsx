"use client";

import { Chip } from "@heroui/chip";

interface ActiveFilter {
  key: string;
  label: string;
  value: string;
  onClear: () => void;
}

interface ExploreFilterChipsProps {
  filters: ActiveFilter[];
}

export function ExploreFilterChips({ filters }: ExploreFilterChipsProps) {
  if (filters.length === 0) return null;

  return (
    <div className="flex flex-wrap gap-2 mb-4">
      {filters.map((f) => (
        <Chip
          key={f.key}
          size="sm"
          variant="flat"
          onClose={f.onClear}
        >
          {f.label}: {f.value}
        </Chip>
      ))}
    </div>
  );
}
