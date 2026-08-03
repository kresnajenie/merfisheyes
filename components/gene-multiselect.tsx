"use client";

import { Input } from "@heroui/input";
import { useMemo, useState } from "react";

/**
 * Compact searchable gene multi-select (checkbox list). Controlled — the parent
 * owns the selected `value` and how it's saved. Caps the selection at `max`.
 * Used for owner "default genes" pickers (pure SM viewer + SC+SM overlay).
 */
export function GeneMultiSelect({
  genes,
  value,
  onChange,
  label = "Genes",
  max = 20,
}: {
  genes: string[];
  value: string[];
  onChange: (v: string[]) => void;
  label?: string;
  max?: number;
}) {
  const [query, setQuery] = useState("");

  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase();
    const list = q ? genes.filter((g) => g.toLowerCase().includes(q)) : genes;

    return list.slice(0, 200);
  }, [query, genes]);

  const toggle = (g: string) => {
    if (value.includes(g)) onChange(value.filter((x) => x !== g));
    else if (value.length < max) onChange([...value, g]);
  };

  return (
    <div>
      <p className="text-[10px] text-default-400 mb-1">
        {label} ({value.length})
      </p>
      <Input
        className="mb-1"
        placeholder="Search genes"
        size="sm"
        value={query}
        onValueChange={setQuery}
      />
      <div className="max-h-36 overflow-y-auto rounded-md border border-white/10 p-1">
        {filtered.map((g) => (
          <label
            key={g}
            className="flex items-center gap-2 text-xs py-0.5 px-1 cursor-pointer hover:bg-white/5 rounded"
          >
            <input
              checked={value.includes(g)}
              type="checkbox"
              onChange={() => toggle(g)}
            />
            <span className="truncate">{g}</span>
          </label>
        ))}
      </div>
    </div>
  );
}
