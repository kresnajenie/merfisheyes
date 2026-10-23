import type { StandardizedDataset } from "../StandardizedDataset";

import {
  CONVENTIONS,
  NOT_DEFAULT_COLUMNS,
} from "@/lib/spatialdata/conventions";

/**
 * Determines the best cluster column to use for visualization
 * Priority:
 * 1. The priority columns, in order (taxonomy names, "leiden", "Cluster")
 * 2. Any column containing "celltype"
 * 3. First available categorical column
 * 4. First available column
 * 5. null if no clusters exist, or only positional ones
 *
 * The priority list and the positional columns (a Visium spot's grid row
 * and column, a MALDI spot's pixel in the ion image: offered in the picker,
 * never the default colouring) are in docs/specs/conventions.json.
 */
const PRIORITY_COLUMNS = CONVENTIONS.tables.defaultColumnPriority;

export function selectBestClusterColumn(
  dataset: StandardizedDataset | null,
): string | null {
  if (!dataset || !dataset.clusters || dataset.clusters.length === 0) {
    return null;
  }

  const clusters = dataset.clusters;
  const columnSet = new Set(clusters.map((c) => c.column));

  // Priority 1-2: Check priority columns in order
  for (const col of PRIORITY_COLUMNS) {
    if (columnSet.has(col)) return col;
  }

  // Priority 3: Find anything with "celltype" in it
  const celltype = clusters.find((c) =>
    c.column.toLowerCase().includes("celltype"),
  );

  if (celltype) return celltype.column;

  // Priority 4: Use the first categorical column
  const candidates = clusters.filter((c) => !NOT_DEFAULT_COLUMNS.has(c.column));
  const categoricalColumn = candidates.find((c) => c.type === "categorical");

  if (categoricalColumn) return categoricalColumn.column;

  // Priority 5: Use the first available column (none: one neutral colour)
  return candidates[0]?.column ?? null;
}

/**
 * Lightweight version that works on column names + types without needing loaded values.
 * Used during deferred cluster loading to pick the priority column before values are fetched.
 */
export function selectBestClusterColumnByName(
  columnNames: string[],
  columnTypes: Record<string, string>,
): string | null {
  if (!columnNames || columnNames.length === 0) return null;

  // Priority 1-2: Check priority columns in order
  for (const col of PRIORITY_COLUMNS) {
    if (columnNames.includes(col)) return col;
  }

  // Priority 3: contains "celltype"
  const celltype = columnNames.find((name) =>
    name.toLowerCase().includes("celltype"),
  );

  if (celltype) return celltype;

  // Priority 4: first categorical column
  const candidates = columnNames.filter((n) => !NOT_DEFAULT_COLUMNS.has(n));
  const categoricalColumn = candidates.find(
    (name) => columnTypes[name] !== "numerical",
  );

  if (categoricalColumn) return categoricalColumn;

  // Priority 5: first available column (none: one neutral colour)
  return candidates[0] ?? null;
}
