import type { StandardizedDataset } from "../StandardizedDataset";

/**
 * Determines the best cluster column to use for visualization
 * Priority:
 * 1. Taxonomy columns (class_name, subclass_name, supertype_name, etc.)
 * 2. "leiden"
 * 3. "Cluster"
 * 4. Any column containing "celltype"
 * 5. First available categorical column
 * 6. First available column
 * 7. null if no clusters exist
 */

const PRIORITY_COLUMNS = [
  "class_name",
  "subclass_name",
  "supertype_name",
  "cluster_name",
  "subcluster_name",
  "super_cluster_name",
  "leiden",
  "Cluster",
];

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
  const categoricalColumn = clusters.find((c) => c.type === "categorical");

  if (categoricalColumn) return categoricalColumn.column;

  // Priority 5: Use the first available column
  return clusters[0].column;
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
  const categoricalColumn = columnNames.find(
    (name) => columnTypes[name] !== "numerical",
  );
  if (categoricalColumn) return categoricalColumn;

  // Priority 5: first available column
  return columnNames[0];
}
