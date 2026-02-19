import type { StandardizedDataset } from "../StandardizedDataset";

/**
 * Determines the best cluster column to use for visualization
 * Priority:
 * 1. "leiden" if it exists
 * 2. Any column containing "celltype"
 * 3. First available categorical column
 * 4. First available column
 * 5. null if no clusters exist
 */
export function selectBestClusterColumn(
  dataset: StandardizedDataset | null,
): string | null {
  if (!dataset || !dataset.clusters || dataset.clusters.length === 0) {
    return null;
  }

  const clusters = dataset.clusters;

  // Priority 1: Check for "class_name"
  const className = clusters.find((c) => c.column === "class_name");

  if (className) return className.column;

  // Priority 2: Check for "leiden"
  const leiden = clusters.find((c) => c.column === "leiden");

  if (leiden) return leiden.column;
  const categoryCluster = clusters.find((c) => c.column === "Cluster");

  if (categoryCluster) return categoryCluster.column;
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

  // Priority 1: "class_name"
  if (columnNames.includes("class_name")) return "class_name";

  // Priority 2: "leiden"
  if (columnNames.includes("leiden")) return "leiden";

  // Priority 3: "Cluster"
  if (columnNames.includes("Cluster")) return "Cluster";

  // Priority 4: contains "celltype"
  const celltype = columnNames.find((name) =>
    name.toLowerCase().includes("celltype"),
  );
  if (celltype) return celltype;

  // Priority 5: first categorical column
  const categoricalColumn = columnNames.find(
    (name) => columnTypes[name] !== "numerical",
  );
  if (categoricalColumn) return categoricalColumn;

  // Priority 6: first available column
  return columnNames[0];
}
