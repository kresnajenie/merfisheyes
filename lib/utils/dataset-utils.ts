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

  // Priority 1: Check for "leiden"
  const leiden = clusters.find((c) => c.column === "leiden");

  if (leiden) return leiden.column;
  const categoryCluster = clusters.find((c) => c.column === "Cluster");

  if (categoryCluster) return categoryCluster.column;
  // Priority 2: Find anything with "celltype" in it
  const celltype = clusters.find((c) =>
    c.column.toLowerCase().includes("celltype"),
  );

  if (celltype) return celltype.column;

  // Priority 3: Use the first categorical column
  const categoricalColumn = clusters.find((c) => c.type === "categorical");

  if (categoricalColumn) return categoricalColumn.column;

  // Priority 4: Use the first available column
  return clusters[0].column;
}
