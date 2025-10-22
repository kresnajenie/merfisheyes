import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { SingleMoleculeDataset } from "@/lib/SingleMoleculeDataset";

/**
 * Generate a content-based fingerprint for a dataset
 * Based on: file contents + dataset structure (NOT name/timestamp)
 *
 * Uses sampling for efficiency:
 * - Every 5th cell for coordinates
 * - Every 3rd gene for gene names
 * - Sample intervals for file contents
 */
export async function generateDatasetFingerprint(
  dataset: StandardizedDataset,
): Promise<string> {
  const parts: string[] = [];

  // 1. Dataset structure
  parts.push(`cells:${dataset.getPointCount()}`);
  parts.push(`genes:${dataset.genes.length}`);
  parts.push(`type:${dataset.type}`);

  // 2. Sample gene names (every 3rd gene)
  const sampledGenes: string[] = [];

  for (let i = 0; i < dataset.genes.length; i += 3) {
    sampledGenes.push(dataset.genes[i]);
  }
  parts.push(`genes:${sampledGenes.join(",")}`);

  // 3. Sample spatial coordinates (every 5th cell)
  if (dataset.spatial && dataset.spatial.coordinates) {
    const coords = dataset.spatial.coordinates;
    const sampledCoords: string[] = [];

    for (let i = 0; i < coords.length; i += 5) {
      const cellCoords = coords[i];

      sampledCoords.push(cellCoords.join(","));
    }
    parts.push(`spatial:${sampledCoords.join(";")}`);
  }

  // 4. Sample UMAP coordinates if available (every 5th cell)
  if (dataset.embeddings && dataset.embeddings.umap) {
    const coords = dataset.embeddings.umap;
    const sampledCoords: string[] = [];

    for (let i = 0; i < coords.length; i += 5) {
      const cellCoords = coords[i];

      sampledCoords.push(cellCoords.join(","));
    }
    parts.push(`umap:${sampledCoords.join(";")}`);
  }

  // 5. Sample expression data (every 5th cell, every 3rd gene)
  const sampledExpression: string[] = [];

  for (let geneIdx = 0; geneIdx < dataset.genes.length; geneIdx += 3) {
    const geneName = dataset.genes[geneIdx];
    const values = await dataset.getGeneExpression(geneName);

    if (values) {
      const sampledValues: number[] = [];

      for (let cellIdx = 0; cellIdx < values.length; cellIdx += 5) {
        sampledValues.push(values[cellIdx]);
      }
      sampledExpression.push(`${geneName}:${sampledValues.join(",")}`);
    }
  }
  parts.push(`expr:${sampledExpression.join(";")}`);

  // 6. Cluster metadata keys (if available)
  if (dataset.clusters && dataset.clusters.length > 0) {
    const clusterColumns = dataset.clusters
      .map((cluster) => cluster.column)
      .sort();

    parts.push(`clusters:${clusterColumns.join(",")}`);
  }

  // Combine all parts and hash
  const combined = parts.join("|");

  return await hashString(combined);
}

/**
 * Generate a content-based fingerprint for a single molecule dataset
 * Based on: gene names + per-gene molecule counts
 *
 * Uses gene names and molecule counts to uniquely identify the dataset
 */
export async function generateSingleMoleculeFingerprint(
  dataset: SingleMoleculeDataset,
): Promise<string> {
  const parts: string[] = [];

  // 1. Dataset structure
  parts.push(`molecules:${dataset.getMoleculeCount()}`);
  parts.push(`genes:${dataset.uniqueGenes.length}`);
  parts.push(`type:${dataset.type}`);
  parts.push(`dimensions:${dataset.dimensions}`);

  // 2. Gene names (sorted for consistency)
  const sortedGenes = [...dataset.uniqueGenes].sort();

  parts.push(`gene_names:${sortedGenes.join(",")}`);

  // 3. Per-gene molecule counts (sorted by gene name for consistency)
  const geneCounts: string[] = [];

  for (const gene of sortedGenes) {
    const coords = dataset.getCoordinatesByGene(gene);
    const moleculeCount = coords.length / dataset.dimensions; // coordinates are flat array

    geneCounts.push(`${gene}:${moleculeCount}`);
  }
  parts.push(`gene_counts:${geneCounts.join(";")}`);

  // Combine all parts and hash
  const combined = parts.join("|");

  return await hashString(combined);
}

/**
 * Simple hash function using Web Crypto API
 */
async function hashString(str: string): Promise<string> {
  const encoder = new TextEncoder();
  const data = encoder.encode(str);
  const hashBuffer = await crypto.subtle.digest("SHA-256", data);
  const hashArray = Array.from(new Uint8Array(hashBuffer));
  const hashHex = hashArray
    .map((b) => b.toString(16).padStart(2, "0"))
    .join("");

  return hashHex;
}
