import type { StandardizedDataset } from "@/lib/StandardizedDataset";

// The browser single-cell upload's manifest + file map, extracted from
// components/upload-settings-modal.tsx so the parity harness (scripts/parity/)
// can produce exactly what the modal uploads, headlessly.

/**
 * Create manifest.json
 */
export async function createManifest(
  dataset: StandardizedDataset,
  datasetId: string,
  datasetName: string,
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  obsMetadata: Record<string, any>,
  deStatsColumns: string[] = [],
): Promise<string> {
  // Determine spatial dimensions
  let spatialDimensions = 2;

  if (
    dataset.spatial &&
    dataset.spatial.coordinates &&
    dataset.spatial.coordinates.length > 0
  ) {
    spatialDimensions = dataset.spatial.dimensions;
  }

  // Get available embeddings (prefer allEmbeddingNames for lazy-loaded datasets)
  const availableEmbeddings: string[] =
    dataset.allEmbeddingNames && dataset.allEmbeddingNames.length > 0
      ? [...dataset.allEmbeddingNames]
      : Object.keys(dataset.embeddings ?? {});

  // Count clusters
  let clusterCount = 0;

  if (dataset.clusters && dataset.clusters.length > 0) {
    clusterCount = new Set(dataset.clusters).size;
  }

  const manifest = {
    version: "2.0",
    normalized: dataset.normalized,
    created_at: new Date().toISOString(),
    dataset_id: datasetId,
    name: datasetName,
    type: dataset.type,
    statistics: {
      total_cells: dataset.getPointCount(),
      total_genes: dataset.genes.length,
      spatial_dimensions: spatialDimensions,
      available_embeddings: availableEmbeddings,
      cluster_count: clusterCount,
    },
    files: {
      coordinates: Object.keys(coordinates),
      expression: {
        num_chunks: chunks.length,
        chunk_size: index.chunk_size,
        total_genes: dataset.genes.length,
      },
      // Keys read by consumers (same names as the server pipeline writes).
      expression_chunks: chunks.length,
      observation_columns: Object.keys(obsMetadata),
      observations: Object.keys(obsMetadata),
      palettes: [], // TODO: Add palette support
      de_stats: deStatsColumns,
    },
    processing: {
      chunk_strategy: "adaptive",
      compression: "gzip",
      sparse_format: true,
      created_by: "MERFISH Visualizer",
      source_file: dataset.name || "unknown",
    },
  };

  return JSON.stringify(manifest, null, 2);
}

/**
 * Prepare files for upload with proper structure
 */
export async function prepareFilesForUpload(
  chunks: any[],
  index: any,
  coordinates: Record<string, Blob>,
  observationFiles: Record<string, Blob>,
  observationMetadata: Record<string, any>,
  paletteFiles: Record<string, Blob>,
  deStatsFiles: Record<string, Blob>,
  manifestJson: string,
): Promise<{ key: string; blob: Blob; size: number; contentType: string }[]> {
  const files: {
    key: string;
    blob: Blob;
    size: number;
    contentType: string;
  }[] = [];

  // Manifest
  const manifestBlob = new Blob([manifestJson], { type: "application/json" });

  files.push({
    key: "manifest.json",
    blob: manifestBlob,
    size: manifestBlob.size,
    contentType: "application/json",
  });

  // Gene chunks
  for (const chunk of chunks) {
    files.push({
      key: `expr/${chunk.filename}`,
      blob: chunk.data,
      size: chunk.data.size,
      contentType: "application/gzip",
    });
  }

  // Gene index
  const indexBlob = new Blob([JSON.stringify(index, null, 2)], {
    type: "application/json",
  });

  files.push({
    key: "expr/index.json",
    blob: indexBlob,
    size: indexBlob.size,
    contentType: "application/json",
  });

  // Coordinates
  for (const [name, blob] of Object.entries(coordinates)) {
    files.push({
      key: `coords/${name}.bin.gz`,
      blob: blob,
      size: blob.size,
      contentType: "application/gzip",
    });
  }

  // Observations (binary obs format, see lib/utils/obs-binary.ts)
  for (const [name, blob] of Object.entries(observationFiles)) {
    files.push({
      key: `obs/${name}.bin.gz`,
      blob: blob,
      size: blob.size,
      contentType: "application/gzip",
    });
  }

  // Observation metadata
  const obsMetadataBlob = new Blob(
    [JSON.stringify(observationMetadata, null, 2)],
    { type: "application/json" },
  );

  files.push({
    key: "obs/metadata.json",
    blob: obsMetadataBlob,
    size: obsMetadataBlob.size,
    contentType: "application/json",
  });

  // Palettes
  for (const [name, blob] of Object.entries(paletteFiles)) {
    files.push({
      key: `palettes/${name}.json`,
      blob: blob,
      size: blob.size,
      contentType: "application/json",
    });
  }

  // DE stats
  for (const [name, blob] of Object.entries(deStatsFiles)) {
    files.push({
      key: `de/${name}.bin.gz`,
      blob: blob,
      size: blob.size,
      contentType: "application/gzip",
    });
  }

  return files;
}
