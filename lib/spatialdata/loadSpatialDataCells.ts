import { SpatialDataTableAdapter } from "./SpatialDataTableAdapter";
import { type SpatialDataStore } from "./store";

import { StandardizedDataset } from "@/lib/StandardizedDataset";
import { selectBestClusterColumnByName } from "@/lib/utils/dataset-utils";

/**
 * Load the cells of a SpatialData zarr store as a StandardizedDataset. Only
 * the table's coordinates, gene list and one obs column are read up front;
 * gene expression is fetched per gene on demand.
 */
export async function loadSpatialDataCells(
  { store, keys, sidecar }: SpatialDataStore,
  datasetId: string,
  name: string,
  onProgress?: (progress: number, message: string) => void,
  // Folder holding the store, for stores opened by URL: the scene looks for
  // a mapping.json there to find the linked molecule dataset.
  customS3BaseUrl?: string,
): Promise<StandardizedDataset> {
  const adapter = new SpatialDataTableAdapter(store, keys, sidecar);

  await adapter.initialize(onProgress);

  onProgress?.(50, "Loading spatial coordinates...");
  const spatial = await adapter.loadSpatialCoordinates();

  onProgress?.(60, "Loading genes...");
  const genes = await adapter.loadGenes();

  onProgress?.(70, "Loading cell metadata...");
  const columnInfo = adapter.getClusterColumnInfo();
  const priorityColumn = selectBestClusterColumnByName(
    columnInfo.names,
    columnInfo.types,
  );
  const clusters = priorityColumn
    ? await adapter.loadClusters([priorityColumn])
    : null;
  const dataInfo = adapter.getDatasetInfo();

  const dataset = new StandardizedDataset({
    id: datasetId,
    name,
    type: "h5ad-zarr",
    spatial,
    embeddings: {},
    genes,
    clusters,
    metadata: {
      numCells: dataInfo.numCells,
      numGenes: dataInfo.numGenes,
      spatialDimensions: dataInfo.spatialDimensions,
      availableEmbeddings: dataInfo.availableEmbeddings,
      clusterCount: dataInfo.clusterCount,
      xFormat: dataInfo.xFormat,
      loadedFrom: "s3_zarr",
      ...(customS3BaseUrl && { customS3BaseUrl }),
    },
    adapter,
    rawData: null,
    normalized: false,
  });

  dataset.allClusterColumnNames = columnInfo.names;
  dataset.allClusterColumnTypes = columnInfo.types;
  dataset.clustersFullyLoaded = columnInfo.names.length <= 1;
  dataset.availableDeStatsColumns = adapter.getAvailableDeStatsColumns();
  dataset.allEmbeddingNames = dataInfo.availableEmbeddings;
  dataset.embeddingsFullyLoaded = false;

  onProgress?.(100, "Dataset loaded");

  return dataset;
}
