import { pointSource, readCellCentroids } from "./centroids";
import { readCellFootprint } from "./footprint";
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
  sd: SpatialDataStore,
  datasetId: string,
  name: string,
  onProgress?: (progress: number, message: string) => void,
  // Folder holding the store, for stores opened by URL: the scene looks for
  // a mapping.json there to find the linked molecule dataset.
  customS3BaseUrl?: string,
): Promise<StandardizedDataset> {
  const { store, keys, sidecar } = sd;
  const adapter = new SpatialDataTableAdapter(store, keys, sidecar);

  await adapter.initialize(onProgress);

  onProgress?.(50, "Loading spatial coordinates...");
  // The cells sit where the region's Point element puts them (MALDI's
  // obsm/spatial is the ion image's pixel index, 25-50x off); with none,
  // obsm/spatial; with neither, the polygons' centres (CosMx)
  const points = await pointSource(sd).catch(() => null);
  const fromObsm = () =>
    adapter.hasSpatialCoordinates()
      ? adapter.loadSpatialCoordinates()
      : readCellCentroids(sd, adapter);
  const spatial = points
    ? await readCellCentroids(sd, adapter, points).catch((error) => {
        console.warn(`Cells not placed from ${points}:`, error);

        return fromObsm();
      })
    : await fromObsm();

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
  // Spots and bins are drawn at their real size; a store that cannot say
  // keeps its dots
  dataset.cellFootprint = await readCellFootprint(
    sd,
    spatial.coordinates,
    spatial.dimensions,
  ).catch((error) => {
    console.warn("Cell footprint not read:", error);

    return null;
  });
  dataset.allEmbeddingNames = dataInfo.availableEmbeddings;
  dataset.embeddingsFullyLoaded = false;

  onProgress?.(100, "Dataset loaded");

  return dataset;
}
