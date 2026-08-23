/**
 * Gene Chunk Processor
 * Processes gene expression data into optimized chunks for efficient storage and retrieval
 * Based on TECH_SPEC_UPLOAD.md
 */

import type { StandardizedDataset } from "../StandardizedDataset";

import { getClusterValue } from "../StandardizedDataset";
import {
  OBS_BINARY_FORMAT,
  encodeCategoricalObs,
  encodeNumericalObs,
} from "./obs-binary";
import { computeDeStats, encodeDeStatsBuffer } from "./de-stats";

/** Same threshold as --de-max-celltypes in scripts/process_spatial_data.py. */
const DE_MAX_CELLTYPES = 20000;

export interface GeneChunkMetadata {
  name: string;
  mean: number;
  percent_expressed: number;
  num_nonzero: number;
}

export interface ChunkInfo {
  id: number;
  filename: string;
  gene_range: [number, number];
  size_compressed: number;
  size_uncompressed: number;
  num_genes: number;
}

export interface GeneInfo {
  index: number;
  name: string;
  chunk_id: number;
  position_in_chunk: number;
  mean_expression: number;
  percent_cells_expressed: number;
  num_nonzero: number;
}

export interface GeneChunkIndex {
  version: number;
  created_at: string;
  total_genes: number;
  total_cells: number;
  chunk_size: number;
  num_chunks: number;
  chunks: ChunkInfo[];
  genes: GeneInfo[];
  gene_lookup: Record<string, number>;
}

export interface ProcessedChunk {
  filename: string;
  data: Blob;
  metadata: {
    genes: GeneChunkMetadata[];
  };
}

export interface ProcessedGenes {
  chunks: ProcessedChunk[];
  index: GeneChunkIndex;
}

export class GeneChunkProcessor {
  private chunkSize: number | null;

  constructor(chunkSize: number | null = null) {
    this.chunkSize = chunkSize;
  }

  /**
   * Genes per expression chunk. The spec for BOTH pipelines (browser and
   * server) is ONE gene per chunk: the viewer fetches exactly the gene it
   * needs, whatever the cell count. An explicit constructor value overrides
   * it for local experiments only.
   */
  determineChunkSize(_numGenes: number): number {
    if (this.chunkSize) return this.chunkSize;

    return 1;
  }

  /**
   * Process genes into chunks with binary format
   */
  async processGenes(
    dataset: StandardizedDataset,
    onProgress?: (progress: number, message: string) => void,
  ): Promise<ProcessedGenes> {
    const genes = dataset.genes;
    const numCells = dataset.getPointCount();
    const chunkSize = this.determineChunkSize(genes.length);

    console.log(`Processing ${genes.length} genes into chunks of ${chunkSize}`);

    // Load expression matrix once for all chunks
    onProgress?.(0, "Loading expression matrix...");

    // Use cached matrix if available (from worker), otherwise fetch via adapter
    let matrix;

    if (dataset.matrix) {
      console.log("Using cached expression matrix from worker");
      matrix = dataset.matrix;
    } else if (dataset.adapter) {
      console.log("Fetching expression matrix via adapter");
      matrix = dataset.adapter.fetchFullMatrix();
    } else {
      throw new Error("No expression matrix or adapter available");
    }

    console.log("Expression matrix loaded");

    const chunks: ProcessedChunk[] = [];
    const index: GeneChunkIndex = {
      version: 2,
      created_at: new Date().toISOString(),
      total_genes: genes.length,
      total_cells: numCells,
      chunk_size: chunkSize,
      num_chunks: Math.ceil(genes.length / chunkSize),
      chunks: [],
      genes: [],
      gene_lookup: {},
    };

    // Process genes in chunks
    for (let chunkId = 0; chunkId < index.num_chunks; chunkId++) {
      const startIdx = chunkId * chunkSize;
      const endIdx = Math.min(startIdx + chunkSize, genes.length);
      const genesInChunk = endIdx - startIdx;

      const progress = (chunkId / index.num_chunks) * 100;

      onProgress?.(
        progress,
        `Processing chunk ${chunkId + 1}/${index.num_chunks} (${genesInChunk} genes)`,
      );

      console.log(
        `Processing chunk ${chunkId}: genes ${startIdx}-${endIdx - 1} (${genesInChunk} genes)`,
      );

      // Create chunk data
      const chunkData = await this.createChunk(
        dataset,
        genes,
        matrix,
        startIdx,
        endIdx,
        chunkId,
        numCells,
      );

      // Compress chunk
      const compressed = await this.compress(chunkData.buffer);

      chunks.push({
        filename: `chunk_${chunkId.toString().padStart(5, "0")}.bin.gz`,
        data: compressed,
        metadata: chunkData.metadata,
      });

      // Update index
      index.chunks.push({
        id: chunkId,
        filename: `chunk_${chunkId.toString().padStart(5, "0")}.bin.gz`,
        gene_range: [startIdx, endIdx - 1],
        size_compressed: compressed.size,
        size_uncompressed: chunkData.buffer.byteLength,
        num_genes: genesInChunk,
      });

      // Add gene entries to index
      for (let i = 0; i < genesInChunk; i++) {
        const geneIdx = startIdx + i;
        const geneName = genes[geneIdx];

        index.genes.push({
          index: geneIdx,
          name: geneName,
          chunk_id: chunkId,
          position_in_chunk: i,
          mean_expression: chunkData.metadata.genes[i].mean,
          percent_cells_expressed:
            chunkData.metadata.genes[i].percent_expressed,
          num_nonzero: chunkData.metadata.genes[i].num_nonzero,
        });

        index.gene_lookup[geneName] = geneIdx;
      }
    }

    onProgress?.(100, "Gene processing complete");

    return { chunks, index };
  }

  /**
   * Create binary chunk data for a range of genes
   */
  private async createChunk(
    dataset: StandardizedDataset,
    geneNames: string[],
    matrix: any,
    startIdx: number,
    endIdx: number,
    chunkId: number,
    numCells: number,
  ): Promise<{
    buffer: ArrayBuffer;
    metadata: { genes: GeneChunkMetadata[] };
  }> {
    const numGenes = endIdx - startIdx;

    // Create header (16 bytes)
    const header = new ArrayBuffer(16);
    const headerView = new DataView(header);

    headerView.setUint32(0, 1, true); // version
    headerView.setUint32(4, numGenes, true); // num_genes
    headerView.setUint32(8, chunkId, true); // chunk_id
    headerView.setUint32(12, numCells, true); // total_cells

    // Gene table (24 bytes per gene)
    const geneTable = new ArrayBuffer(numGenes * 24);
    const geneTableView = new DataView(geneTable);

    // Process each gene
    const geneDataBuffers: ArrayBuffer[] = [];
    const metadata: { genes: GeneChunkMetadata[] } = { genes: [] };
    let currentOffset = 16 + numGenes * 24; // After header and table

    for (let i = 0; i < numGenes; i++) {
      const geneIdx = startIdx + i;
      const geneName = geneNames[geneIdx];

      // Extract gene column from matrix directly
      let geneData: number[];

      // Use dataset's method if adapter is available, otherwise extract directly
      if (dataset.adapter) {
        geneData = dataset.adapter.fetchColumn(matrix, geneIdx);
      } else {
        // Extract column directly from cached matrix (same logic as extractColumnFromMatrix)
        geneData = this.extractColumnFromMatrix(
          matrix,
          geneIdx,
          geneName,
          geneNames,
        );
      }

      if (!geneData || geneData.length === 0) {
        throw new Error(`Failed to get expression data for gene: ${geneName}`);
      }

      // Create sparse representation
      const geneSparseData = this.extractAndEncodeSparseGene(
        geneData,
        numCells,
      );

      // Write gene table entry
      const tableOffset = i * 24;

      geneTableView.setUint32(tableOffset, geneIdx, true);
      geneTableView.setUint32(tableOffset + 4, currentOffset, true);
      geneTableView.setUint32(tableOffset + 8, geneSparseData.byteLength, true);
      geneTableView.setUint32(
        tableOffset + 12,
        geneSparseData.uncompressedSize,
        true,
      );
      geneTableView.setUint32(
        tableOffset + 16,
        geneSparseData.numNonZero,
        true,
      );
      geneTableView.setUint32(tableOffset + 20, 0, true); // reserved

      geneDataBuffers.push(geneSparseData.buffer);
      currentOffset += geneSparseData.byteLength;

      // Collect metadata
      metadata.genes.push({
        name: geneName,
        mean: geneSparseData.meanExpression,
        percent_expressed: (geneSparseData.numNonZero / numCells) * 100,
        num_nonzero: geneSparseData.numNonZero,
      });
    }

    // Combine all buffers
    const totalSize =
      16 +
      numGenes * 24 +
      geneDataBuffers.reduce((sum, buf) => sum + buf.byteLength, 0);
    const chunkBuffer = new ArrayBuffer(totalSize);
    const chunkArray = new Uint8Array(chunkBuffer);

    // Copy header
    chunkArray.set(new Uint8Array(header), 0);

    // Copy gene table
    chunkArray.set(new Uint8Array(geneTable), 16);

    // Copy gene data
    let offset = 16 + numGenes * 24;

    for (const geneBuffer of geneDataBuffers) {
      chunkArray.set(new Uint8Array(geneBuffer), offset);
      offset += geneBuffer.byteLength;
    }

    return {
      buffer: chunkBuffer,
      metadata: metadata,
    };
  }

  /**
   * Convert dense gene data to sparse binary format
   */
  private extractAndEncodeSparseGene(
    geneData: number[],
    numCells: number,
  ): {
    buffer: ArrayBuffer;
    byteLength: number;
    uncompressedSize: number;
    numNonZero: number;
    meanExpression: number;
  } {
    // Find non-zero values
    const nonZeroIndices: number[] = [];
    const nonZeroValues: number[] = [];
    let sum = 0;

    for (let i = 0; i < geneData.length; i++) {
      if (geneData[i] !== 0 && !isNaN(geneData[i])) {
        nonZeroIndices.push(i);
        nonZeroValues.push(geneData[i]);
        sum += geneData[i];
      }
    }

    const numNonZero = nonZeroIndices.length;
    const meanExpression = numNonZero > 0 ? sum / numCells : 0;

    // Create sparse buffer
    const bufferSize = 8 + numNonZero * 8; // header + indices + values
    const buffer = new ArrayBuffer(bufferSize);
    const view = new DataView(buffer);

    // Write header
    view.setUint32(0, numCells, true);
    view.setUint32(4, numNonZero, true);

    // Write indices
    let offset = 8;

    for (const idx of nonZeroIndices) {
      view.setUint32(offset, idx, true);
      offset += 4;
    }

    // Write values
    for (const val of nonZeroValues) {
      view.setFloat32(offset, val, true);
      offset += 4;
    }

    return {
      buffer: buffer,
      byteLength: bufferSize,
      uncompressedSize: numCells * 4, // If it were dense
      numNonZero: numNonZero,
      meanExpression: meanExpression,
    };
  }

  /**
   * Extract a column from a matrix without requiring an adapter
   * Handles different matrix formats (Map, Array, TypedArray)
   */
  private extractColumnFromMatrix(
    matrix: any,
    columnIndex: number,
    geneName: string,
    allGenes: string[],
  ): number[] {
    // Case 1: Map<string, Float32Array> (Xenium/MERSCOPE format)
    if (matrix instanceof Map) {
      const gene = allGenes[columnIndex];

      if (!gene || !matrix.has(gene)) {
        return [];
      }

      return Array.from(matrix.get(gene)!);
    }

    // Case 2: Array of arrays (row-major)
    if (Array.isArray(matrix) && Array.isArray(matrix[0])) {
      return matrix.map((row: any) => row[columnIndex]);
    }

    // Case 3: TypedArray (flattened row-major - H5AD format)
    if (ArrayBuffer.isView(matrix)) {
      const typedArray = matrix as any;
      const numCells = typedArray.length / allGenes.length;
      const numGenes = allGenes.length;

      if (columnIndex >= numGenes) {
        throw new Error("Column index out of bounds");
      }

      return Array.from(
        { length: numCells },
        (_, i) => typedArray[i * numGenes + columnIndex],
      );
    }

    throw new Error("Unsupported matrix format");
  }

  /**
   * Compress binary data using gzip
   */
  private async compress(buffer: ArrayBuffer): Promise<Blob> {
    const blob = new Blob([buffer]);
    const stream = blob.stream().pipeThrough(new CompressionStream("gzip"));

    return new Response(stream).blob();
  }

  /**
   * Process coordinates into binary format
   */
  async processCoordinates(
    dataset: StandardizedDataset,
  ): Promise<Record<string, Blob>> {
    const coords: Record<string, Blob> = {};

    console.log("Processing coordinates");

    // Process spatial coordinates
    if (dataset.spatial && dataset.spatial.coordinates) {
      coords.spatial = await this.encodeCoordinates(
        dataset.spatial.coordinates,
        dataset.spatial.dimensions,
      );
    }

    // Process embeddings
    if (dataset.embeddings) {
      for (const [name, coordinates] of Object.entries(dataset.embeddings)) {
        if (coordinates && coordinates.length > 0) {
          coords[name] = await this.encodeCoordinates(coordinates);
        }
      }
    }

    return coords;
  }

  /**
   * Encode coordinates as binary format
   */
  private async encodeCoordinates(
    coordinates: Float32Array | number[][],
    dimensionsHint?: number,
  ): Promise<Blob> {
    const isFlat = coordinates instanceof Float32Array;
    const dimensions = isFlat
      ? (dimensionsHint || 3)
      : (coordinates as number[][])[0]?.length || 0;
    const numPoints = isFlat
      ? coordinates.length / dimensions
      : (coordinates as number[][]).length;

    // Create buffer: header + coordinate data
    const bufferSize = 8 + numPoints * dimensions * 4;
    const buffer = new ArrayBuffer(bufferSize);
    const view = new DataView(buffer);

    // Write header
    view.setUint32(0, numPoints, true);
    view.setUint32(4, dimensions, true);

    // Write coordinates
    let offset = 8;

    if (isFlat) {
      for (let i = 0; i < coordinates.length; i++) {
        view.setFloat32(offset, coordinates[i], true);
        offset += 4;
      }
    } else {
      const nested = coordinates as number[][];

      for (let i = 0; i < numPoints; i++) {
        for (let d = 0; d < dimensions; d++) {
          view.setFloat32(offset, nested[i][d], true);
          offset += 4;
        }
      }
    }

    // Compress
    return await this.compress(buffer);
  }

  /**
   * Process observations into separate files
   */
  async processObservations(dataset: StandardizedDataset): Promise<{
    files: Record<string, Blob>;
    metadata: Record<string, any>;
  }> {
    const files: Record<string, Blob> = {};
    const metadata: Record<string, any> = {};

    if (!dataset.clusters || dataset.clusters.length === 0) {
      return { files, metadata };
    }

    for (const cluster of dataset.clusters) {
      // Binary obs format (lib/utils/obs-binary.ts) — identical to what the
      // server pipeline writes: sorted dictionary + per-cell codes for
      // categorical columns, float32 per cell for numerical ones.
      const type = cluster.type === "numerical" ? "numerical" : "categorical";
      const n = cluster.valueIndices?.length ?? cluster.values.length;
      let payload: Uint8Array;
      let uniqueCount: number;

      if (type === "numerical") {
        const values = new Float32Array(n);
        const seen = new Set<number>();

        for (let i = 0; i < n; i++) {
          const raw = getClusterValue(cluster, i);
          const num = raw === "" || raw === "null" ? NaN : Number(raw);

          values[i] = num;
          seen.add(Number.isNaN(num) ? NaN : num);
        }
        payload = encodeNumericalObs(values);
        uniqueCount = seen.size;
      } else {
        // Dictionary in sorted label order; missing values are the label "".
        const labelOf = (i: number) => {
          const v = getClusterValue(cluster, i);

          return v == null || v === "null" || v === "undefined" ? "" : v;
        };
        const labels = new Set<string>();

        for (let i = 0; i < n; i++) labels.add(labelOf(i));
        const dict = Array.from(labels).sort();
        const codeOf = new Map(dict.map((label, idx) => [label, idx]));
        const Codes =
          dict.length <= 0xff ? Uint8Array : dict.length <= 0xffff ? Uint16Array : Uint32Array;
        const codes = new Codes(n);

        for (let i = 0; i < n; i++) codes[i] = codeOf.get(labelOf(i))!;
        payload = encodeCategoricalObs(dict, codes);
        uniqueCount = dict.length;
      }

      files[cluster.column] = await this.compress(
        payload.buffer.slice(payload.byteOffset, payload.byteOffset + payload.byteLength) as ArrayBuffer,
      );

      metadata[cluster.column] = {
        type,
        unique_values: uniqueCount,
        format: OBS_BINARY_FORMAT,
      };
    }

    return { files, metadata };
  }

  /**
   * Per-celltype DE stats (de/{col}.bin.gz) for every categorical column —
   * the same files the server pipeline writes, so browser-uploaded datasets
   * get the DEG panel too. Columns with more than DE_MAX_CELLTYPES values
   * are skipped (ID-like columns misclassified as categorical).
   */
  async processDeStats(
    dataset: StandardizedDataset,
    onProgress?: (progress: number, message: string) => void,
  ): Promise<Record<string, Blob>> {
    const files: Record<string, Blob> = {};
    const matrix = dataset.matrix;

    if (!matrix || !dataset.clusters || dataset.clusters.length === 0) {
      return files;
    }

    const categorical = dataset.clusters.filter(
      (c) => c.type === "categorical" && c.valueIndices && c.uniqueValues,
    );

    for (let i = 0; i < categorical.length; i++) {
      const cluster = categorical[i];

      if ((cluster.uniqueValues?.length ?? 0) > DE_MAX_CELLTYPES) continue;
      onProgress?.((i / categorical.length) * 100, `DE stats: ${cluster.column}`);
      const stats = await computeDeStats(cluster, dataset.genes, matrix);

      if (!stats) continue;
      const bytes = encodeDeStatsBuffer(stats);

      files[cluster.column] = await this.compress(
        bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength) as ArrayBuffer,
      );
    }

    return files;
  }

  /**
   * Process palettes for categorical columns
   */
  async processPalettes(
    dataset: StandardizedDataset,
  ): Promise<Record<string, Blob>> {
    const files: Record<string, Blob> = {};

    if (!dataset.clusters || dataset.clusters.length === 0) {
      return files;
    }

    for (const cluster of dataset.clusters) {
      // Only save palettes for categorical columns that have a palette
      if (cluster.type === "categorical" && cluster.palette) {
        const json = JSON.stringify(cluster.palette, null, 2);
        const blob = new Blob([json], { type: "application/json" });

        files[cluster.column] = blob;
      }
    }

    return files;
  }

  /**
   * Process metadata for obs/metadata.json
   */
  async processMetadata(dataset: StandardizedDataset): Promise<Blob> {
    const metadata = {
      id: dataset.id,
      name: dataset.name,
      type: dataset.type,
      metadata: dataset.metadata,
      clusters: dataset.clusters,
    };

    const json = JSON.stringify(metadata, null, 2);

    return await this.compressText(json);
  }

  /**
   * Compress text data
   */
  private async compressText(text: string): Promise<Blob> {
    const blob = new Blob([text], { type: "application/json" });
    const stream = blob.stream().pipeThrough(new CompressionStream("gzip"));

    return new Response(stream).blob();
  }
}
